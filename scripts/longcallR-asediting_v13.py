#!/usr/bin/env python
r"""
longcallR-asediting.py

Allele-specific A-to-I RNA editing analysis using haplotype-tagged long-read
RNA-seq BAM files produced by longcallR.

Architecture
------------
Rather than opening the BAM thousands of times (once per site), this script:

  1. Loads the entire REDIportal database for BAM-covered chromosomes into
     memory (~few hundred MB for the full hg38 table).
  2. Builds repeat-element clusters from the in-memory site table (no BAM
     needed at this stage).
  3. Spawns one worker process per chromosome. Each worker does a single
     linear scan through that chromosome's reads and accumulates:
       - Per-site base counts (for the site-level analysis)
       - Per-read edited/unedited status per cluster (for cluster analysis)
  4. Merges counts across chromosomes, runs statistical tests, writes output.

This eliminates the per-site BAM seeking that made the previous architecture
slow at genome-wide scale.

Multiple testing correction strategy
-------------------------------------
Naive genome-wide BH correction across all sites is overly conservative because
RNA editing sites are strongly correlated — sites within the same Alu element
are edited by the same ADAR event acting on the same dsRNA structure and are
supported by the same set of reads. Correcting 100,000+ correlated tests as if
they were independent inflates the adjusted p-value floor far beyond what the
raw p-values warrant.

We instead use a two-stage hierarchical correction analogous to the approach
used for eQTL analysis (e.g. FastQTL):

  Stage 1 — within element/gene (reducing correlated tests to one per unit):
    - For sites in repeat elements (repeat_type != "NONREP"): the element-level
      p-value (P_element) is taken directly from the cluster-level test, which
      already aggregates reads across all sites in the element. This avoids
      re-applying Bonferroni on top of an already-aggregated test.
    - For sites NOT in repeat elements (NONREP sites): sites are grouped by
      gene_name and a Bonferroni correction is applied within each gene:
        P_element = min(P_site) * n_sites_in_gene
      This is conservative but computationally trivial and much better
      calibrated than genome-wide correction alone.
    - Sites with no gene annotation (gene_name == ".") are treated as their
      own group (one site = one element), so P_element = P_site.

  Stage 2 — across elements/genes (BH on the element-level p-values):
    BH correction is applied to the full set of element-level p-values,
    producing P_element_adj. This is the recommended value for genome-wide
    significance calling.

Output columns added by this procedure:
  P_element      — element/gene-level p-value from Stage 1
  P_element_adj  — BH-adjusted element p-value from Stage 2
  Sig_hierarchical — TRUE if BOTH P_element_adj < 0.05 (element-level FDR)
                    AND P_value < 0.05 (nominal site-level threshold).
                    A site is only flagged if it has evidence at both levels:
                    the element shows allele-specific editing overall, AND this
                    specific site shows nominal allele-specific signal.

The original naive P_value_adj column is retained for reference so that the
inflation from genome-wide correction can be directly observed.

REDIportal column layout (TABLE1 TSV, bgzipped + tabix-indexed with -S 1 -s 2 -b 3 -e 3):
    col0=Accession  col1=Region   col2=Position(1-based)  col3=Ref  col4=Ed
    col5=Strand     col6=db       col7=type               col8=dbsnp col9=repeat

Usage:
    python longcallR-asediting.py \
        -b phased.bam \
        -r TABLE1_hg38_v3.txt.gz \
        -o output_prefix \
        -t 32 \
        [--annotation annotation.gtf.gz] \
        [--min_coverage 10] \
        [--min_editing_rate 0.0] \
        [--overdispersion 0.001] \
        [--max_cluster_gap 300] \
        [--min_cluster_sites 2] \
        [--gene_types protein_coding lncRNA] \
        [--dna_vcf dna.vcf] \
        [--rna_vcf longcallR_phased.vcf]
"""

import argparse
import concurrent.futures
import gzip
import math
import os
import sys
from collections import defaultdict

import pysam
from scipy.stats import betabinom
from statsmodels.stats.multitest import multipletests


# ---------------------------------------------------------------------------
# Beta-binomial test
# ---------------------------------------------------------------------------

def convert_mu_rho_to_alpha_beta(mu, rho):
    phi = (1 - rho) / rho - 1
    alpha = mu * phi
    beta  = (1 - mu) * phi
    return alpha, beta


def beta_binomial_p_value(k_obs, n, mu, rho, alternative="two-sided"):
    """
    Beta-binomial test for allele-specific signal.
    mu is clamped to (1e-6, 1-1e-6) to prevent degenerate distributions.
    Returns 1.0 for degenerate inputs.
    """
    if n == 0:
        return 1.0
    mu = max(1e-6, min(1.0 - 1e-6, mu))
    alpha, beta_param = convert_mu_rho_to_alpha_beta(mu, rho)
    if alpha <= 0 or beta_param <= 0:
        return 1.0
    bb = betabinom(n, alpha, beta_param)
    p_obs = bb.pmf(k_obs)
    if alternative == "two-sided":
        p_value = sum(p for p in [bb.pmf(k) for k in range(n + 1)]
                      if p <= p_obs + 1e-15)
    elif alternative == "greater":
        p_value = float(bb.sf(k_obs - 1))
    elif alternative == "less":
        p_value = float(bb.cdf(k_obs))
    else:
        raise ValueError("Invalid alternative.")
    return float(p_value)


def site_p_value(h1_edit, h1_total, h2_edit, h2_total, overdispersion):
    """Compute p-value for a single site or cluster, handling invariant cases."""
    total = h1_total + h2_total
    if h1_total == 0 or h2_total == 0 or total == 0:
        return 1.0
    if h1_edit == 0 and h2_edit == 0:
        return 1.0
    if (h1_total - h1_edit) == 0 and (h2_total - h2_edit) == 0:
        return 1.0
    expected_rate = (h1_edit + h2_edit) / total
    return beta_binomial_p_value(
        h1_edit, h1_total, expected_rate, overdispersion, alternative="two-sided"
    )


# ---------------------------------------------------------------------------
# BAM helpers
# ---------------------------------------------------------------------------

def get_bam_chroms(bam_file):
    """Return set of chromosomes with mapped reads (uses index, instant)."""
    chroms = set()
    with pysam.AlignmentFile(bam_file, "rb") as bam:
        for stat in bam.get_index_statistics():
            if stat.mapped > 0:
                chroms.add(stat.contig)
    print(f"  BAM covers {len(chroms)} chromosome(s): {', '.join(sorted(chroms))}",
          flush=True)
    return chroms


# ---------------------------------------------------------------------------
# REDIportal loading (entire database into memory)
# ---------------------------------------------------------------------------

def load_rediportal(rediportal_file, bam_chroms):
    """
    Load the full REDIportal TABLE1 TSV (bgzipped) into memory, restricted to
    chromosomes present in the BAM.

    Returns:
        sites: dict (chrom, pos_0based) -> {
                   "strand": str,
                   "repeat_type": str,   # "ALU", "REP", "NONREP"
                   "repeat_name": str,   # e.g. "SINE/AluSx1"
               }
    """
    print(f"Loading REDIportal database from {rediportal_file}...", flush=True)
    sites = {}
    skipped = 0

    open_func = gzip.open if rediportal_file.endswith(".gz") else open
    with open_func(rediportal_file, "rt") as f:
        for line in f:
            if not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            # Skip header (Accession starts the header line)
            if parts[0] == "Accession":
                continue
            if len(parts) < 10:
                continue
            chrom = parts[1]
            if chrom not in bam_chroms:
                skipped += 1
                continue
            strand = parts[5]
            if strand not in ("+", "-"):
                continue
            pos_0based   = int(parts[2]) - 1   # convert 1-based -> 0-based
            repeat_type  = parts[7]
            repeat_name  = parts[9]
            sites[(chrom, pos_0based)] = {
                "strand":      strand,
                "repeat_type": repeat_type,
                "repeat_name": repeat_name,
            }

    n_rep = sum(1 for m in sites.values() if m["repeat_type"] != "NONREP")
    print(f"  Loaded {len(sites):,} sites ({skipped:,} skipped, not in BAM); "
          f"{n_rep:,} in repeat elements", flush=True)
    return sites


# ---------------------------------------------------------------------------
# Repeat-element cluster building (pure in-memory, no BAM)
# ---------------------------------------------------------------------------

def build_repeat_clusters(sites, max_gap, min_sites):
    """
    Group REDIportal sites into repeat-element clusters.

    Clustering rules:
      - Same chromosome, same strand, same repeat_name
      - repeat_type != "NONREP"
      - Consecutive sites within max_gap bp

    Returns list of cluster dicts:
        chrom, strand, rep_type, rep_name,
        site_positions (sorted list of 0-based ints)
    """
    # Collect repeat sites
    repeat_sites = []
    for (chrom, pos_0based), meta in sites.items():
        if meta["repeat_type"] == "NONREP":
            continue
        repeat_sites.append((
            chrom,
            meta["strand"],
            meta["repeat_name"],
            meta["repeat_type"],
            pos_0based,
        ))

    # Sort for greedy clustering: chrom, strand, repeat_name, position
    repeat_sites.sort(key=lambda x: (x[0], x[1], x[2], x[4]))

    clusters = []
    i = 0
    while i < len(repeat_sites):
        chrom, strand, rep_name, rep_type, pos = repeat_sites[i]
        cluster_pos = [pos]
        j = i + 1
        while j < len(repeat_sites):
            c2, s2, rn2, rt2, p2 = repeat_sites[j]
            if (c2 == chrom and s2 == strand and rn2 == rep_name
                    and (p2 - cluster_pos[-1]) <= max_gap):
                cluster_pos.append(p2)
                j += 1
            else:
                break
        if len(cluster_pos) >= min_sites:
            clusters.append({
                "chrom":    chrom,
                "strand":   strand,
                "rep_type": rep_type,
                "rep_name": rep_name,
                "sites":    cluster_pos,   # 0-based, sorted
            })
        i = j if j > i else i + 1

    return clusters


# ---------------------------------------------------------------------------
# Per-chromosome BAM scan worker
# ---------------------------------------------------------------------------

# ---------------------------------------------------------------------------
# Picklable factory functions for defaultdict accumulators
# ---------------------------------------------------------------------------
# ProcessPoolExecutor serialises worker arguments with pickle. Lambda functions
# defined inside other functions (nested lambdas) cannot be pickled because
# pickle needs a top-level name to reference. These module-level functions
# replace all nested lambdas used as defaultdict factories in _scan_chromosome.

def _site_count_factory():
    """Factory for per-site base count accumulator."""
    return {"ref": 0, "edit": 0}


def _cluster_hp_factory():
    """Factory for per-(phase_set, haplotype) cluster count accumulator."""
    return {"edited": 0, "unedited": 0, "total": 0}


def _cluster_count_factory():
    """Factory for per-cluster count accumulator (outer defaultdict)."""
    return defaultdict(_cluster_hp_factory)


def _scan_chromosome(args):
    """
    Worker function: scan all reads on one chromosome in a single linear pass.

    Accumulates:
      site_counts:    (chrom, pos_0based, ps, hp) -> {"ref": int, "edit": int}
      cluster_counts: cluster_id -> (ps, hp) -> {"edited": int, "unedited": int}

    A read contributes to a cluster count at most once per cluster, even if it
    covers multiple sites in that cluster (avoiding double-counting).
    """
    (bam_file, chrom, chrom_sites, chrom_clusters,
     min_baseq, min_mapq) = args

    # Build a sorted position list and site->index map for fast overlap queries
    # chrom_sites: dict pos_0based -> {"strand", "repeat_type", "repeat_name"}
    sorted_positions = sorted(chrom_sites.keys())
    pos_set = set(sorted_positions)

    # For each cluster, build a set of its site positions for O(1) lookup
    # cluster_id -> frozenset of 0-based positions
    cluster_site_sets = {
        cid: frozenset(c["sites"])
        for cid, c in enumerate(chrom_clusters)
    }
    # pos_0based -> list of cluster_ids that contain it
    pos_to_clusters = defaultdict(list)
    for cid, c in enumerate(chrom_clusters):
        for p in c["sites"]:
            pos_to_clusters[p].append(cid)

    # Accumulators
    # site_counts[(pos_0based, ps, hp)] = {"ref": 0, "edit": 0}
    # cluster_counts[cid][(ps, hp)] = {"edited": int, "unedited": int, "total": int}
    #
    # NOTE: we use module-level factory functions (_site_count_factory,
    # _cluster_count_factory, _cluster_hp_factory) instead of lambda
    # expressions here. ProcessPoolExecutor pickles worker arguments, and
    # nested lambdas (lambda inside lambda) cannot be pickled because pickle
    # requires a top-level importable name. Module-level functions are
    # picklable. See the factory function definitions above _scan_chromosome.
    site_counts    = defaultdict(_site_count_factory)
    cluster_counts = defaultdict(_cluster_count_factory)
    # Track which clusters each read overlaps (for total count), separate from
    # classification (which requires seeing ref/edit base at a site).
    # cluster_span[cid] = (fetch_start_0based, fetch_end_0based)
    cluster_spans = {
        cid: (min(c["sites"]), max(c["sites"]) + 1)
        for cid, c in enumerate(chrom_clusters)
    }

    with pysam.AlignmentFile(bam_file, "rb") as bam:
        for read in bam.fetch(chrom):
            if (read.is_unmapped or read.is_secondary
                    or read.is_supplementary or read.reference_end is None):
                continue
            if not read.has_tag("HP") or not read.has_tag("PS"):
                continue
            if read.mapping_quality < min_mapq:
                continue

            hp = read.get_tag("HP")
            ps = read.get_tag("PS")
            if hp not in (1, 2):
                continue

            seq = read.query_sequence
            quals = read.query_qualities
            if seq is None:
                continue

            # Build aligned_pairs dict: ref_pos -> query_pos
            # Only for positions within read span that overlap our site list
            read_start = read.reference_start
            read_end   = read.reference_end

            # Quick range check: any sites in [read_start, read_end)?
            # Use binary search boundaries to avoid iterating all positions
            import bisect
            lo = bisect.bisect_left(sorted_positions, read_start)
            hi = bisect.bisect_right(sorted_positions, read_end - 1)
            if lo >= hi:
                continue

            # Build aligned pairs only once per read
            aligned_pairs = {
                ref_pos: q_idx
                for q_idx, ref_pos in read.get_aligned_pairs(matches_only=True)
                if read_start <= ref_pos < read_end
            }

            # Track which clusters this read has already been classified for
            cluster_classified = {}  # cid -> "edited" or "unedited"

            for ref_pos in sorted_positions[lo:hi]:
                if ref_pos not in pos_set:
                    continue
                q_idx = aligned_pairs.get(ref_pos)
                if q_idx is None:
                    continue
                if quals is not None and quals[q_idx] < min_baseq:
                    continue

                meta   = chrom_sites[ref_pos]
                strand = meta["strand"]
                ref_base  = "A" if strand == "+" else "T"
                edit_base = "G" if strand == "+" else "C"
                base = seq[q_idx].upper()

                # --- Per-site accumulation ---
                if base == ref_base:
                    site_counts[(ref_pos, ps, hp)]["ref"] += 1
                elif base == edit_base:
                    site_counts[(ref_pos, ps, hp)]["edit"] += 1

                # --- Per-cluster accumulation (one count per read per cluster) ---
                for cid in pos_to_clusters.get(ref_pos, []):
                    if cid in cluster_classified:
                        # Already classified for this cluster; upgrade to edited
                        # if we see an edit (edited trumps unedited)
                        if base == edit_base:
                            cluster_classified[cid] = "edited"
                    else:
                        if base == edit_base:
                            cluster_classified[cid] = "edited"
                        elif base == ref_base:
                            cluster_classified[cid] = "unedited"

            # Commit cluster classifications (edited/unedited)
            for cid, status in cluster_classified.items():
                cluster_counts[cid][(ps, hp)][status] += 1

            # Count total overlapping reads per cluster (the true denominator).
            # This includes reads that covered the cluster interval but showed
            # neither the ref nor edit base at any site (e.g. due to mismatches,
            # base quality failures, or reads not reaching any specific site).
            # We use read_start/read_end for the overlap check since aligned_pairs
            # is already built at this point.
            for cid, (cspan_start, cspan_end) in cluster_spans.items():
                if read_start < cspan_end and read_end > cspan_start:
                    cluster_counts[cid][(ps, hp)]["total"] += 1

    # Convert all defaultdicts to plain dicts before returning to the main
    # process. ProcessPoolExecutor pickles return values, and defaultdict
    # objects with lambda/function factories can fail to unpickle correctly
    # depending on Python version and how the factory was defined. Plain dicts
    # are always safely picklable. We must convert all nested levels.
    plain_site_counts = dict(site_counts)
    plain_cluster_counts = {
        cid: {ps_hp: dict(counts) for ps_hp, counts in hp_dict.items()}
        for cid, hp_dict in cluster_counts.items()
    }
    return chrom, plain_site_counts, plain_cluster_counts


# ---------------------------------------------------------------------------
# Count aggregation and statistical testing
# ---------------------------------------------------------------------------

def aggregate_site_counts(all_site_counts, sites):
    """
    Merge per-chromosome site count dicts into per-site haplotype counts.

    Returns list of site result dicts (before p-value computation).
    """
    # merged[(chrom, pos_0based)][(ps, hp)] = {"ref": int, "edit": int}
    merged = defaultdict(lambda: defaultdict(lambda: {"ref": 0, "edit": 0}))
    for chrom, site_counts in all_site_counts:
        for (pos_0based, ps, hp), counts in site_counts.items():
            merged[(chrom, pos_0based)][(ps, hp)]["ref"]  += counts["ref"]
            merged[(chrom, pos_0based)][(ps, hp)]["edit"] += counts["edit"]
    return merged


def aggregate_cluster_counts(all_cluster_counts, chrom_cluster_offsets):
    """
    Merge per-chromosome cluster count dicts, adjusting for per-chromosome
    cluster ID offsets.

    Returns dict: global_cluster_id -> (ps, hp) ->
        {"edited": int, "unedited": int, "total": int}

    "total" is the count of all HP-tagged reads overlapping the cluster
    interval, used as the coverage filter denominator. "edited" and "unedited"
    are subsets of reads that additionally showed the edit/ref base at a site.
    """
    # Use module-level factory functions to avoid pickle errors (same reason
    # as in _scan_chromosome — nested lambdas are not picklable).
    merged = defaultdict(_cluster_count_factory)
    for chrom, cluster_counts, offset in all_cluster_counts:
        for local_cid, ps_hp_counts in cluster_counts.items():
            global_cid = local_cid + offset
            for (ps, hp), counts in ps_hp_counts.items():
                merged[global_cid][(ps, hp)]["edited"]   += counts.get("edited", 0)
                merged[global_cid][(ps, hp)]["unedited"] += counts.get("unedited", 0)
                merged[global_cid][(ps, hp)]["total"]    += counts.get("total", 0)
    return merged


def pick_best_phase_set(ps_hp_counts, count_keys):
    """
    Given a dict of (ps, hp) -> counts, pick the phase set with the most
    total reads (summing across both haplotypes).

    count_keys: tuple of count dict keys to sum, e.g. ("ref", "edit") or
                ("total",) — supports single-key tuples for the cluster total
                denominator, or multi-key tuples for the edited+unedited sum.
    """
    ps_totals = defaultdict(int)
    for (ps, hp), counts in ps_hp_counts.items():
        ps_totals[ps] += sum(counts.get(k, 0) for k in count_keys)
    if not ps_totals:
        return None
    return max(ps_totals, key=ps_totals.get)


# ---------------------------------------------------------------------------
# Optional: gene annotation
# ---------------------------------------------------------------------------

def load_gene_annotation(annotation_file, gene_types, bam_chroms):
    if annotation_file is None:
        return None

    print(f"Loading gene annotation from {annotation_file}...", flush=True)
    from intervaltree import Interval, IntervalTree
    trees = defaultdict(IntervalTree)
    n_genes = 0

    open_func = gzip.open if annotation_file.endswith(".gz") else open
    file_type = "gff3" if ".gff3" in annotation_file else "gtf"

    with open_func(annotation_file, "rt") as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.strip().split("\t")
            if len(parts) < 9 or parts[2] != "gene":
                continue
            chrom = parts[0]
            if chrom not in bam_chroms:
                continue
            attr_dict = {}
            if file_type == "gtf":
                for attr in parts[8].strip().split(";"):
                    attr = attr.strip()
                    if not attr:
                        continue
                    kv = attr.split(" ", 1)
                    if len(kv) == 2:
                        attr_dict[kv[0]] = kv[1].replace('"', '')
            else:
                for attr in parts[8].strip().split(";"):
                    attr = attr.strip()
                    if "=" not in attr:
                        continue
                    k, v = attr.split("=", 1)
                    attr_dict[k] = v.replace('"', '')
            gene_type = attr_dict.get("gene_type", attr_dict.get("gene_biotype", ""))
            if gene_type not in gene_types:
                continue
            gene_name = attr_dict.get("gene_name", attr_dict.get("gene_id", "."))
            start = int(parts[3])
            end   = int(parts[4])
            trees[chrom].add(Interval(start, end + 1, gene_name))
            n_genes += 1

    print(f"  Loaded {n_genes:,} gene record(s) on BAM chromosome(s)", flush=True)
    return trees


def get_gene_name(trees, chrom, pos_0based):
    if trees is None:
        return "."
    hits = trees[chrom].overlap(pos_0based + 1, pos_0based + 2)
    if not hits:
        return "."
    return min(hits, key=lambda iv: iv.end - iv.begin).data


# ---------------------------------------------------------------------------
# Optional: DNA VCF filter
# ---------------------------------------------------------------------------

def load_dna_vcf(vcf_file):
    dna_vcfs = {}
    with pysam.VariantFile(vcf_file) as vcf:
        for record in vcf.fetch():
            gt = record.samples[0]["GT"]
            if len(record.ref) != 1 or any(len(a) != 1 for a in record.alts):
                continue
            if gt in [(0, 1), (1, 0)]:
                dna_vcfs[f"{record.contig}:{record.pos}"] = True
    return dna_vcfs


def load_longcallR_phased_vcf(vcf_file):
    rna_vcfs = defaultdict(list)
    with pysam.VariantFile(vcf_file) as vcf:
        for record in vcf.fetch():
            if "PASS" not in record.filter.keys():
                continue
            gt = record.samples[0]["GT"]
            if len(record.ref) != 1 or any(len(a) != 1 for a in record.alts):
                continue
            if gt in [(0, 1), (1, 0)] and record.samples[0].phased:
                ps = record.samples[0].get("PS", None)
                if ps and ps != ".":
                    rna_vcfs[ps].append(f"{record.contig}:{record.pos}")
    return rna_vcfs


def phase_set_has_dna_anchor(ps, rna_vcfs, dna_vcfs):
    """Return True if phase set contains at least one DNA-confirmed SNP."""
    return any(
        f"{s.split(':')[0]}:{s.split(':')[1]}" in dna_vcfs
        for s in rna_vcfs.get(ps, [])
    )


# ---------------------------------------------------------------------------
# Output helpers
# ---------------------------------------------------------------------------

SITE_HEADER = (
    "#Chrom\tPos\tStrand\tGene_name\tRepeat_type\tPhase_set\t"
    "H1_unedited\tH1_edited\tH2_unedited\tH2_edited\t"
    "H1_editing_rate\tH2_editing_rate\tOverall_editing_rate\t"
    "LogFC_H1_vs_H2\tP_value\tP_value_adj\tSignificant\t"
    "P_element\tP_element_adj\tSig_hierarchical\n"
)

CLUSTER_HEADER = (
    "#Chrom\tCluster_start\tCluster_end\tStrand\tRepeat_type\tRepeat_name\t"
    "N_sites\tSites\tPhase_set\t"
    "H1_unedited\tH1_edited\tH2_unedited\tH2_edited\t"
    "H1_editing_rate\tH2_editing_rate\tOverall_editing_rate\t"
    "LogFC_H1_vs_H2\tP_value\tP_value_adj\tSignificant\t"
    "P_element\tP_element_adj\tSig_hierarchical\n"
)


def fmt_rate(r):
    return f"{r:.4f}" if not math.isnan(r) else "NA"


def compute_rates_and_logfc(h1_edit, h1_total, h2_edit, h2_total):
    h1_rate = h1_edit / h1_total if h1_total > 0 else float("nan")
    h2_rate = h2_edit / h2_total if h2_total > 0 else float("nan")
    logfc   = (math.log2((h1_edit + 1) / (h1_total + 1))
               - math.log2((h2_edit + 1) / (h2_total + 1)))
    overall = ((h1_edit + h2_edit) / (h1_total + h2_total)
               if (h1_total + h2_total) > 0 else 0.0)
    return h1_rate, h2_rate, overall, logfc


# ---------------------------------------------------------------------------
# Two-stage hierarchical multiple testing correction
# ---------------------------------------------------------------------------

def hierarchical_correction(site_results, cluster_results):
    """
    Apply two-stage hierarchical multiple testing correction to site-level
    results, analogous to the approach used in eQTL analysis.

    Stage 1 — reduce correlated tests to one p-value per biological unit:

      For sites in repeat elements (repeat_type != "NONREP"):
        P_element is taken directly from the corresponding cluster's p-value.
        Rationale: the cluster test already aggregates reads across all sites in
        the element using a per-read (not per-observation) model, so it is the
        correct element-level summary statistic. Applying Bonferroni on top of
        it would double-correct.

      For NONREP sites:
        Sites are grouped by gene_name. Within each gene, Bonferroni correction
        is applied: P_element = min(P_site_in_gene) * n_sites_in_gene.
        Sites with no gene annotation (gene_name == ".") are treated
        individually: P_element = P_site.
        Rationale: NONREP sites lack the shared-read structure of repeat
        elements, so within-gene Bonferroni is the appropriate conservative
        correction for the number of independent-ish tests per gene.

    Stage 2 — BH correction across all element/gene-level p-values:
        P_element_adj is computed from BH applied to all P_element values.
        This controls FDR at the level of biological units (elements/genes)
        rather than individual sites, which is the appropriate level given
        the within-element correlation structure.

    Parameters
    ----------
    site_results : list of dicts
        Per-site results, each must have keys: p_value, gene_name, repeat_type,
        pos_1based, chrom. Must already have naive p_value_adj populated.
    cluster_results : list of dicts
        Per-cluster results, each must have keys: chrom, start, end, p_value.

    Returns
    -------
    site_results : list of dicts, with P_element and P_element_adj added in-place
    """
    # Build a lookup from cluster (chrom, start_1based, end_1based) -> p_value
    # so we can join clusters back to their member sites.
    cluster_p = {}
    for cr in cluster_results:
        key = (cr["chrom"], cr["start"], cr["end"])
        cluster_p[key] = cr["p_value"]

    # For each repeat site, find its parent cluster by checking which cluster
    # interval contains its position. We build a simple interval lookup per chrom.
    # cluster_intervals[chrom] = list of (start_1based, end_1based, p_value)
    cluster_intervals = defaultdict(list)
    for cr in cluster_results:
        cluster_intervals[cr["chrom"]].append(
            (cr["start"], cr["end"], cr["p_value"])
        )

    def find_cluster_p(chrom, pos_1based):
        """Return the cluster p-value for a site, or None if not in any cluster."""
        for cstart, cend, cp in cluster_intervals.get(chrom, []):
            if cstart <= pos_1based <= cend:
                return cp
        return None

    # Stage 1: assign P_element to each site
    # For NONREP sites: group by gene, compute Bonferroni within gene
    nonrep_by_gene = defaultdict(list)   # gene -> list of site indices
    for i, r in enumerate(site_results):
        if r["repeat_type"] == "NONREP":
            nonrep_by_gene[r["gene_name"]].append(i)
        else:
            # Repeat site: P_element = parent cluster p-value
            cp = find_cluster_p(r["chrom"], r["pos_1based"])
            if cp is not None:
                r["p_element"] = cp
            else:
                # Site in a repeat element but no cluster passed filters
                # (e.g. cluster failed coverage threshold). Fall back to
                # within-element Bonferroni — but we only have one site here
                # so P_element = P_site.
                r["p_element"] = r["p_value"]

    # Bonferroni within each gene for NONREP sites
    for gene, indices in nonrep_by_gene.items():
        n_in_gene = len(indices)
        min_p = min(site_results[i]["p_value"] for i in indices)
        p_element = min(min_p * n_in_gene, 1.0)   # cap at 1.0
        for i in indices:
            site_results[i]["p_element"] = p_element

    # Stage 2: BH across all element-level p-values
    all_p_element = [r["p_element"] for r in site_results]
    if all_p_element:
        reject, adj_p, _, _ = multipletests(all_p_element, alpha=0.05, method="fdr_bh")
        for i, r in enumerate(site_results):
            r["p_element_adj"] = adj_p[i]
            # Sig_hierarchical requires BOTH conditions:
            #   1. The element-level test passes FDR < 0.05 (P_element_adj < 0.05)
            #   2. The site's own p-value passes a nominal threshold (P_value < 0.05)
            #
            # Rationale: P_element for a repeat site is the cluster p-value, which
            # reflects the aggregate allele-specific editing of the whole element.
            # A site with P_value=1.0 (no site-level signal) should not be marked
            # significant just because its parent cluster is significant — the cluster
            # result is already reported in the cluster output file, and the site-level
            # output should only flag sites with genuine site-level evidence.
            # The nominal threshold of 0.05 is intentionally uncorrected since
            # the element-level BH correction has already controlled FDR globally.
            r["sig_hierarchical"] = reject[i] and r["p_value"] < 0.05
    else:
        for r in site_results:
            r["p_element_adj"] = 1.0
            r["sig_hierarchical"] = False

    return site_results


def hierarchical_correction_clusters(cluster_results):
    """
    Apply BH correction across cluster p-values only (Stage 2 for clusters).

    For clusters, Stage 1 is trivially already done — each cluster IS the
    element-level unit, so P_element = P_value. Stage 2 is BH across clusters,
    which is what P_value_adj already computes. We add P_element_adj as an
    explicit alias to make the output consistent with the site-level file and
    to document that clusters undergo the same Stage 2 correction.

    The P_element column is set equal to P_value (the cluster test p-value)
    to be explicit that no within-element aggregation step was needed.
    """
    for r in cluster_results:
        # P_element = P_value for clusters: the cluster test IS the
        # element-level summary, no further within-element reduction needed.
        r["p_element"] = r["p_value"]
        # P_element_adj = P_value_adj: BH across clusters is Stage 2.
        r["p_element_adj"] = r.get("p_value_adj", r["p_value"])
        r["sig_hierarchical"] = r.get("significant", False)
    return cluster_results


# ---------------------------------------------------------------------------
# Main analysis
# ---------------------------------------------------------------------------

def analyze(bam_file, sites, clusters, annotation_trees, output_prefix,
            min_coverage, min_cluster_coverage, min_editing_rate, overdispersion,
            threads, min_baseq, min_mapq,
            dna_vcfs, rna_vcfs):

    bam_chroms = set(sites_chrom for sites_chrom, _ in
                     [(k[0], k) for k in sites.keys()])

    # Group sites and clusters by chromosome for per-chromosome workers
    sites_by_chrom = defaultdict(dict)
    for (chrom, pos_0based), meta in sites.items():
        sites_by_chrom[chrom][pos_0based] = meta

    # Assign global cluster IDs and split clusters by chromosome
    # Each cluster belongs to exactly one chromosome
    clusters_by_chrom = defaultdict(list)
    cluster_offsets = {}   # chrom -> global ID of first cluster on that chrom
    global_id = 0
    for chrom in sorted(sites_by_chrom.keys()):
        chrom_clusters = [c for c in clusters if c["chrom"] == chrom]
        cluster_offsets[chrom] = global_id
        clusters_by_chrom[chrom] = chrom_clusters
        global_id += len(chrom_clusters)
    total_clusters = global_id

    chroms_to_scan = sorted(sites_by_chrom.keys())
    print(f"Scanning {len(chroms_to_scan)} chromosome(s) with {threads} worker(s)...",
          flush=True)

    worker_args = [
        (bam_file, chrom,
         sites_by_chrom[chrom],
         clusters_by_chrom.get(chrom, []),
         min_baseq, min_mapq)
        for chrom in chroms_to_scan
    ]

    all_site_counts    = []   # list of (chrom, site_counts_dict)
    all_cluster_counts = []   # list of (chrom, cluster_counts_dict, offset)

    with concurrent.futures.ProcessPoolExecutor(max_workers=threads) as executor:
        futures = {executor.submit(_scan_chromosome, a): a[1]
                   for a in worker_args}
        n_done = 0
        for future in concurrent.futures.as_completed(futures):
            chrom = futures[future]
            try:
                chrom_out, site_counts, cluster_counts = future.result()
                all_site_counts.append((chrom_out, site_counts))
                all_cluster_counts.append(
                    (chrom_out, cluster_counts, cluster_offsets.get(chrom_out, 0))
                )
            except Exception as e:
                print(f"  Warning: error on {chrom}: {e}", file=sys.stderr, flush=True)
            n_done += 1
            print(f"  {n_done}/{len(chroms_to_scan)} chromosomes done...", flush=True)

    # ---- Per-site results ----
    print("Computing per-site statistics...", flush=True)
    merged_sites = aggregate_site_counts(all_site_counts, sites)

    site_results = []
    for (chrom, pos_0based), ps_hp_counts in merged_sites.items():
        meta   = sites[(chrom, pos_0based)]
        strand = meta["strand"]

        best_ps = pick_best_phase_set(ps_hp_counts, ("ref", "edit"))
        if best_ps is None:
            continue

        # Filter to chosen phase set
        h1 = ps_hp_counts.get((best_ps, 1), {"ref": 0, "edit": 0})
        h2 = ps_hp_counts.get((best_ps, 2), {"ref": 0, "edit": 0})
        h1_ref, h1_edit = h1["ref"], h1["edit"]
        h2_ref, h2_edit = h2["ref"], h2["edit"]
        h1_total = h1_ref + h1_edit
        h2_total = h2_ref + h2_edit
        total    = h1_total + h2_total

        if total < min_coverage:
            continue
        overall_rate = (h1_edit + h2_edit) / total
        if overall_rate < min_editing_rate:
            continue
        if (dna_vcfs is not None
                and not phase_set_has_dna_anchor(best_ps, rna_vcfs, dna_vcfs)):
            continue

        p_val = site_p_value(h1_edit, h1_total, h2_edit, h2_total, overdispersion)
        h1_rate, h2_rate, overall, logfc = compute_rates_and_logfc(
            h1_edit, h1_total, h2_edit, h2_total
        )
        gene_name = get_gene_name(annotation_trees, chrom, pos_0based)

        site_results.append({
            "chrom":        chrom,
            "pos_1based":   pos_0based + 1,
            "strand":       strand,
            "gene_name":    gene_name,
            "phase_set":    best_ps,
            "h1_ref":       h1_ref,
            "h1_edit":      h1_edit,
            "h2_ref":       h2_ref,
            "h2_edit":      h2_edit,
            "h1_rate":      h1_rate,
            "h2_rate":      h2_rate,
            "overall_rate": overall,
            "logfc":        logfc,
            "p_value":      p_val,
            "repeat_type":  meta["repeat_type"],
        })

    print(f"  {len(site_results):,} site(s) passing filters", flush=True)

    # -----------------------------------------------------------------------
    # Cluster results are computed BEFORE site writing because the hierarchical
    # correction for sites uses cluster p-values to annotate repeat-element
    # sites (P_element = cluster p-value for repeat sites). Computing clusters
    # first ensures cluster_results is available when hierarchical_correction
    # is called during site writing.
    # -----------------------------------------------------------------------

    # ---- Cluster results ----
    print("Computing per-cluster statistics...", flush=True)
    merged_clusters = aggregate_cluster_counts(all_cluster_counts, cluster_offsets)

    # Diagnostic: report distribution of cluster total read counts so the user
    # can see where coverage is falling relative to min_cluster_coverage.
    # We do this before filtering so the full distribution is visible.
    # "Total reads" here means all HP-tagged reads overlapping the cluster
    # interval, which is the correct denominator for the coverage filter (see
    # note in _scan_chromosome about the distinction between total and
    # edited/unedited counts).
    all_totals = []
    for global_cid, ps_hp_counts in merged_clusters.items():
        # Sum total reads across both haplotypes for the best phase set
        best_ps_diag = pick_best_phase_set(ps_hp_counts, ("total",))
        if best_ps_diag is None:
            all_totals.append(0)
            continue
        h1t = ps_hp_counts.get((best_ps_diag, 1), {"total": 0})["total"]
        h2t = ps_hp_counts.get((best_ps_diag, 2), {"total": 0})["total"]
        all_totals.append(h1t + h2t)

    if all_totals:
        all_totals_sorted = sorted(all_totals)
        n = len(all_totals_sorted)
        def pct(p):
            return all_totals_sorted[min(int(n * p / 100), n - 1)]
        print(f"  Cluster total read count distribution (n={n:,} clusters with any reads):",
              flush=True)
        print(f"    min={all_totals_sorted[0]}  p10={pct(10)}  p25={pct(25)}  "
              f"median={pct(50)}  p75={pct(75)}  p90={pct(90)}  max={all_totals_sorted[-1]}",
              flush=True)
        n_passing = sum(1 for t in all_totals if t >= min_cluster_coverage)
        print(f"  Clusters with total reads >= {min_cluster_coverage} "
              f"(min_cluster_coverage): {n_passing:,}/{n:,}", flush=True)

    cluster_results = []
    for global_cid, ps_hp_counts in merged_clusters.items():
        c = clusters[global_cid]

        # Use total overlapping reads (not just edited+unedited) as the
        # denominator for the coverage filter. This avoids rejecting clusters
        # where many reads cover the interval but happen to show neither the
        # ref nor edit base (e.g. due to SNPs, base quality failures, or
        # reads spanning the cluster but not reaching any specific site).
        # The edited/unedited counts are still used for the statistical test.
        best_ps = pick_best_phase_set(ps_hp_counts, ("total",))
        if best_ps is None:
            continue

        h1 = ps_hp_counts.get((best_ps, 1), {"edited": 0, "unedited": 0, "total": 0})
        h2 = ps_hp_counts.get((best_ps, 2), {"edited": 0, "unedited": 0, "total": 0})
        h1_edit, h1_unedited, h1_total = h1["edited"], h1["unedited"], h1["total"]
        h2_edit, h2_unedited, h2_total = h2["edited"], h2["unedited"], h2["total"]
        total = h1_total + h2_total

        # Coverage filter uses total overlapping reads as the denominator.
        if total < min_cluster_coverage:
            continue

        # Editing rate filter uses edited/(edited+unedited), i.e. reads that
        # actually voted at a site. This is correct: reads that overlapped the
        # cluster but didn't cover any site don't inform the editing rate.
        classifiable = (h1_edit + h1_unedited) + (h2_edit + h2_unedited)
        if classifiable == 0:
            continue
        overall_rate = (h1_edit + h2_edit) / classifiable
        if overall_rate < min_editing_rate:
            continue

        if (dna_vcfs is not None
                and not phase_set_has_dna_anchor(best_ps, rna_vcfs, dna_vcfs)):
            continue

        # Statistical test uses edited vs unedited (classifiable reads only).
        # h1_total/h2_total used here are the classifiable counts, not the
        # total overlapping counts, because only classifiable reads contribute
        # evidence about the editing rate.
        h1_classifiable = h1_edit + h1_unedited
        h2_classifiable = h2_edit + h2_unedited
        p_val = site_p_value(
            h1_edit, h1_classifiable, h2_edit, h2_classifiable, overdispersion
        )
        h1_rate, h2_rate, overall, logfc = compute_rates_and_logfc(
            h1_edit, h1_classifiable, h2_edit, h2_classifiable
        )

        cluster_results.append({
            "chrom":         c["chrom"],
            "start":         min(c["sites"]) + 1,   # 1-based
            "end":           max(c["sites"]) + 1,
            "strand":        c["strand"],
            "rep_type":      c["rep_type"],
            "rep_name":      c["rep_name"],
            "n_sites":       len(c["sites"]),
            "sites":         ",".join(str(p + 1) for p in c["sites"]),  # 1-based
            "phase_set":     best_ps,
            "h1_unedited":   h1_unedited,
            "h1_edit":       h1_edit,
            "h2_unedited":   h2_unedited,
            "h2_edit":       h2_edit,
            "h1_total":      h1_total,   # total overlapping reads (incl. unclassifiable)
            "h2_total":      h2_total,
            "h1_rate":       h1_rate,
            "h2_rate":       h2_rate,
            "overall_rate":  overall,
            "logfc":         logfc,
            "p_value":       p_val,
        })

    print(f"  {len(cluster_results):,} cluster(s) passing filters", flush=True)

    # ---- Per-site output ----
    # Written here (after cluster_results is computed) so that
    # hierarchical_correction can use cluster p-values for repeat-element sites.
    site_out = output_prefix + ".asediting.tsv"
    print(f"Writing per-site results to {site_out}...", flush=True)
    with open(site_out, "w") as f:
        f.write(SITE_HEADER)
        if site_results:
            # Naive BH across all sites (retained for reference).
            p_values = [r["p_value"] for r in site_results]
            reject, adj_p, _, _ = multipletests(p_values, alpha=0.05, method="fdr_bh")
            for i, r in enumerate(site_results):
                r["p_value_adj"] = adj_p[i]
                r["significant"] = reject[i]

            # Two-stage hierarchical correction (see module docstring).
            # cluster_results is now populated; repeat sites get their cluster
            # p-value as P_element. NONREP sites get Bonferroni-within-gene.
            site_results = hierarchical_correction(site_results, cluster_results)

            site_results.sort(key=lambda r: r["p_element_adj"])
            for r in site_results:
                f.write(
                    f"{r['chrom']}\t{r['pos_1based']}\t{r['strand']}\t"
                    f"{r['gene_name']}\t{r['repeat_type']}\t{r['phase_set']}\t"
                    f"{r['h1_ref']}\t{r['h1_edit']}\t{r['h2_ref']}\t{r['h2_edit']}\t"
                    f"{fmt_rate(r['h1_rate'])}\t{fmt_rate(r['h2_rate'])}\t"
                    f"{r['overall_rate']:.4f}\t{r['logfc']:.4f}\t"
                    f"{r['p_value']:.6e}\t{r['p_value_adj']:.6e}\t"
                    f"{'TRUE' if r['significant'] else 'FALSE'}\t"
                    f"{r['p_element']:.6e}\t{r['p_element_adj']:.6e}\t"
                    f"{'TRUE' if r['sig_hierarchical'] else 'FALSE'}\n"
                )
    n_sig_naive = sum(1 for r in site_results if r.get("significant", False))
    n_sig_hier  = sum(1 for r in site_results if r.get("sig_hierarchical", False))
    print(f"  Significant sites — naive BH: {n_sig_naive} | "
          f"hierarchical: {n_sig_hier} (of {len(site_results)} total)", flush=True)

    # Write cluster output (always, even if empty)
    cluster_out = output_prefix + ".asediting_clusters.tsv"
    print(f"Writing cluster results to {cluster_out}...", flush=True)
    with open(cluster_out, "w") as f:
        f.write(CLUSTER_HEADER)
        if cluster_results:
            # Stage 2 naive correction: BH across all clusters.
            # For clusters, Stage 1 is trivially complete — each cluster IS the
            # element-level unit (P_element = P_value). Stage 2 BH across
            # clusters is therefore equivalent to the hierarchical correction,
            # so P_element_adj = P_value_adj for clusters.
            p_values = [r["p_value"] for r in cluster_results]
            reject, adj_p, _, _ = multipletests(p_values, alpha=0.05, method="fdr_bh")
            for i, r in enumerate(cluster_results):
                r["p_value_adj"] = adj_p[i]
                r["significant"] = reject[i]

            # Add explicit P_element / P_element_adj / Sig_hierarchical columns
            # for consistency with the site-level output. See
            # hierarchical_correction_clusters() docstring for rationale.
            cluster_results = hierarchical_correction_clusters(cluster_results)

            cluster_results.sort(key=lambda r: r["p_element_adj"])
            for r in cluster_results:
                f.write(
                    f"{r['chrom']}\t{r['start']}\t{r['end']}\t{r['strand']}\t"
                    f"{r['rep_type']}\t{r['rep_name']}\t{r['n_sites']}\t"
                    f"{r['sites']}\t{r['phase_set']}\t"
                    f"{r['h1_unedited']}\t{r['h1_edit']}\t"
                    f"{r['h2_unedited']}\t{r['h2_edit']}\t"
                    f"{fmt_rate(r['h1_rate'])}\t{fmt_rate(r['h2_rate'])}\t"
                    f"{r['overall_rate']:.4f}\t{r['logfc']:.4f}\t"
                    f"{r['p_value']:.6e}\t{r['p_value_adj']:.6e}\t"
                    f"{'TRUE' if r['significant'] else 'FALSE'}\t"
                    f"{r['p_element']:.6e}\t{r['p_element_adj']:.6e}\t"
                    f"{'TRUE' if r['sig_hierarchical'] else 'FALSE'}\n"
                )
    n_sig = sum(1 for r in cluster_results if r.get("significant", False))
    print(f"  Significant clusters (FDR<0.05): {n_sig}/{len(cluster_results)}", flush=True)


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Allele-specific A-to-I RNA editing analysis using longcallR "
                    "haplotype-tagged BAM files. Scans the BAM once per chromosome "
                    "for genome-wide scalability."
    )
    parser.add_argument("-b", "--bam", required=True,
                        help="Haplotype-tagged BAM file from longcallR (must be indexed)")
    parser.add_argument("-r", "--rediportal", required=True,
                        help="REDIportal TABLE1 TSV file (plain or bgzipped). "
                             "Column layout: Accession Region Position(1-based) Ref Ed "
                             "Strand db type dbsnp repeat")
    parser.add_argument("-o", "--output_prefix", required=True,
                        help="Prefix for output files")
    parser.add_argument("-a", "--annotation", default=None,
                        help="Annotation GTF/GFF3 for gene name lookup (optional)")
    parser.add_argument("--gene_types", nargs="+", default=["protein_coding", "lncRNA"],
                        help="Gene biotypes to include. Default: protein_coding lncRNA")
    parser.add_argument("-t", "--threads", type=int, default=1,
                        help="Worker processes (one per chromosome). Default: 1")
    parser.add_argument("--min_coverage", type=int, default=10,
                        help="Minimum total phased reads (H1+H2) at a site/cluster. "
                             "Default: 10")
    parser.add_argument("--min_editing_rate", type=float, default=0.0,
                        help="Minimum overall editing rate to retain a site/cluster. "
                             "Default: 0.0")
    parser.add_argument("--min_baseq", type=int, default=20,
                        help="Minimum base quality. Default: 20")
    parser.add_argument("--min_mapq", type=int, default=20,
                        help="Minimum mapping quality. Default: 20")
    parser.add_argument("--overdispersion", type=float, default=0.001,
                        help="Beta-binomial overdispersion parameter rho. Default: 0.001")
    parser.add_argument("--max_cluster_gap", type=int, default=300,
                        help="Maximum gap (bp) between sites in the same repeat-element "
                             "cluster. Default: 300")
    parser.add_argument("--min_cluster_sites", type=int, default=2,
                        help="Minimum sites per cluster to report. Default: 2")
    parser.add_argument("--min_cluster_coverage", type=int, default=5,
                        help="Minimum total HP-tagged reads overlapping a cluster "
                             "interval (H1+H2) to retain it for testing. This uses "
                             "the total overlapping read count as the denominator, "
                             "not just reads classifiable as edited/unedited, to avoid "
                             "rejecting clusters where reads cover the element but don't "
                             "reach a specific REDIportal site. Lower than --min_coverage "
                             "by default because cluster reads are a union not a sum. "
                             "Default: 5")
    parser.add_argument("--dna_vcf", default=None,
                        help="DNA VCF. With --rna_vcf, restricts to phase blocks "
                             "anchored by at least one genomic SNP.")
    parser.add_argument("--rna_vcf", default=None,
                        help="LongcallR phased RNA VCF. Must be paired with --dna_vcf.")

    args = parser.parse_args()

    if (args.dna_vcf is None) ^ (args.rna_vcf is None):
        parser.error("--dna_vcf and --rna_vcf must be provided together or not at all.")

    # 1. Get BAM chromosomes
    print("Reading BAM index...", flush=True)
    bam_chroms = get_bam_chroms(args.bam)

    # 2. Load full REDIportal database into memory (filtered to BAM chromosomes)
    sites = load_rediportal(args.rediportal, bam_chroms)
    if not sites:
        print("Error: no REDIportal sites found for BAM chromosomes.\n"
              "Check chromosome naming (chr1 vs 1).", file=sys.stderr)
        sys.exit(1)

    # 3. Build repeat-element clusters (pure in-memory, instant)
    print("Building repeat-element clusters...", flush=True)
    clusters = build_repeat_clusters(
        sites, args.max_cluster_gap, args.min_cluster_sites
    )
    print(f"  {len(clusters):,} cluster(s) with >= {args.min_cluster_sites} site(s)",
          flush=True)

    # 4. Load gene annotation (optional)
    annotation_trees = load_gene_annotation(
        args.annotation, set(args.gene_types), bam_chroms
    )

    # 5. Load optional VCF filters
    dna_vcfs = rna_vcfs = None
    if args.dna_vcf and args.rna_vcf:
        print("Loading VCF files for phase-block filtering...", flush=True)
        dna_vcfs = load_dna_vcf(args.dna_vcf)
        rna_vcfs = load_longcallR_phased_vcf(args.rna_vcf)
        print(f"  DNA VCF: {len(dna_vcfs):,} SNPs | RNA VCF: {len(rna_vcfs):,} phase sets",
              flush=True)

    # 6. Run analysis
    analyze(
        bam_file=args.bam,
        sites=sites,
        clusters=clusters,
        annotation_trees=annotation_trees,
        output_prefix=args.output_prefix,
        min_coverage=args.min_coverage,
        min_cluster_coverage=args.min_cluster_coverage,
        min_editing_rate=args.min_editing_rate,
        overdispersion=args.overdispersion,
        threads=args.threads,
        min_baseq=args.min_baseq,
        min_mapq=args.min_mapq,
        dna_vcfs=dna_vcfs,
        rna_vcfs=rna_vcfs,
    )
