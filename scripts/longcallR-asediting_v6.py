#!/usr/bin/env python
r"""
longcallR-asediting.py

Allele-specific A-to-I RNA editing analysis using haplotype-tagged long-read RNA-seq BAM files
produced by longcallR. Tests whether the editing frequency at known REDIportal sites differs
significantly between haplotype 1 (H1) and haplotype 2 (H2).

For each site in the REDIportal database:
  - On the '+' strand: count A (unedited) and G (edited) bases from each haplotype
  - On the '-' strand: count T (unedited) and C (edited) bases from each haplotype
  - Apply a beta-binomial test to detect haplotype-specific editing frequency differences
  - Apply Benjamini-Hochberg multiple testing correction

REDIportal indexing
-------------------
For fast region-based queries the REDIportal file must be bgzipped and tabix-indexed.
If you pass a plain BED or TSV the script will do this for you automatically (requires
bgzip and tabix on PATH), producing <file>.gz and <file>.gz.tbi alongside the original.
Subsequent runs reuse the index — the indexing is a one-time cost.

Usage:
    python longcallR-asediting.py \
        -b phased.bam \
        -r REDIportal.bed \  # or .bed.gz / TABLE1_hg38.txt / TABLE1_hg38.txt.gz
        -f reference.fa \
        -o output_prefix \
        -t 8 \
        [--annotation annotation.gtf] \
        [--min_coverage 10] \
        [--min_editing_rate 0.0] \
        [--overdispersion 0.001] \
        [--gene_types protein_coding lncRNA] \
        [--dna_vcf dna.vcf] \
        [--rna_vcf longcallR_phased.vcf]
"""

import argparse
import concurrent.futures
import gzip
import math
import sys
from collections import defaultdict

import os
import pysam
from scipy.stats import betabinom
from statsmodels.stats.multitest import multipletests


# ---------------------------------------------------------------------------
# Beta-binomial test
# ---------------------------------------------------------------------------

def convert_mu_rho_to_alpha_beta(mu, rho):
    phi = (1 - rho) / rho - 1
    alpha = mu * phi
    beta = (1 - mu) * phi
    return alpha, beta


def beta_binomial_p_value(k_obs, n, mu, rho, alternative="two-sided"):
    """
    Beta-binomial test for allele-specific signal.

    k_obs: observed count for one allele (e.g. H1 edited reads)
    n:     total reads for that haplotype
    mu:    expected proportion under null
    rho:   overdispersion parameter

    mu is clamped to (1e-6, 1-1e-6) to prevent degenerate distributions.
    When mu->0 or mu->1 (fully unedited or fully edited on both haplotypes),
    the beta distribution becomes undefined. Those sites are invariant and the
    caller should return p=1.0 before calling this function.
    """
    if n == 0:
        return 1.0

    # Clamp mu strictly away from 0 and 1
    mu = max(1e-6, min(1.0 - 1e-6, mu))

    alpha, beta_param = convert_mu_rho_to_alpha_beta(mu, rho)

    if alpha <= 0 or beta_param <= 0:
        return 1.0

    bb = betabinom(n, alpha, beta_param)
    p_obs = bb.pmf(k_obs)

    if alternative == "two-sided":
        pmf_values = [bb.pmf(k) for k in range(n + 1)]
        p_value = sum(p for p in pmf_values if p <= p_obs + 1e-15)
    elif alternative == "greater":
        p_value = float(bb.sf(k_obs - 1))
    elif alternative == "less":
        p_value = float(bb.cdf(k_obs))
    else:
        raise ValueError("Invalid alternative. Choose 'two-sided', 'greater', or 'less'.")

    return float(p_value)


# ---------------------------------------------------------------------------
# BAM region helpers
# ---------------------------------------------------------------------------

def get_bam_chroms(bam_file):
    """
    Return the set of chromosome names with at least one mapped read,
    using the BAM index (near-instantaneous).
    """
    chroms = set()
    with pysam.AlignmentFile(bam_file, "rb") as bam:
        for stat in bam.get_index_statistics():
            if stat.mapped > 0:
                chroms.add(stat.contig)
    print(f"  BAM covers {len(chroms)} chromosome(s): {', '.join(sorted(chroms))}", flush=True)
    return chroms


def get_covered_intervals(bam_file):
    """
    Single lightweight pass through the BAM (coordinates only, no sequence or
    tags loaded) to build merged covered intervals per chromosome.

    Returns dict: chrom -> list of (start, end) tuples, 0-based half-open,
    sorted and non-overlapping.
    """
    print("Scanning BAM to find read-covered intervals...", flush=True)
    raw = defaultdict(list)

    with pysam.AlignmentFile(bam_file, "rb") as bam:
        for read in bam.fetch():
            if read.is_unmapped or read.reference_end is None:
                continue
            raw[read.reference_name].append(
                (read.reference_start, read.reference_end)
            )

    merged = {}
    total = 0
    for chrom, intervals in raw.items():
        intervals.sort()
        m = []
        for start, end in intervals:
            if m and start <= m[-1][1]:
                m[-1] = (m[-1][0], max(m[-1][1], end))
            else:
                m.append((start, end))
        merged[chrom] = m
        total += len(m)

    print(f"  {total} merged interval(s) across {len(merged)} chromosome(s)", flush=True)
    return merged


# ---------------------------------------------------------------------------
# REDIportal tabix indexing and loading
# ---------------------------------------------------------------------------

def load_rediportal_tabix(rediportal_file, covered_intervals):
    """
    Use tabix to fetch only REDIportal sites overlapping BAM-covered intervals.
    The file must already be bgzipped and tabix-indexed (.gz + .gz.tbi).

    Expects REDIportal TABLE1 TSV column layout:
        col0=Accession  col1=Region   col2=Position  col3=Ref  col4=Ed
        col5=Strand     col6=db       col7=type       col8=dbsnp col9=repeat
    Position is 1-based; converted to 0-based internally.

    Returns:
        sites: dict (chrom, pos_0based) -> strand
        site_meta: dict (chrom, pos_0based) -> {"repeat_type": str, "repeat_name": str}
            repeat_type: col7 e.g. "ALU", "REP", "NONREP"
            repeat_name: col9 e.g. "SINE/AluSx1"
    """
    print("Fetching REDIportal sites overlapping read coverage via tabix...", flush=True)
    sites = {}
    site_meta = {}
    tabix = pysam.TabixFile(rediportal_file)

    for chrom, intervals in covered_intervals.items():
        for start, end in intervals:
            try:
                for row in tabix.fetch(chrom, start, end):
                    parts = row.split("\t")
                    if len(parts) < 10:
                        continue
                    strand = parts[5]
                    if strand not in ("+", "-"):
                        continue
                    pos_0based = int(parts[2]) - 1
                    sites[(chrom, pos_0based)] = strand
                    site_meta[(chrom, pos_0based)] = {
                        "repeat_type": parts[7],
                        "repeat_name": parts[9],
                    }
            except ValueError:
                continue

    tabix.close()
    n_rep = sum(1 for m in site_meta.values() if m["repeat_type"] != "NONREP")
    print(f"  {len(sites)} REDIportal site(s) overlap read coverage "
          f"({n_rep} in repeat elements)", flush=True)
    return sites, site_meta


# ---------------------------------------------------------------------------
# Optional: gene annotation for site labelling
# ---------------------------------------------------------------------------

def load_gene_annotation(annotation_file, gene_types, bam_chroms):
    """
    Build chrom -> IntervalTree of gene intervals, restricted to BAM chromosomes.
    """
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
            attributes = parts[8]
            attr_dict = {}
            if file_type == "gtf":
                for attr in attributes.strip().split(";"):
                    attr = attr.strip()
                    if not attr:
                        continue
                    kv = attr.split(" ", 1)
                    if len(kv) == 2:
                        attr_dict[kv[0]] = kv[1].replace('"', '')
            else:
                for attr in attributes.strip().split(";"):
                    attr = attr.strip()
                    if "=" not in attr:
                        continue
                    k, v = attr.split("=", 1)
                    attr_dict[k] = v.replace('"', '')

            gene_type = attr_dict.get("gene_type", attr_dict.get("gene_biotype", ""))
            if gene_type not in gene_types:
                continue
            gene_name = attr_dict.get("gene_name", attr_dict.get("gene_id", "."))
            start = int(parts[3])   # 1-based inclusive
            end   = int(parts[4])   # 1-based inclusive
            trees[chrom].add(Interval(start, end + 1, (gene_name,)))
            n_genes += 1

    print(f"  Loaded {n_genes} gene record(s) on BAM chromosome(s)", flush=True)
    return trees


def get_gene_name_for_site(trees, chrom, pos_0based):
    if trees is None:
        return "."
    hits = trees[chrom].overlap(pos_0based + 1, pos_0based + 2)
    if not hits:
        return "."
    return min(hits, key=lambda iv: iv.end - iv.begin).data[0]


# ---------------------------------------------------------------------------
# Optional: DNA VCF filter
# ---------------------------------------------------------------------------

def load_dna_vcf(vcf_file):
    dna_vcfs = {}
    with pysam.VariantFile(vcf_file) as vcf:
        for record in vcf.fetch():
            gt = record.samples[0]["GT"]
            if len(record.ref) != 1 or any(len(alt) != 1 for alt in record.alts):
                continue
            if gt in [(0, 1), (1, 0)]:
                dna_vcfs[f"{record.contig}:{record.pos}"] = {
                    "ref": record.ref, "alt": record.alts[0],
                }
    return dna_vcfs


def load_longcallR_phased_vcf(vcf_file):
    rna_vcfs = defaultdict(list)
    with pysam.VariantFile(vcf_file) as vcf:
        for record in vcf.fetch():
            if "PASS" not in record.filter.keys():
                continue
            gt = record.samples[0]["GT"]
            if len(record.ref) != 1 or any(len(alt) != 1 for alt in record.alts):
                continue
            if gt in [(0, 1), (1, 0)] and record.samples[0].phased:
                ps = record.samples[0].get("PS", None)
                if ps and ps != ".":
                    rna_vcfs[ps].append(f"{record.contig}:{record.pos}")
    return rna_vcfs


# ---------------------------------------------------------------------------
# Core: pileup at a single editing site
# ---------------------------------------------------------------------------

def pileup_site(bam_file, chrom, pos_0based, strand,
                min_baseq=20, min_mapq=20, max_depth=100000):
    """
    Pileup one A-to-I site and return per-haplotype base counts.

    '+' strand: A = unedited, G = edited
    '-' strand: T = unedited, C = edited  (A->I on minus strand reads as T->C)

    Returns dict: {phase_set -> {haplotype -> {"ref": int, "edit": int}}}
    """
    ref_base  = "A" if strand == "+" else "T"
    edit_base = "G" if strand == "+" else "C"

    counts = defaultdict(lambda: defaultdict(lambda: {"ref": 0, "edit": 0}))

    with pysam.AlignmentFile(bam_file, "rb") as bam:
        for col in bam.pileup(
            chrom, pos_0based, pos_0based + 1,
            min_base_quality=min_baseq,
            min_mapping_quality=min_mapq,
            max_depth=max_depth,
            truncate=True,
            ignore_overlaps=True,
        ):
            if col.reference_pos != pos_0based:
                continue
            for pr in col.pileups:
                if pr.is_del or pr.is_refskip:
                    continue
                aln = pr.alignment
                if not aln.has_tag("HP") or not aln.has_tag("PS"):
                    continue
                hp = aln.get_tag("HP")
                ps = aln.get_tag("PS")
                if hp not in (1, 2):
                    continue
                base = aln.query_sequence[pr.query_position].upper()
                if base == ref_base:
                    counts[ps][hp]["ref"] += 1
                elif base == edit_base:
                    counts[ps][hp]["edit"] += 1

    return counts


# ---------------------------------------------------------------------------
# Per-site analysis
# ---------------------------------------------------------------------------

def analyze_site(bam_file, chrom, pos_0based, strand, gene_name,
                 min_coverage, min_editing_rate, overdispersion,
                 dna_vcfs, rna_vcfs):
    """
    Analyse one REDIportal site for allele-specific editing.
    Returns a result dict or None if the site fails filters.
    """
    counts = pileup_site(bam_file, chrom, pos_0based, strand)
    if not counts:
        return None

    # Pick the phase set with the most phased reads at this site
    best_ps = max(
        counts,
        key=lambda ps: sum(hc["ref"] + hc["edit"] for hc in counts[ps].values())
    )
    hap = counts[best_ps]

    h1 = hap.get(1, {"ref": 0, "edit": 0})
    h2 = hap.get(2, {"ref": 0, "edit": 0})
    h1_ref, h1_edit = h1["ref"], h1["edit"]
    h2_ref, h2_edit = h2["ref"], h2["edit"]
    h1_total = h1_ref + h1_edit
    h2_total = h2_ref + h2_edit
    total    = h1_total + h2_total

    if total < min_coverage:
        return None

    overall_edit_rate = (h1_edit + h2_edit) / total
    if overall_edit_rate < min_editing_rate:
        return None

    # Optional: require at least one DNA-anchored SNP in the same phase block
    if dna_vcfs is not None and rna_vcfs is not None:
        ps_variants = rna_vcfs.get(best_ps, [])
        if not any(
            f"{s.split(':')[0]}:{s.split(':')[1]}" in dna_vcfs
            for s in ps_variants
        ):
            return None

    # Early-exit for invariant sites to avoid degenerate beta distribution:
    #   - Fully unedited on both haplotypes -> p = 1.0
    #   - Fully edited on both haplotypes   -> p = 1.0
    #   - One haplotype has no reads         -> p = 1.0
    if h1_total == 0 or h2_total == 0:
        p_val = 1.0
    elif h1_edit == 0 and h2_edit == 0:
        p_val = 1.0
    elif h1_ref == 0 and h2_ref == 0:
        p_val = 1.0
    else:
        # Pooled editing rate as expected proportion under H0.
        # beta_binomial_p_value clamps this away from 0/1 internally.
        expected_rate = (h1_edit + h2_edit) / total
        p_val = beta_binomial_p_value(
            h1_edit, h1_total, expected_rate, overdispersion, alternative="two-sided"
        )

    h1_rate = h1_edit / h1_total if h1_total > 0 else float("nan")
    h2_rate = h2_edit / h2_total if h2_total > 0 else float("nan")
    logfc   = (math.log2((h1_edit + 1) / (h1_total + 1))
               - math.log2((h2_edit + 1) / (h2_total + 1)))

    return {
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
        "overall_rate": overall_edit_rate,
        "p_value":      p_val,
        "logfc":        logfc,
    }


# ---------------------------------------------------------------------------
# Worker for parallel execution
# ---------------------------------------------------------------------------

def _worker(args):
    (bam_file, chrom, pos_0based, strand, gene_name,
     min_coverage, min_editing_rate, overdispersion,
     dna_vcfs, rna_vcfs) = args
    try:
        return analyze_site(
            bam_file, chrom, pos_0based, strand, gene_name,
            min_coverage, min_editing_rate, overdispersion,
            dna_vcfs, rna_vcfs,
        )
    except Exception as e:
        print(f"  Warning: error at {chrom}:{pos_0based+1} - {e}",
              file=sys.stderr, flush=True)
        return None



# ---------------------------------------------------------------------------
# Repeat-element clustering
# ---------------------------------------------------------------------------

def cluster_repeat_elements(site_results, site_meta, overdispersion, max_gap=300, min_sites=2):
    """
    Group per-site results into clusters where:
      - Sites are on the same chromosome, same strand, same repeat_name
      - Consecutive sites are within max_gap bp of each other
      - repeat_type != "NONREP"

    For each cluster, pool H1 and H2 counts across all member sites,
    using the dominant phase set (most reads) per site.

    Returns list of cluster dicts, one per cluster passing min_sites.
    """
    # Index results by (chrom, pos_1based) for fast lookup
    result_by_pos = {}
    for r in site_results:
        result_by_pos[(r["chrom"], r["pos_1based"])] = r

    # Collect repeat sites sorted for clustering
    repeat_sites = []
    for (chrom, pos_0based), meta in site_meta.items():
        if meta["repeat_type"] == "NONREP":
            continue
        pos_1based = pos_0based + 1
        # Only include sites that passed per-site filters (have a result)
        if (chrom, pos_1based) not in result_by_pos:
            continue
        repeat_sites.append((
            chrom,
            pos_1based,
            result_by_pos[(chrom, pos_1based)]["strand"],
            meta["repeat_type"],
            meta["repeat_name"],
        ))

    # Sort by chrom, strand, repeat_name, position for grouping
    repeat_sites.sort(key=lambda x: (x[0], x[2], x[4], x[1]))

    clusters = []
    i = 0
    while i < len(repeat_sites):
        chrom, pos, strand, rep_type, rep_name = repeat_sites[i]
        # Start a new cluster
        cluster_sites = [pos]
        j = i + 1
        while j < len(repeat_sites):
            c2, p2, s2, rt2, rn2 = repeat_sites[j]
            # Must match chrom, strand, repeat_name and be within max_gap
            if c2 == chrom and s2 == strand and rn2 == rep_name and (p2 - cluster_sites[-1]) <= max_gap:
                cluster_sites.append(p2)
                j += 1
            else:
                break

        if len(cluster_sites) >= min_sites:
            # Pool counts across sites in this cluster
            h1_ref_total = h1_edit_total = 0
            h2_ref_total = h2_edit_total = 0
            phase_sets = []
            for site_pos in cluster_sites:
                r = result_by_pos[(chrom, site_pos)]
                h1_ref_total  += r["h1_ref"]
                h1_edit_total += r["h1_edit"]
                h2_ref_total  += r["h2_ref"]
                h2_edit_total += r["h2_edit"]
                phase_sets.append(str(r["phase_set"]))

            h1_total = h1_ref_total + h1_edit_total
            h2_total = h2_ref_total + h2_edit_total
            total    = h1_total + h2_total

            h1_rate = h1_edit_total / h1_total if h1_total > 0 else float("nan")
            h2_rate = h2_edit_total / h2_total if h2_total > 0 else float("nan")
            overall = (h1_edit_total + h2_edit_total) / total if total > 0 else 0.0
            logfc   = (math.log2((h1_edit_total + 1) / (h1_total + 1))
                       - math.log2((h2_edit_total + 1) / (h2_total + 1)))

            # Dominant phase set (most common across sites in cluster)
            dominant_ps = max(set(phase_sets), key=phase_sets.count)

            # p-value on pooled counts
            if h1_total == 0 or h2_total == 0:
                p_val = 1.0
            elif h1_edit_total == 0 and h2_edit_total == 0:
                p_val = 1.0
            elif h1_ref_total == 0 and h2_ref_total == 0:
                p_val = 1.0
            else:
                expected_rate = (h1_edit_total + h2_edit_total) / total
                p_val = beta_binomial_p_value(
                    h1_edit_total, h1_total, expected_rate,
                    overdispersion, alternative="two-sided"
                )

            clusters.append({
                "chrom":       chrom,
                "start":       min(cluster_sites),
                "end":         max(cluster_sites),
                "strand":      strand,
                "repeat_type": rep_type,
                "repeat_name": rep_name,
                "n_sites":     len(cluster_sites),
                "sites":       ",".join(str(p) for p in cluster_sites),
                "phase_set":   dominant_ps,
                "h1_ref":      h1_ref_total,
                "h1_edit":     h1_edit_total,
                "h2_ref":      h2_ref_total,
                "h2_edit":     h2_edit_total,
                "h1_rate":     h1_rate,
                "h2_rate":     h2_rate,
                "overall_rate": overall,
                "logfc":       logfc,
                "p_value":     p_val,
            })

        i = j if j > i else i + 1

    return clusters


def write_cluster_results(clusters, output_prefix, overdispersion):
    """
    Apply BH correction across clusters and write the cluster TSV.
    """
    if not clusters:
        print("No repeat-element clusters found.", flush=True)
        return

    p_values = [c["p_value"] for c in clusters]
    reject, adjusted_p, _, _ = multipletests(p_values, alpha=0.05, method="fdr_bh")
    for i, c in enumerate(clusters):
        c["p_value_adj"] = adjusted_p[i]
        c["significant"] = reject[i]

    clusters.sort(key=lambda c: c["p_value_adj"])

    out_path = output_prefix + ".asediting_clusters.tsv"
    print(f"Writing cluster results to {out_path}...", flush=True)
    with open(out_path, "w") as f:
        f.write(
            "#Chrom	Cluster_start	Cluster_end	Strand	Repeat_type	Repeat_name	"
            "N_sites	Sites	Phase_set	"
            "H1_unedited	H1_edited	H2_unedited	H2_edited	"
            "H1_editing_rate	H2_editing_rate	Overall_editing_rate	"
            "LogFC_H1_vs_H2	P_value	P_value_adj	Significant\n"
        )
        for c in clusters:
            h1r = f"{c['h1_rate']:.4f}" if not math.isnan(c["h1_rate"]) else "NA"
            h2r = f"{c['h2_rate']:.4f}" if not math.isnan(c["h2_rate"]) else "NA"
            f.write(
                f"{c['chrom']}\t{c['start']}\t{c['end']}\t{c['strand']}\t"
                f"{c['repeat_type']}\t{c['repeat_name']}\t{c['n_sites']}\t"
                f"{c['sites']}\t{c['phase_set']}\t"
                f"{c['h1_ref']}\t{c['h1_edit']}\t{c['h2_ref']}\t{c['h2_edit']}\t"
                f"{h1r}\t{h2r}\t{c['overall_rate']:.4f}\t{c['logfc']:.4f}\t"
                f"{c['p_value']:.6e}\t{c['p_value_adj']:.6e}\t"
                f"{'TRUE' if c['significant'] else 'FALSE'}\n"
            )

    n_sig = sum(1 for c in clusters if c["significant"])
    print(f"  Significant clusters (FDR<0.05): {n_sig}/{len(clusters)}", flush=True)

# ---------------------------------------------------------------------------
# Main analysis driver
# ---------------------------------------------------------------------------

def analyze(bam_file, rediportal_sites, site_meta, annotation_trees, output_prefix,
            min_coverage, min_editing_rate, overdispersion, threads,
            dna_vcfs, rna_vcfs, max_cluster_gap, min_cluster_sites):

    print(f"Building task list for {len(rediportal_sites)} site(s)...", flush=True)
    tasks = []
    for (chrom, pos_0based), strand in sorted(rediportal_sites.items()):
        gene_name = get_gene_name_for_site(annotation_trees, chrom, pos_0based)
        tasks.append((
            bam_file, chrom, pos_0based, strand, gene_name,
            min_coverage, min_editing_rate, overdispersion,
            dna_vcfs, rna_vcfs,
        ))

    print(f"Running pileup on {len(tasks)} site(s) using {threads} thread(s)...", flush=True)
    results = []
    n_done = 0
    report_every = max(1, len(tasks) // 10)
    with concurrent.futures.ProcessPoolExecutor(max_workers=threads) as executor:
        for result in executor.map(_worker, tasks):
            n_done += 1
            if n_done % report_every == 0 or n_done == len(tasks):
                print(f"  {n_done}/{len(tasks)} sites processed...", flush=True)
            if result is not None:
                results.append(result)

    print(f"Sites passing filters: {len(results)}", flush=True)

    if not results:
        print("No sites passed filters. Exiting.", flush=True)
        return

    print("Applying Benjamini-Hochberg correction...", flush=True)
    p_values = [r["p_value"] for r in results]
    reject, adjusted_p, _, _ = multipletests(p_values, alpha=0.05, method="fdr_bh")
    for i, r in enumerate(results):
        r["p_value_adj"] = adjusted_p[i]
        r["significant"] = reject[i]

    results.sort(key=lambda r: r["p_value_adj"])

    out_path = output_prefix + ".asediting.tsv"
    print(f"Writing results to {out_path}...", flush=True)
    with open(out_path, "w") as f:
        f.write(
            "#Chrom\tPos\tStrand\tGene_name\tPhase_set\t"
            "H1_unedited\tH1_edited\tH2_unedited\tH2_edited\t"
            "H1_editing_rate\tH2_editing_rate\tOverall_editing_rate\t"
            "LogFC_H1_vs_H2\tP_value\tP_value_adj\tSignificant\n"
        )
        for r in results:
            h1r = f"{r['h1_rate']:.4f}" if not math.isnan(r["h1_rate"]) else "NA"
            h2r = f"{r['h2_rate']:.4f}" if not math.isnan(r["h2_rate"]) else "NA"
            f.write(
                f"{r['chrom']}\t{r['pos_1based']}\t{r['strand']}\t{r['gene_name']}\t"
                f"{r['phase_set']}\t"
                f"{r['h1_ref']}\t{r['h1_edit']}\t{r['h2_ref']}\t{r['h2_edit']}\t"
                f"{h1r}\t{h2r}\t{r['overall_rate']:.4f}\t{r['logfc']:.4f}\t"
                f"{r['p_value']:.6e}\t{r['p_value_adj']:.6e}\t"
                f"{'TRUE' if r['significant'] else 'FALSE'}\n"
            )

    n_sig = sum(1 for r in results if r["significant"])
    print(f"Done. Significant AS editing sites (FDR<0.05): {n_sig}/{len(results)}", flush=True)

    # Repeat-element cluster analysis
    print("Clustering sites by repeat element...", flush=True)
    clusters = cluster_repeat_elements(
        results, site_meta,
        overdispersion=overdispersion,
        max_gap=max_cluster_gap,
        min_sites=min_cluster_sites,
    )
    print(f"  {len(clusters)} cluster(s) with >= {min_cluster_sites} site(s)", flush=True)
    write_cluster_results(clusters, output_prefix, overdispersion)


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Allele-specific A-to-I RNA editing analysis using longcallR haplotype-tagged BAM files."
    )
    parser.add_argument("-b", "--bam", required=True,
                        help="Haplotype-tagged BAM file from longcallR (must be indexed)")
    parser.add_argument("-r", "--rediportal", required=True,
                        help="REDIportal bgzipped BED file with tabix index "
                             "(.gz and .gz.tbi must both exist). "
                             "BED format: chrom start(0-based) end . . strand")
    parser.add_argument("-f", "--reference", required=True,
                        help="Reference genome FASTA (indexed with samtools faidx)")
    parser.add_argument("-o", "--output_prefix", required=True,
                        help="Prefix for output files")
    parser.add_argument("-a", "--annotation", default=None,
                        help="Annotation GTF/GFF3 for gene name lookup (optional)")
    parser.add_argument("--gene_types", nargs="+", default=["protein_coding", "lncRNA"],
                        help="Gene biotypes to include. Default: protein_coding lncRNA")
    parser.add_argument("-t", "--threads", type=int, default=1,
                        help="Parallel worker processes. Default: 1")
    parser.add_argument("--min_coverage", type=int, default=10,
                        help="Minimum total phased reads (H1+H2) at a site. Default: 10")
    parser.add_argument("--min_editing_rate", type=float, default=0.0,
                        help="Minimum overall editing rate to retain a site. Default: 0.0")
    parser.add_argument("--min_baseq", type=int, default=20,
                        help="Minimum base quality for pileup. Default: 20")
    parser.add_argument("--min_mapq", type=int, default=20,
                        help="Minimum mapping quality for pileup. Default: 20")
    parser.add_argument("--overdispersion", type=float, default=0.001,
                        help="Beta-binomial overdispersion parameter rho. Default: 0.001")
    parser.add_argument("--dna_vcf", default=None,
                        help="DNA VCF. With --rna_vcf, restricts to phase blocks "
                             "anchored by at least one genomic SNP.")
    parser.add_argument("--rna_vcf", default=None,
                        help="LongcallR phased RNA VCF. Must be paired with --dna_vcf.")

    parser.add_argument("--max_cluster_gap", type=int, default=300,
                        help="Maximum distance (bp) between consecutive sites "
                             "to be merged into the same repeat-element cluster. "
                             "Default: 300")
    parser.add_argument("--min_cluster_sites", type=int, default=2,
                        help="Minimum number of sites required to report a cluster. "
                             "Default: 2")

    args = parser.parse_args()

    if (args.dna_vcf is None) ^ (args.rna_vcf is None):
        parser.error("--dna_vcf and --rna_vcf must be provided together or not at all.")

    # 1. Get BAM chromosomes from index (instant)
    print("Reading BAM index...", flush=True)
    bam_chroms = get_bam_chroms(args.bam)

    # 2. Check the REDIportal file is bgzipped and tabix-indexed
    tbi_path = args.rediportal + ".tbi"
    if not args.rediportal.endswith(".gz"):
        print("Error: REDIportal file must be bgzipped (.gz).", file=sys.stderr)
        print("  Run: bgzip REDIportal.bed && tabix -p bed REDIportal.bed.gz", file=sys.stderr)
        sys.exit(1)
    if not os.path.exists(tbi_path):
        print(f"Error: tabix index not found: {tbi_path}", file=sys.stderr)
        print(f"  Run: tabix -p bed {args.rediportal}", file=sys.stderr)
        sys.exit(1)

    # 3. Get covered intervals (one fast coordinate-only read pass)
    covered = get_covered_intervals(args.bam)

    # 4. Fetch REDIportal sites via tabix restricted to covered intervals only
    rediportal_sites, site_meta = load_rediportal_tabix(args.rediportal, covered)

    if not rediportal_sites:
        print(
            "Error: no REDIportal sites overlap read coverage.\n"
            "Check chromosome naming (chr1 vs 1) and that the BAM has mapped reads.",
            file=sys.stderr,
        )
        sys.exit(1)

    # 5. Load gene annotation restricted to BAM chromosomes (optional)
    annotation_trees = load_gene_annotation(args.annotation, set(args.gene_types), bam_chroms)

    # 6. Load optional VCF filters
    dna_vcfs = rna_vcfs = None
    if args.dna_vcf and args.rna_vcf:
        print("Loading VCF files for phase-block filtering...", flush=True)
        dna_vcfs = load_dna_vcf(args.dna_vcf)
        rna_vcfs = load_longcallR_phased_vcf(args.rna_vcf)
        print(f"  DNA VCF: {len(dna_vcfs)} SNPs | RNA VCF: {len(rna_vcfs)} phase sets", flush=True)

    analyze(
        bam_file=args.bam,
        rediportal_sites=rediportal_sites,
        site_meta=site_meta,
        annotation_trees=annotation_trees,
        output_prefix=args.output_prefix,
        min_coverage=args.min_coverage,
        min_editing_rate=args.min_editing_rate,
        overdispersion=args.overdispersion,
        threads=args.threads,
        dna_vcfs=dna_vcfs,
        rna_vcfs=rna_vcfs,
        max_cluster_gap=args.max_cluster_gap,
        min_cluster_sites=args.min_cluster_sites,
    )
