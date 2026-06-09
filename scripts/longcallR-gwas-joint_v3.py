#!/usr/bin/env python
r"""
longcallR-gwas-joint.py

Joint analysis of GWAS risk allele effects on allele-specific expression,
splicing, or RNA editing across multiple longcallR samples.

For each sample, this script runs the same per-sample annotation as
longcallR-gwas-annotate.py, then pools results across all samples to ask:
for each (gene x GWAS SNP) pair, does the risk allele consistently show
higher or lower expression/splicing/editing?

Statistical approach
--------------------
For each (gene x GWAS SNP) pair observed in >= min_samples samples:

  - n_samples:      number of samples where this gene was in a phase set
                    containing this GWAS SNP
  - n_concordant:   samples where risk_logfc > 0 (risk allele -> higher signal)
  - n_discordant:   samples where risk_logfc < 0 (risk allele -> lower signal)
  - mean_risk_logfc: mean of risk_logfc across samples (meta-analytic effect)
  - binomial_p:     two-sided binomial test p-value (H0: P(concordant) = 0.5)

No multiple testing correction is applied at this stage — the output is
intended for exploratory review. The GWAS loci are already pre-filtered at
genome-wide significance (P < 5e-8), and phase set sizes are typically
gene-sized so within-locus redundancy is low.

File structure expected
-----------------------
    results/<sample>/phased/<sample>.vcf
    results/<sample>/phased/<sample>.ase.tsv      (or asj.tsv / asediting.tsv)

Usage
-----
    python longcallR-gwas-joint.py \
        --results-dir results \
        --mode ase \
        --gwas gwas_sumstats.tsv.gz \
        --out output_prefix \
        [--gwas-p 5e-8] \
        [--min-samples 3] \
        [--keep-ambiguous] \
        [--samples sample1 sample2 ...]   # optional: subset of samples
"""

import argparse
import gzip
import math
import os
import sys
from collections import defaultdict
from scipy.stats import binomtest

# ---------------------------------------------------------------------------
# Import shared functions from longcallR-gwas-annotate.py
# ---------------------------------------------------------------------------
# Rather than duplicating code, we import the annotation functions directly.
# This requires longcallR-gwas-annotate.py to be in the same directory or
# on PYTHONPATH. We import it as a module by manipulating sys.path.

def _import_annotate():
    script_dir = os.path.dirname(os.path.abspath(__file__))
    if script_dir not in sys.path:
        sys.path.insert(0, script_dir)
    try:
        import importlib.util
        spec = importlib.util.spec_from_file_location(
            "longcallR_gwas_annotate",
            os.path.join(script_dir, "longcallR-gwas-annotate.py")
        )
        mod = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(mod)
        return mod
    except Exception as e:
        print(f"Error: could not import longcallR-gwas-annotate.py: {e}",
              file=sys.stderr)
        print("Make sure longcallR-gwas-annotate.py is in the same directory.",
              file=sys.stderr)
        sys.exit(1)


# ---------------------------------------------------------------------------
# Sample discovery
# ---------------------------------------------------------------------------

def discover_samples(results_dir, mode, sample_list=None):
    """
    Find all samples that have both a VCF and a longcallR result file.

    Expected file structure:
        results/<sample>/phased/<sample>.vcf
        results/<sample>/phased/<sample>.<mode>.tsv

    Returns list of (sample_name, vcf_path, result_path) tuples.
    """
    found = []
    if sample_list:
        candidates = sample_list
    else:
        try:
            candidates = sorted(os.listdir(results_dir))
        except FileNotFoundError:
            print(f"Error: results directory not found: {results_dir}",
                  file=sys.stderr)
            sys.exit(1)

    for sample in candidates:
        phased_dir = os.path.join(results_dir, sample, "phased")
        vcf_path    = os.path.join(phased_dir, f"{sample}.vcf")
        vcf_gz_path = vcf_path + ".gz"
        result_path = os.path.join(phased_dir, f"{sample}.{mode}.tsv")

        # Accept plain or bgzipped VCF
        actual_vcf = vcf_gz_path if os.path.exists(vcf_gz_path) else vcf_path

        if not os.path.exists(actual_vcf):
            print(f"  Skipping {sample}: no VCF found at {actual_vcf}",
                  flush=True)
            continue
        if not os.path.exists(result_path):
            print(f"  Skipping {sample}: no {mode} result at {result_path}",
                  flush=True)
            continue

        found.append((sample, actual_vcf, result_path))

    return found


# ---------------------------------------------------------------------------
# Per-sample annotation (wraps the single-sample logic)
# ---------------------------------------------------------------------------

def annotate_sample(sample, vcf_path, result_path, gwas_hits,
                    keep_ambiguous, ann):
    """
    Run the per-sample GWAS annotation for one sample.

    Returns list of dicts with keys:
        sample, gene_name, chrom, phase_set,
        gwas_snp, gwas_effect, gwas_pvalue,
        risk_haplotype, risk_logfc, ase_logfc, ase_pvalue,
        direction_concordant, strand_ambiguous, n_gwas_snps_in_ps
    """
    # Map GWAS hits to phase sets via VCF
    phase_assignments, phase_set_to_hits = ann.load_vcf_phase_assignments(
        vcf_path, gwas_hits
    )

    if not phase_assignments:
        return []

    # Load longcallR results
    header, rows, mode = ann.load_longcallr_results(result_path)

    # Annotate
    annotated = ann.annotate_results(
        rows, phase_set_to_hits, phase_assignments, keep_ambiguous
    )

    # Extract relevant rows (those with a GWAS annotation)
    results = []
    for row in annotated:
        if row.get("gwas_snp", "NA") == "NA":
            continue
        if row.get("risk_logfc", "NA") == "NA":
            continue

        # Extract gene/junction name robustly.
        # For ASJ, include both junction ID and gene name so results are
        # identifiable in the joint output (Gene_name alone is not unique
        # since one gene can have many junctions).
        junction = row.get("Junction", "")
        gene     = row.get("Gene_name", row.get("gene_name", ""))
        if junction and gene:
            gene_name = f"{gene}:{junction}"
        elif junction:
            gene_name = junction
        elif gene:
            gene_name = gene
        else:
            gene_name = "."

        # Extract chrom
        chrom = row.get("Chr", row.get("Chrom", row.get("#Chrom", ".")))

        # Extract site position for editing mode — used to make the grouping
        # key site-specific rather than gene-specific. Without this, all
        # editing sites in the same gene x GWAS SNP combination collapse into
        # one row with n_samples = n_sites x n_samples, not n_samples.
        # For ASE/ASJ this is "." since they are already gene/junction-level.
        site_pos = row.get("Pos", row.get("pos", "."))

        # Extract ASE p-value (use hierarchical if available, else raw)
        ase_pvalue = "NA"
        for col in ("P_element_adj", "P_value_adj", "P_value"):
            if col in row and row[col] not in ("NA", ""):
                ase_pvalue = row[col]
                break

        ase_logfc = ann.get_logfc(row)
        risk_logfc_str = row.get("risk_logfc", "NA")

        try:
            risk_logfc = float(risk_logfc_str)
        except ValueError:
            continue

        results.append({
            "sample":             sample,
            "gene_name":          gene_name,
            "chrom":              chrom,
            "site_pos":           str(site_pos),
            "phase_set":          str(ann.get_phase_set(row)),
            "gwas_snp":           row.get("gwas_snp", "NA"),
            "gwas_rsid":          row.get("gwas_rsid", "NA"),
            "gwas_allele1":       row.get("gwas_allele1", "NA"),
            "gwas_allele2":       row.get("gwas_allele2", "NA"),
            "gwas_effect":        row.get("gwas_effect", "NA"),
            "gwas_pvalue":        row.get("gwas_pvalue", "NA"),
            "risk_haplotype":     row.get("risk_haplotype", "NA"),
            "risk_logfc":         risk_logfc,
            "ase_logfc":          str(ase_logfc) if ase_logfc is not None else "NA",
            "ase_pvalue":         ase_pvalue,
            "direction_concordant": row.get("direction_concordant", "NA"),
            "strand_ambiguous":   row.get("strand_ambiguous", "NA"),
            "n_gwas_snps_in_ps":  row.get("n_gwas_snps_in_ps", "NA"),
        })

    return results


# ---------------------------------------------------------------------------
# Joint aggregation and testing
# ---------------------------------------------------------------------------

def compute_se_from_p_and_logfc(logfc, p_value):
    """
    Derive the standard error of a logFC estimate from its p-value using the
    relationship: Z = logfc / SE, where Z = sign(logfc) * |Phi^{-1}(p/2)|.

    This avoids needing to compute SE from raw read counts and works uniformly
    across ASE, ASJ, and ASediting outputs.

    Edge cases:
      - p_value == 0: clamp to a very small value to avoid infinite Z
      - p_value == 1: Z = 0, SE = inf (uninformative, weight = 0)
      - logfc == 0: SE = inf (uninformative, weight = 0)

    Returns SE (float), or inf if the estimate is uninformative.
    """
    import math
    from scipy.stats import norm

    try:
        logfc = float(logfc)
        p     = float(p_value)
    except (ValueError, TypeError):
        return float("inf")

    if p <= 0:
        p = 1e-300        # clamp to avoid log(0)
    if p >= 1:
        return float("inf")   # p=1 gives Z=0, uninformative
    if logfc == 0.0:
        return float("inf")   # zero effect, uninformative

    # Two-sided: Z = Phi^{-1}(1 - p/2)
    z = norm.ppf(1.0 - p / 2.0)
    if z <= 0:
        return float("inf")

    # SE = |logfc| / Z  (sign of logfc handled separately in weighted mean)
    return abs(logfc) / z


def inverse_variance_meta(risk_logfcs, se_values):
    """
    Inverse-variance weighted meta-analysis across samples.

    Each sample i contributes:
        weight_i  = 1 / SE_i^2
        beta_i    = risk_logfc_i

    Combined estimate:
        beta_IVW  = sum(beta_i * weight_i) / sum(weight_i)
        SE_IVW    = 1 / sqrt(sum(weight_i))
        Z_IVW     = beta_IVW / SE_IVW
        p_IVW     = two-sided p-value from N(0,1)

    Samples with SE = inf (p=1 or logfc=0) get weight 0 and are excluded.

    Returns (beta_ivw, se_ivw, z_ivw, p_ivw, n_weighted) where n_weighted is
    the number of samples with finite weight that contributed to the estimate.
    """
    from scipy.stats import norm

    weights = []
    weighted_betas = []
    n_weighted = 0

    for beta, se in zip(risk_logfcs, se_values):
        if not math.isfinite(se) or se <= 0:
            continue
        w = 1.0 / (se ** 2)
        weights.append(w)
        weighted_betas.append(beta * w)
        n_weighted += 1

    if not weights:
        return float("nan"), float("nan"), float("nan"), 1.0, 0

    sum_w  = sum(weights)
    beta_ivw = sum(weighted_betas) / sum_w
    se_ivw   = 1.0 / math.sqrt(sum_w)
    z_ivw    = beta_ivw / se_ivw
    p_ivw    = 2.0 * norm.sf(abs(z_ivw))   # two-sided

    return beta_ivw, se_ivw, z_ivw, p_ivw, n_weighted


def aggregate_and_test(all_results, min_samples):
    """
    Group results by (gene_name, site_pos, gwas_snp) across samples and
    compute two joint statistics:

    1. Inverse-variance weighted (IVW) meta-analysis of risk_logfc values.
       Each sample's risk_logfc is weighted by 1/SE^2 where SE is derived
       from the per-sample ASE p-value: SE = |logfc| / Z, Z = Phi^{-1}(1-p/2).
       This down-weights imprecise estimates (few reads, large logfc) relative
       to precise ones (many reads, small logfc), addressing the problem that
       a logFC of 0.1 with p=0.001 is more informative than logFC=2.0 with
       p=0.4. The combined Z-score and p-value (ivw_p) are the primary output.

    2. Binomial sign test on direction concordance (retained for comparison
       and for cases where p-values are unavailable).

    The grouping key includes site_pos to handle ASEditing mode (see below).
    For ASE and ASJ, site_pos is "." so the key is effectively (gene, gwas_snp).
    The grouping key uses gwas_snp (CHR:POS) rather than rsID since rsIDs
    can be NA.

    min_samples=1 is valid: for n=1 the IVW estimate is just the single sample's
    logfc and Z-score, directly interpretable. The binomial test is NA for n=1.

    Returns list of aggregated result dicts, sorted by ivw_p.
    """
    # groups[(gene, site_pos, gwas_snp)] = list of per-sample dicts
    groups = defaultdict(list)
    for r in all_results:
        key = (r["gene_name"], r.get("site_pos", "."), r["gwas_snp"])
        groups[key].append(r)

    aggregated = []
    for (gene, site_pos, gwas_snp), sample_results in groups.items():
        n = len(sample_results)
        if n < min_samples:
            continue

        risk_logfcs = [r["risk_logfc"] for r in sample_results]
        ase_pvalues = [r["ase_pvalue"] for r in sample_results]

        # Derive per-sample SE from p-value and logfc
        se_values = [
            compute_se_from_p_and_logfc(logfc, p)
            for logfc, p in zip(risk_logfcs, ase_pvalues)
        ]

        # Inverse-variance weighted meta-analysis (primary statistic)
        beta_ivw, se_ivw, z_ivw, ivw_p, n_weighted = inverse_variance_meta(
            risk_logfcs, se_values
        )

        # Binomial sign test (secondary, retained for comparison)
        n_concordant = sum(1 for v in risk_logfcs if v > 0)
        n_discordant = sum(1 for v in risk_logfcs if v < 0)
        n_zero       = sum(1 for v in risk_logfcs if v == 0)
        n_testable   = n_concordant + n_discordant

        if n_testable >= 2:
            binom_result = binomtest(n_concordant, n_testable, p=0.5,
                                     alternative="two-sided")
            binom_p = binom_result.pvalue
        else:
            binom_p = float("nan")   # not meaningful for n=1

        mean_risk_logfc = sum(risk_logfcs) / n if n > 0 else float("nan")

        # Collect metadata
        chroms     = [r["chrom"] for r in sample_results]
        chrom      = max(set(chroms), key=chroms.count)
        gwas_rsids = [r["gwas_rsid"] for r in sample_results if r["gwas_rsid"] != "NA"]
        gwas_rsid  = gwas_rsids[0] if gwas_rsids else "NA"
        gwas_a1s   = [r["gwas_allele1"] for r in sample_results]
        gwas_a1    = max(set(gwas_a1s), key=gwas_a1s.count)
        gwas_a2s   = [r["gwas_allele2"] for r in sample_results]
        gwas_a2    = max(set(gwas_a2s), key=gwas_a2s.count)
        gwas_effect = sample_results[0]["gwas_effect"]
        gwas_pvalue = sample_results[0]["gwas_pvalue"]

        # Per-sample detail columns
        sample_names  = ";".join(r["sample"] for r in sample_results)
        sample_logfcs = ";".join(f"{r['risk_logfc']:.4f}" for r in sample_results)
        sample_pvals  = ";".join(r["ase_pvalue"] for r in sample_results)
        sample_ses    = ";".join(
            f"{se:.4f}" if math.isfinite(se) else "NA"
            for se in se_values
        )

        aggregated.append({
            "gene_name":        gene,
            "site_pos":         site_pos,
            "chrom":            chrom,
            "gwas_snp":         gwas_snp,
            "gwas_rsid":        gwas_rsid,
            "gwas_allele1":     gwas_a1,
            "gwas_allele2":     gwas_a2,
            "gwas_effect":      gwas_effect,
            "gwas_pvalue":      gwas_pvalue,
            "n_samples":        n,
            "n_weighted":       n_weighted,
            "n_concordant":     n_concordant,
            "n_discordant":     n_discordant,
            "n_zero":           n_zero,
            "ivw_beta":         beta_ivw,
            "ivw_se":           se_ivw,
            "ivw_z":            z_ivw,
            "ivw_p":            ivw_p,
            "mean_risk_logfc":  mean_risk_logfc,
            "binomial_p":       binom_p,
            "samples":          sample_names,
            "per_sample_risk_logfc": sample_logfcs,
            "per_sample_se":    sample_ses,
            "per_sample_ase_pvalue": sample_pvals,
        })

    # Sort by IVW p-value (primary statistic)
    aggregated.sort(key=lambda r: r["ivw_p"] if math.isfinite(r["ivw_p"]) else 1.0)
    return aggregated


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Joint multi-sample GWAS risk allele annotation for "
                    "longcallR ASE/ASJ/ASediting results."
    )
    parser.add_argument("--results-dir", required=True,
                        help="Root results directory containing per-sample "
                             "subdirectories (results/<sample>/phased/...)")
    parser.add_argument("--mode", required=True,
                        choices=["ase", "asj", "asediting"],
                        help="longcallR output mode to analyse")
    parser.add_argument("--gwas", required=True,
                        help="GWAS summary statistics file (plain or gzipped). "
                             "Required columns: Allele1, Allele2, Effect, "
                             "PVALUE, chr, pos")
    parser.add_argument("--out", required=True,
                        help="Output file prefix")
    parser.add_argument("--gwas-p", type=float, default=5e-8,
                        help="GWAS p-value threshold. Default: 5e-8")
    parser.add_argument("--min-samples", type=int, default=1,
                        help="Minimum samples a gene x SNP pair must appear "
                             "in to be included in the joint output. Default: 1 "
                             "(all pairs reported; for n=1 binomial_p is NA)")
    parser.add_argument("--keep-ambiguous", action="store_true",
                        help="Include strand-ambiguous SNPs (A/T and C/G)")
    parser.add_argument("--samples", nargs="+", default=None,
                        help="Optional: restrict to specific sample names. "
                             "Default: all samples found in --results-dir")

    args = parser.parse_args()

    # Import shared functions from longcallR-gwas-annotate.py
    ann = _import_annotate()

    # Discover samples
    print(f"Discovering samples in {args.results_dir}...", flush=True)
    samples = discover_samples(args.results_dir, args.mode, args.samples)
    print(f"  Found {len(samples)} sample(s) with VCF + {args.mode} results",
          flush=True)
    if not samples:
        print("Error: no samples found. Check --results-dir and --mode.",
              file=sys.stderr)
        sys.exit(1)

    # Load GWAS hits once (shared across all samples)
    gwas_hits = ann.load_gwas(args.gwas, args.gwas_p)
    if not gwas_hits:
        print(f"Error: no GWAS SNPs at P < {args.gwas_p:.2e}.", file=sys.stderr)
        sys.exit(1)

    # Per-sample annotation
    all_results = []
    n_total_annotated = 0
    print(f"\nAnnotating {len(samples)} sample(s)...", flush=True)
    for i, (sample, vcf_path, result_path) in enumerate(samples):
        print(f"  [{i+1}/{len(samples)}] {sample}", flush=True)
        sample_results = annotate_sample(
            sample, vcf_path, result_path, gwas_hits,
            args.keep_ambiguous, ann
        )
        all_results.extend(sample_results)
        n_total_annotated += len(sample_results)
        print(f"    {len(sample_results)} annotated result(s)", flush=True)

    print(f"\nTotal annotated results across all samples: {n_total_annotated:,}",
          flush=True)

    # Write per-sample results (all samples combined, one row per annotation)
    persample_out = args.out + f".{args.mode}_gwas_persample.tsv"
    print(f"Writing per-sample results to {persample_out}...", flush=True)
    persample_cols = [
        "sample", "gene_name", "site_pos", "chrom", "phase_set",
        "gwas_snp", "gwas_rsid", "gwas_allele1", "gwas_allele2",
        "gwas_effect", "gwas_pvalue",
        "risk_haplotype", "risk_logfc", "ase_logfc", "ase_pvalue",
        "direction_concordant", "strand_ambiguous", "n_gwas_snps_in_ps",
    ]
    with open(persample_out, "w") as f:
        f.write("#" + "\t".join(persample_cols) + "\n")
        for r in all_results:
            f.write("\t".join(str(r.get(col, "NA")) for col in persample_cols) + "\n")

    # Joint aggregation and testing
    print(f"Aggregating across samples (min_samples={args.min_samples})...",
          flush=True)
    aggregated = aggregate_and_test(all_results, args.min_samples)
    print(f"  {len(aggregated)} gene x SNP pair(s) with >= {args.min_samples} samples",
          flush=True)

    # Summary of how many genes x SNPs were seen at each sample count
    from collections import Counter
    count_dist = Counter(r["n_samples"] for r in aggregated)
    print("  Distribution of sample counts:", flush=True)
    for n in sorted(count_dist):
        print(f"    {n} sample(s): {count_dist[n]} pair(s)", flush=True)

    # Write joint results
    joint_out = args.out + f".{args.mode}_gwas_joint.tsv"
    print(f"Writing joint results to {joint_out}...", flush=True)
    joint_cols = [
        "gene_name", "site_pos", "chrom", "gwas_snp", "gwas_rsid",
        "gwas_allele1", "gwas_allele2", "gwas_effect", "gwas_pvalue",
        "n_samples", "n_weighted",
        "n_concordant", "n_discordant", "n_zero",
        "ivw_beta", "ivw_se", "ivw_z", "ivw_p",
        "mean_risk_logfc", "binomial_p",
        "samples", "per_sample_risk_logfc", "per_sample_se",
        "per_sample_ase_pvalue",
    ]
    with open(joint_out, "w") as f:
        f.write("#" + "\t".join(joint_cols) + "\n")
        for r in aggregated:
            f.write(
                "\t".join(
                    f"{r[col]:.6e}" if isinstance(r[col], float) and col != "mean_risk_logfc"
                    else f"{r[col]:.4f}" if isinstance(r[col], float)
                    else str(r.get(col, "NA"))
                    for col in joint_cols
                ) + "\n"
            )

    # Print top hits
    if aggregated:
        print(f"\nTop 10 gene x SNP pairs by IVW p-value:", flush=True)
        print(f"  {'Gene':<20} {'GWAS SNP':<22} {'N':<4} {'Nwt':<4} "
              f"{'IVW_beta':<10} {'IVW_Z':<8} {'IVW_P':<12} {'Binom_P':<12}",
              flush=True)
        for r in aggregated[:10]:
            binom_str = f"{r['binomial_p']:.4e}" if math.isfinite(r['binomial_p']) else "NA"
            print(f"  {r['gene_name']:<20} {r['gwas_snp']:<22} "
                  f"{r['n_samples']:<4} {r['n_weighted']:<4} "
                  f"{r['ivw_beta']:<10.4f} {r['ivw_z']:<8.3f} "
                  f"{r['ivw_p']:<12.4e} {binom_str:<12}", flush=True)

    print("Done.", flush=True)
