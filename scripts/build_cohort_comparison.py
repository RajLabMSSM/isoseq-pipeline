#!/usr/bin/env python3
"""Comprehensive cross-cohort comparison of the three TDP-43 isoseq cohorts
(tdpkd / sun / tanaka): the DATA (design, chemistry, read-level QC, alignment) AND
the ISOFORMS they produce (per-tool novel candidates, long-read-specificity, consensus,
SQANTI, final reference, cryptic-exon sanity, leafcutter recovery).

Reads the canonical per-cohort files (tracking_numbers, longread_specific,
nanostat_collated, sanity events, leafcutter recovery summary) + a small verified
metadata table (design/chemistry/kit/knockdown, checked against the configs + HANDOFF).

Emits two deliverables next to the cross-cohort reference:
  cross_cohort_comparison.tsv       WIDE  : category, metric, description, tdpkd, sun, tanaka
  cross_cohort_comparison_long.tsv  TIDY  : cohort, category, metric, value  (for plotting)

Numbers are pulled live; only design facts not stored in any file are hard-coded (and
each is annotated with its source). Stdlib only.
"""
import os
import statistics as st
from collections import defaultdict

BASE = "/sc/arion/projects/als-omics/brooke-phd-thesis"
OUT_WIDE = f"{BASE}/plots/cross_cohort_comparison.tsv"
OUT_LONG = f"{BASE}/plots/cross_cohort_comparison_long.tsv"

COHORTS = ["tdpkd", "sun", "tanaka"]
DC = {"tdpkd": "tdpkd_nanopore",
      "sun": "sun_tdp_overexpression_nanopore",
      "tanaka": "tanaka_nanopore"}
V5 = {"tdpkd": f"{BASE}/cohorts/tdpkd/isoseq-pipeline/results/v5",
      "sun": f"{BASE}/cohorts/sun_tdp_overexpression/isoseq-pipeline/results/v5",
      "tanaka": f"{BASE}/cohorts/tanaka/isoseq-pipeline/results/v5"}

# --- verified design/metadata (source noted per row; not stored in any pipeline file) ---
# platform/chemistry/kit/preset from each cohort's *_config.yaml; design + knockdown from
# HANDOFF.md / MEMORY (short-read TARDBP levels, genotyped sample exclusions).
META = {
    "study":            {"tdpkd": "MSSM TDP-43 KD (Humphrey)", "sun": "Sun TDP-43 OE", "tanaka": "Tanaka 2025 (PRJDB19918)"},
    "perturbation":     {"tdpkd": "TARDBP knockdown", "sun": "TDP-43 overexpression", "tanaka": "TARDBP knockdown"},
    "groups":           {"tdpkd": "TDP vs CTRL", "sun": "TDP vs GFP", "tanaka": "KD vs SCR"},
    "ont_platform":     {"tdpkd": "PromethION", "sun": "PromethION", "tanaka": "GridION"},
    "ont_chemistry":    {"tdpkd": "R9.4.1", "sun": "R10.4.1", "tanaka": "R9.4.1"},
    "pychopper_kit":    {"tdpkd": "PCS109", "sun": "PCS111", "tanaka": "PCS111"},
    "minimap_preset":   {"tdpkd": "splice k14", "sun": "splice:hq k15", "tanaka": "splice k14"},
    "n_longread_samples": {"tdpkd": "11", "sun": "8", "tanaka": "9"},
    "n_shortread_samples": {"tdpkd": "12", "sun": "8", "tanaka": "9"},
    "tardbp_knockdown_pct": {"tdpkd": "~86", "sun": "n/a (OE)", "tanaka": "~59"},
    "has_unc13a_cryptic_exon": {"tdpkd": "yes", "sun": "no (OE)", "tanaka": "no"},
    "has_stmn2_cryptic_exon": {"tdpkd": "yes", "sun": "no (OE)", "tanaka": "yes"},
    "excluded_samples": {"tdpkd": "TDP-5 (aberrant diff.)", "sun": "none", "tanaka": "none"},
}
META_DESC = {
    "study": "cohort / dataset",
    "perturbation": "TDP-43 perturbation direction",
    "groups": "contrast (case vs control)",
    "ont_platform": "ONT sequencer",
    "ont_chemistry": "ONT pore chemistry (config mmi/preset)",
    "pychopper_kit": "cDNA kit primers used by pychopper",
    "minimap_preset": "minimap2 splice preset / k-mer",
    "n_longread_samples": "long-read (ONT) samples in v5",
    "n_shortread_samples": "short-read (Illumina) samples",
    "tardbp_knockdown_pct": "TARDBP depletion (short-read TPM)",
    "has_unc13a_cryptic_exon": "UNC13A cryptic exon present in data",
    "has_stmn2_cryptic_exon": "STMN2 cryptic exon present in data",
    "excluded_samples": "samples excluded from analysis",
}


def read_tracking(cohort):
    d = {}
    with open(f"{V5[cohort]}/consensus/{DC[cohort]}_tracking_numbers.tsv") as fh:
        next(fh)
        for ln in fh:
            k, v = ln.rstrip("\n").split("\t")
            d[k] = v
    return d


def read_nanostat(cohort):
    """cohort-level aggregates from the per-sample x stage nanostat_collated."""
    rows = []
    with open(f"{V5[cohort]}/{DC[cohort]}_nanostat_collated.tsv") as fh:
        hdr = fh.readline().rstrip("\n").split("\t")
        idx = {c: i for i, c in enumerate(hdr)}
        for ln in fh:
            f = ln.rstrip("\n").split("\t")
            rows.append(f)

    def col(f, name):
        v = f[idx[name]]
        try:
            return float(v)
        except (ValueError, IndexError):
            return None

    raw = {r[idx["sample"]]: r for r in rows if r[idx["stage"]] == "raw"}
    aln = {r[idx["sample"]]: r for r in rows if r[idx["stage"]] == "aligned"}

    def agg(stage_map, name, fn=st.median):
        vals = [col(f, name) for f in stage_map.values() if col(f, name) is not None]
        return fn(vals) if vals else None

    # read alignment fraction per sample = aligned reads / raw reads
    aln_frac = []
    for s in raw:
        r = col(raw[s], "number_of_reads")
        a = col(aln[s], "number_of_reads") if s in aln else None
        if r and a:
            aln_frac.append(100 * a / r)

    return {
        "n_reads_raw_total": sum(col(f, "number_of_reads") for f in raw.values()),
        "n_bases_raw_total_gb": sum(col(f, "number_of_bases") for f in raw.values()) / 1e9,
        "median_read_length_raw": agg(raw, "median_read_length"),
        "mean_read_length_raw": agg(raw, "mean_read_length", st.mean),
        "n50_raw": agg(raw, "n50"),
        "mean_qual_raw": agg(raw, "mean_qual", st.mean),
        "pct_reads_aligned": st.mean(aln_frac) if aln_frac else None,
        "median_read_length_aligned": agg(aln, "median_read_length"),
        "average_identity_aligned": agg(aln, "average_identity", st.mean),
    }


def read_events(cohort):
    d = {}
    p = f"{V5[cohort]}/consensus/{DC[cohort]}_sanity_genes.tsv.events.tsv"
    with open(p) as fh:
        hdr = fh.readline().rstrip("\n").split("\t")
        i = {c: n for n, c in enumerate(hdr)}
        for ln in fh:
            f = ln.rstrip("\n").split("\t")
            d[f[i["gene"]]] = f[i["verdict"]]
    return d


COHORT_DIR = {"tdpkd": "tdpkd", "sun": "sun_tdp_overexpression", "tanaka": "tanaka"}


def read_leafcutter(cohort):
    """DS-significant + DS-novel any-tool recovery %, if a summary exists."""
    p = (f"{BASE}/cohorts/{COHORT_DIR[cohort]}/isoseq-pipeline/analysis/"
         f"{cohort}_leafcutter_longread_recovery_v5_summary.tsv")
    if not os.path.exists(p):
        return {}
    out = {}
    with open(p) as fh:
        hdr = fh.readline().rstrip("\n").split("\t")
        i = {c: n for n, c in enumerate(hdr)}
        for ln in fh:
            f = ln.rstrip("\n").split("\t")
            js = f[i["junction_set"]]
            out[f"lc_{js}_n"] = f[i["n_junctions"]]
            out[f"lc_{js}_any_pct"] = f[i["any_pct"]]
            out[f"lc_{js}_stringtie_pct"] = f[i["stringtie_pct"]]
            out[f"lc_{js}_bambu_pct"] = f[i["bambu_pct"]]
    return out


def main():
    trk = {c: read_tracking(c) for c in COHORTS}
    nano = {c: read_nanostat(c) for c in COHORTS}
    ev = {c: read_events(c) for c in COHORTS}
    lc = {c: read_leafcutter(c) for c in COHORTS}

    def fmt(v, dp=1):
        if v is None:
            return "NA"
        if isinstance(v, float):
            return f"{v:.{dp}f}"
        return str(v)

    # rows = (category, metric, description, {cohort: value})
    rows = []

    # --- 1. DESIGN / METADATA ---
    for m, desc in META_DESC.items():
        rows.append(("1_design", m, desc, {c: META[m][c] for c in COHORTS}))

    # --- 2. READ-LEVEL DATA (NanoStat) ---
    nano_rows = [
        ("n_reads_raw_total", "total raw ONT reads", 0),
        ("n_bases_raw_total_gb", "total raw bases (Gb)", 1),
        ("median_read_length_raw", "median raw read length (bp)", 0),
        ("mean_read_length_raw", "mean raw read length (bp)", 0),
        ("n50_raw", "raw read N50 (bp)", 0),
        ("mean_qual_raw", "mean raw read quality (Q)", 1),
        ("pct_reads_aligned", "% reads aligned (minimap2)", 1),
        ("median_read_length_aligned", "median aligned read length (bp)", 0),
        ("average_identity_aligned", "mean aligned identity (%)", 1),
    ]
    for key, desc, dp in nano_rows:
        rows.append(("2_reads", key, desc, {c: fmt(nano[c][key], dp) for c in COHORTS}))

    # --- 3. PER-TOOL NOVEL CANDIDATES ---
    for t in ["bambu", "isoquant", "stringtie"]:
        rows.append(("3_pertool_novel", f"novel_candidates_{t}",
                     f"{t} novel candidate isoforms (class != =,c)",
                     {c: trk[c].get(f"novel_candidates_{t}", "NA") for c in COHORTS}))
    rows.append(("3_pertool_novel", "novel_candidates_total", "all tools, pre-consensus",
                 {c: trk[c].get("novel_candidates_total", "NA") for c in COHORTS}))

    # --- 4. LONG-READ-SPECIFIC (>=1 junction absent from short-read STAR SJ) ---
    for t in ["bambu", "isoquant", "stringtie"]:
        rows.append(("4_longread_specific", f"lr_specific_{t}",
                     f"{t} novel isoforms with >=1 junction NOT in short reads",
                     {c: trk[c].get(f"lr_specific_{t}", "NA") for c in COHORTS}))
    rows.append(("4_longread_specific", "lr_specific_total", "all tools long-read-specific",
                 {c: trk[c].get("lr_specific_total", "NA") for c in COHORTS}))
    rows.append(("4_longread_specific", "sr_supported_total",
                 "novel isoforms fully short-read-supported",
                 {c: trk[c].get("sr_supported_total", "NA") for c in COHORTS}))

    # --- 5. CONSENSUS + OVERLAP ---
    rows.append(("5_consensus", "consensus_novel_total", "after end-aware collapse (union)",
                 {c: trk[c].get("consensus_novel_total", "NA") for c in COHORTS}))
    rows.append(("5_consensus", "overlap_all_three_tools",
                 "novel isoforms supported by all 3 tools",
                 {c: trk[c].get("overlap_bambu+isoquant+stringtie", "NA") for c in COHORTS}))

    # --- 6. SQANTI + FINAL REFERENCE ---
    for key, desc in [("sqanti_pass", "novel passing SQANTI filter"),
                      ("sqanti_removed", "novel removed by SQANTI"),
                      ("final_novel", "final novel isoforms in reference"),
                      ("reference_transcripts", "GENCODE v50 transcripts (context)"),
                      ("final_total", "final reference transcripts (GENCODE + novel)")]:
        rows.append(("6_sqanti_final", key, desc,
                     {c: trk[c].get(key, "NA") for c in COHORTS}))

    # --- 7. CRYPTIC-EXON SANITY ---
    for gene in ["STMN2", "UNC13A"]:
        rows.append(("7_cryptic_sanity", f"{gene}_cryptic_exon_verdict",
                     f"{gene} cryptic exon status in final reference",
                     {c: ev[c].get(gene, "NA") for c in COHORTS}))

    # --- 8. LEAFCUTTER SHORT-READ -> LONG-READ RECOVERY (tdpkd, sun) ---
    lc_rows = [
        ("lc_DS_significant_n", "differentially-spliced junctions tested"),
        ("lc_DS_significant_any_pct", "% DS junctions recovered by any tool"),
        ("lc_DS_novel_n", "DS junctions that are novel (not GENCODE)"),
        ("lc_DS_novel_any_pct", "% DS-novel recovered by any tool"),
        ("lc_DS_novel_stringtie_pct", "% DS-novel recovered by StringTie --mix"),
        ("lc_DS_novel_bambu_pct", "% DS-novel recovered by Bambu (pure long-read)"),
    ]
    for key, desc in lc_rows:
        rows.append(("8_leafcutter_recovery", key, desc,
                     {c: lc[c].get(key, "NA (no leafcutter run)") for c in COHORTS}))

    # write WIDE
    with open(OUT_WIDE, "w") as out:
        out.write("category\tmetric\tdescription\t" + "\t".join(COHORTS) + "\n")
        for cat, metric, desc, vals in rows:
            out.write(f"{cat}\t{metric}\t{desc}\t" + "\t".join(str(vals[c]) for c in COHORTS) + "\n")

    # write LONG (tidy)
    with open(OUT_LONG, "w") as out:
        out.write("cohort\tcategory\tmetric\tdescription\tvalue\n")
        for cat, metric, desc, vals in rows:
            for c in COHORTS:
                out.write(f"{c}\t{cat}\t{metric}\t{desc}\t{vals[c]}\n")

    print(f"wrote {OUT_WIDE}")
    print(f"wrote {OUT_LONG}")
    print(f"  {len(rows)} metrics x {len(COHORTS)} cohorts")


if __name__ == "__main__":
    main()
