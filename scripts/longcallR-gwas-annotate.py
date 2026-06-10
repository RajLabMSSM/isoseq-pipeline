#!/usr/bin/env python
r"""
longcallR-gwas-annotate.py

Annotates longcallR allele-specific results (ASE, ASJ, or ASediting) with
GWAS risk allele information to determine whether the haplotype showing
higher expression/splicing/editing carries the GWAS risk allele.

Background
----------
longcallR assigns reads to haplotype 1 (H1) or H2 within each phase block
(identified by the PS tag in the BAM and VCF). H1/H2 labels are arbitrary
and LOCAL — H1 in one gene has no relationship to H1 in another gene or
phase block. The phased VCF records which allele at each heterozygous SNP
belongs to H1 vs H2 within each phase block.

This script joins three pieces of information:
  1. GWAS summary statistics: which SNPs are genome-wide significant, and
     which allele increases the trait (positive Effect/beta = Allele1
     is the risk-increasing allele)
  2. Phased VCF: maps each significant GWAS SNP to a phase set and
     determines which haplotype (H1 or H2) carries the risk allele
  3. longcallR results: for each gene/junction/site, what is the H1 vs H2
     log fold change (LogFC_H1_vs_H2 > 0 means H1 has higher signal)

The key output is risk_logfc: the log fold change of the RISK haplotype
relative to the non-risk haplotype. Positive risk_logfc means the risk
allele haplotype has higher expression/splicing/editing.

GWAS file format
----------------
This script expects the following column layout (header required):
    MARKER_ID  Allele1  Allele2  Freq1  Effect  StdErr  PVALUE  ...  chr  pos  RSID

Where:
    Allele1 = effect allele (positive Effect means Allele1 increases trait)
    Allele2 = non-effect allele
    Effect  = beta (log odds ratio or linear effect)
    PVALUE  = p-value
    chr     = chromosome number (without "chr" prefix, e.g. 1 not chr1)
    pos     = 1-based position

Indels (where len(Allele1) != 1 or len(Allele2) != 1) are skipped since
haplotype phasing in longcallR is SNP-based.

VCF format
----------
Expects the longcallR phased VCF output, which uses pipe-separated GT for
phased variants (e.g. 0|1) and slash-separated for unphased (0/1).
The PS field in FORMAT gives the phase set ID, matching the Phase_set
column in longcallR output files.

Phasing convention:
    0|1 -> REF on H1, ALT on H2
    1|0 -> ALT on H1, REF on H2

Allele matching
---------------
GWAS Allele1/Allele2 are matched to VCF REF/ALT. If they don't match
directly, strand complementing is attempted (A<->T, C<->G). SNPs where
Allele1/Allele2 are a complementary pair (A/T or C/G) are flagged as
strand_ambiguous since direction cannot be confirmed without allele
frequency data.

Output columns (appended to each input longcallR row)
------------------------------------------------------
    gwas_snp        : CHR:POS of the GWAS SNP in the same phase set
    gwas_rsid       : RSID (may be NA)
    gwas_allele1    : effect allele from GWAS (risk-increasing)
    gwas_allele2    : non-effect allele from GWAS
    gwas_effect     : beta from GWAS
    gwas_pvalue     : p-value from GWAS
    risk_haplotype  : which haplotype (H1 or H2) carries the risk allele
    risk_logfc      : LogFC of risk haplotype vs non-risk haplotype
                      (positive = risk haplotype has higher signal)
    direction_concordant : TRUE if risk_logfc > 0 (risk allele associated
                      with higher expression/splicing/editing)
    strand_ambiguous: TRUE if the SNP is A/T or C/G (allele direction
                      cannot be confirmed from sequence alone)
    n_gwas_snps_in_ps : number of genome-wide significant GWAS SNPs in
                      this phase set (useful for weighting/filtering)

If multiple GWAS SNPs fall in the same phase set, the most significant
(lowest PVALUE) is used as the representative SNP. All hits are reported
in a separate summary table.

Usage
-----
    python longcallR-gwas-annotate.py \
        --input  sample.ase.tsv \
        --mode   ase \
        --vcf    longcallR_phased.vcf.gz \
        --gwas   gwas_sumstats.txt.gz \
        --out    output_prefix \
        [--gwas-p 5e-8] \
        [--keep-ambiguous]

Modes
-----
    ase       : longcallR ASE output (gene-level, Phase_set column present)
    asj       : longcallR ASJ output (junction-level)
    asediting : longcallR ASediting output (site-level)

All three modes use the same core columns: Phase_set and LogFC_H1_vs_H2.
"""

import argparse
import gzip
import math
import sys
from collections import defaultdict


# ---------------------------------------------------------------------------
# Strand complement helpers
# ---------------------------------------------------------------------------

COMPLEMENT = {"A": "T", "T": "A", "C": "G", "G": "C"}
AMBIGUOUS_PAIRS = {frozenset(["A", "T"]), frozenset(["C", "G"])}


def complement(allele):
    """Return the strand complement of a single-character allele."""
    return COMPLEMENT.get(allele.upper(), allele.upper())


def is_strand_ambiguous(a1, a2):
    """
    Return True if the allele pair is A/T or C/G.
    These pairs cannot be strand-disambiguated without allele frequency data.
    """
    return frozenset([a1.upper(), a2.upper()]) in AMBIGUOUS_PAIRS


def alleles_match(gwas_a1, gwas_a2, vcf_ref, vcf_alt):
    """
    Determine whether GWAS Allele1/Allele2 match VCF REF/ALT, with optional
    strand flip.

    Returns one of:
        "direct"   : A1=REF, A2=ALT (or A1=ALT, A2=REF) on forward strand
        "flipped"  : A1=REF, A2=ALT after strand complement
        None       : alleles don't match even after complement
    Also returns a boolean: a1_is_alt (True if A1 corresponds to VCF ALT)
    """
    a1, a2 = gwas_a1.upper(), gwas_a2.upper()
    ref, alt = vcf_ref.upper(), vcf_alt.upper()

    # Direct match
    if a1 == alt and a2 == ref:
        return "direct", True    # A1 is ALT
    if a1 == ref and a2 == alt:
        return "direct", False   # A1 is REF

    # Strand-flipped match
    c1, c2 = complement(a1), complement(a2)
    if c1 == alt and c2 == ref:
        return "flipped", True
    if c1 == ref and c2 == alt:
        return "flipped", False

    return None, None


# ---------------------------------------------------------------------------
# GWAS loading
# ---------------------------------------------------------------------------

# ---------------------------------------------------------------------------
# GWAS column name synonyms
# ---------------------------------------------------------------------------
# Maps canonical internal name -> set of known column name variants
# (case-insensitive matching is applied at runtime).
# Compiled from common GWAS summary statistic formats:
#   Bellenguez, LDSC munged, PLINK2, SAIGE, REGENIE, BOLT-LMM,
#   MRC-IEU harmonised (hm_*), FinnGen, and APOE-stratified VA/UKB formats.

GWAS_SYNONYMS = {
    "chrom": {
        "chr", "chrom", "chromosome", "hg19chrc", "hm_chrom",
    },
    "pos": {
        "bp", "pos", "position", "genpos",
        "base_pair_location", "base_pair_position", "hm_pos",
    },
    "allele1": {
        # Effect allele
        "allele1", "a1", "eff", "effect_allele", "effect_allele1",
        "hm_effect_allele", "allele0",
        # Named explicitly
        "effect_allele", "non_effect_allele".replace("non_", ""),  # belt-and-braces
        "ea",
    },
    "allele2": {
        # Non-effect / other allele
        "allele2", "a2", "noneff", "non_effect_allele", "other_allele",
        "hm_other_allele", "nea",
    },
    "effect": {
        # Beta / log-OR — all treated as log scale EXCEPT bare "or"/"OR"
        # (handled separately: those get math.log() applied)
        "beta", "effect", "hm_beta",
        # OR already on log scale (e.g. LogOR column)
        "logor",
        # OR on natural scale — flagged for conversion
        "or",
        # APOE-stratified columns (VA/UKB format)
        "phase1+2_apoe4pos_m1a_beta", "phase1+2_apoe4neg_m1a_beta",
    },
    "pvalue": {
        "p", "pval", "pvalue", "p_value", "p-value",
        "phase1+2_apoe4pos_m1a_p", "phase1+2_apoe4neg_m1a_p",
    },
    "rsid": {
        "snpid", "markername", "snp", "rsid", "id", "hm_rsid", "rsid", "snp_id",
    },
}

# Column names where the value is an odds ratio (not log-OR) and needs log()
OR_SCALE_NAMES = {"or"}


def detect_gwas_columns(header_cols, col_overrides=None):
    """
    Map canonical column names to their index in the header.

    header_cols: list of raw column name strings from the file header
    col_overrides: dict of canonical_name -> raw_column_name (from --col-* args)

    Returns:
        col_idx:    dict canonical_name -> int index
        effect_is_or: bool — True if effect column is on OR (not log) scale
        warnings:   list of warning strings

    Raises ValueError if a required canonical column cannot be mapped.
    """
    # Build a case-insensitive lookup: lowercase_raw_name -> index
    lower_to_idx = {col.strip().lower(): i for i, col in enumerate(header_cols)}
    required = {"chrom", "pos", "allele1", "allele2", "effect", "pvalue"}

    col_idx = {}
    effect_is_or = False
    warnings = []

    # Apply explicit overrides first
    if col_overrides:
        for canonical, raw in col_overrides.items():
            raw_lower = raw.strip().lower()
            if raw_lower in lower_to_idx:
                col_idx[canonical] = lower_to_idx[raw_lower]
            else:
                raise ValueError(
                    f"Column override --col-{canonical}={raw!r} not found in file. "
                    f"Available columns: {list(header_cols)}"
                )

    # Auto-detect remaining columns from synonym table
    for canonical, synonyms in GWAS_SYNONYMS.items():
        if canonical in col_idx:
            continue   # already set by override

        matches = []
        for syn in synonyms:
            if syn.lower() in lower_to_idx:
                matches.append((syn, lower_to_idx[syn.lower()]))

        if not matches:
            if canonical in required:
                raise ValueError(
                    f"Cannot find required column '{canonical}' in GWAS file.\n"
                    f"  Tried synonyms: {sorted(synonyms)}\n"
                    f"  Found columns:  {list(header_cols)}\n"
                    f"  Use --col-{canonical}=<name> to specify it explicitly."
                )
            # Optional column (rsid) — skip silently
            continue

        if len(matches) > 1:
            # Multiple synonyms matched — pick the first and warn
            warnings.append(
                f"  Warning: multiple synonyms matched for '{canonical}': "
                f"{[m[0] for m in matches]}. Using '{matches[0][0]}'. "
                f"Use --col-{canonical}=<name> to override."
            )

        matched_name, matched_idx = matches[0]
        col_idx[canonical] = matched_idx

        # Check if effect column is on OR scale
        if canonical == "effect" and matched_name.lower() in OR_SCALE_NAMES:
            effect_is_or = True
            warnings.append(
                f"  Note: effect column '{matched_name}' detected as odds ratio "
                f"(not log scale). Values will be log-transformed automatically."
            )

    return col_idx, effect_is_or, warnings


def get_liftover(from_build, to_build="hg38"):
    """
    Initialise a pyliftover LiftOver object for on-the-fly coordinate
    conversion. Chain files are downloaded and cached automatically by
    pyliftover on first use (~50MB per chain file).

    from_build: source genome build string, e.g. "hg19", "hg18", "GRCh37"
    to_build:   target build (default "hg38")

    Returns a LiftOver object, or raises ImportError if pyliftover is not
    installed, or RuntimeError if the chain file cannot be loaded.

    Build name normalisation:
        GRCh37 / grch37 -> hg19
        GRCh38 / grch38 -> hg38
    These are the two most common alternative names for hg19/hg38.
    """
    # Normalise common alternative build names
    _aliases = {
        "grch37": "hg19", "grch38": "hg38",
        "hg37":   "hg19",   # non-standard but seen in the wild
    }
    from_build = _aliases.get(from_build.lower(), from_build.lower())
    to_build   = _aliases.get(to_build.lower(),   to_build.lower())

    if from_build == to_build:
        return None   # no liftover needed

    try:
        from pyliftover import LiftOver
    except ImportError:
        raise ImportError(
            "pyliftover is required for --liftover but is not installed.\n"
            "Install it with: pip install pyliftover --break-system-packages"
        )

    print(f"  Initialising liftover: {from_build} -> {to_build} "
          f"(chain file will be downloaded if not cached)...", flush=True)
    try:
        lo = LiftOver(from_build, to_build)
    except Exception as e:
        raise RuntimeError(
            f"Failed to load liftover chain file ({from_build} -> {to_build}): {e}\n"
            f"Check your internet connection or provide a chain file path with "
            f"--liftover-chain."
        )
    print(f"  Liftover ready.", flush=True)
    return lo


def lift_position(lo, chrom, pos_1based):
    """
    Lift a 1-based genomic position to the target build using pyliftover.

    pyliftover uses 0-based coordinates internally.

    Returns (new_chrom, new_pos_1based) or None if the position cannot be
    lifted (unmapped region, deleted sequence, or one-to-many mapping).

    One-to-many mappings (len(result) > 1) are discarded conservatively —
    this affects a tiny fraction of positions in segmental duplications.
    """
    if lo is None:
        return chrom, pos_1based   # no liftover requested

    result = lo.convert_coordinate(chrom, pos_1based - 1)  # 0-based input

    if not result:
        return None   # unmapped

    if len(result) > 1:
        return None   # one-to-many: discard conservatively

    new_chrom, new_pos_0based, strand, score = result[0]

    # Normalise new chrom to chr prefix
    if not new_chrom.startswith("chr"):
        new_chrom = "chr" + new_chrom

    return new_chrom, new_pos_0based + 1   # back to 1-based


def load_gwas(gwas_file, p_threshold, col_overrides=None, liftover_build=None,
              liftover_chain=None):
    """
    Load genome-wide significant SNPs from the GWAS file.

    Columns are detected automatically from a synonym table covering common
    GWAS summary statistic formats (Bellenguez, LDSC, PLINK2, SAIGE, REGENIE,
    BOLT-LMM, MRC-IEU harmonised, FinnGen, APOE-stratified VA/UKB).

    col_overrides:   dict canonical_name -> raw_column_name (from --col-* args)
    liftover_build:  source genome build of the GWAS file, e.g. "hg19" or
                     "GRCh37". If provided, positions are lifted to hg38 on
                     the fly using pyliftover. Positions that fail to lift
                     (unmapped or one-to-many) are silently dropped.
    liftover_chain:  path to a local chain file (overrides auto-download).
                     If None and liftover_build is set, pyliftover downloads
                     the chain file and caches it in ~/.pyliftover.

    Effect values in OR columns (detected by column name "or") are
    automatically log-transformed to log-OR scale.

    Skips indels and rows where pvalue > p_threshold.

    Returns:
        dict: (chrom_with_prefix, pos_1based_hg38) -> {
            "allele1": str, "allele2": str,
            "effect": float, "pvalue": float, "rsid": str,
            "strand_ambiguous": bool
        }
    """
    print(f"Loading GWAS summary statistics from {gwas_file}...", flush=True)
    open_func = gzip.open if gwas_file.endswith(".gz") else open

    # Initialise liftover if a source build is specified
    if liftover_chain:
        try:
            from pyliftover import LiftOver
            lo = LiftOver(liftover_chain)
            print(f"  Liftover: using chain file {liftover_chain}", flush=True)
        except Exception as e:
            print(f"Error loading chain file {liftover_chain}: {e}", file=sys.stderr)
            sys.exit(1)
    elif liftover_build:
        lo = get_liftover(liftover_build, "hg38")
    else:
        lo = None

    hits = {}
    n_total = 0
    n_indel = 0
    n_or_zero = 0   # OR=0 rows skipped (log(0) undefined)
    n_lifted = 0    # positions successfully lifted
    n_lift_fail = 0 # positions that failed to lift
    delimiter = None
    col_idx = None
    effect_is_or = False

    with open_func(gwas_file, "rt") as f:
        for line in f:
            if not line.strip():
                continue
            stripped = line.rstrip("\n")

            # ---- Header line ----
            if col_idx is None:
                delimiter = "\t" if "\t" in stripped else None
                header_cols = [c.strip() for c in
                               (stripped.split(delimiter) if delimiter
                                else stripped.split())]

                try:
                    col_idx, effect_is_or, warnings = detect_gwas_columns(
                        header_cols, col_overrides
                    )
                except ValueError as e:
                    print(f"Error: {e}", file=sys.stderr)
                    sys.exit(1)

                for w in warnings:
                    print(w, flush=True)

                print(f"  Delimiter: {'tab' if delimiter else 'whitespace'}",
                      flush=True)
                print(f"  Column mapping:", flush=True)
                for canonical, idx in sorted(col_idx.items()):
                    print(f"    {canonical:<12} -> [{idx}] {header_cols[idx]}",
                          flush=True)
                if effect_is_or:
                    print(f"  Effect scale: odds ratio (will log-transform)",
                          flush=True)
                else:
                    print(f"  Effect scale: log scale (beta / log-OR)",
                          flush=True)
                continue

            # ---- Data line ----
            parts = stripped.split(delimiter) if delimiter else stripped.split()
            n_total += 1

            try:
                a1     = parts[col_idx["allele1"]].strip()
                a2     = parts[col_idx["allele2"]].strip()
                effect = float(parts[col_idx["effect"]])
                pvalue = float(parts[col_idx["pvalue"]])
                chrom  = str(parts[col_idx["chrom"]]).strip()
                pos    = int(float(parts[col_idx["pos"]]))
                rsid   = (parts[col_idx["rsid"]].strip()
                          if "rsid" in col_idx else "NA")
            except (ValueError, IndexError):
                continue

            # Lift position to hg38 if source build differs
            if lo is not None:
                lifted = lift_position(lo, "chr" + chrom if not chrom.startswith("chr") else chrom, pos)
                if lifted is None:
                    n_lift_fail += 1
                    continue
                chrom, pos = lifted
                n_lifted += 1

            # Convert OR to log-OR
            if effect_is_or:
                if effect <= 0:
                    n_or_zero += 1
                    continue
                effect = math.log(effect)

            # Skip indels
            if len(a1) != 1 or len(a2) != 1:
                n_indel += 1
                continue

            if pvalue > p_threshold:
                continue

            # Normalise chromosome to "chr" prefix
            if not chrom.startswith("chr"):
                chrom = "chr" + chrom

            key = (chrom, pos)
            if key in hits and hits[key]["pvalue"] <= pvalue:
                continue

            hits[key] = {
                "allele1":          a1.upper(),
                "allele2":          a2.upper(),
                "effect":           effect,
                "pvalue":           pvalue,
                "rsid":             rsid,
                "strand_ambiguous": is_strand_ambiguous(a1, a2),
            }

    print(f"  {n_total:,} total SNPs; {n_indel:,} indels skipped; "
          f"{len(hits):,} SNPs at P < {p_threshold:.2e}", flush=True)
    if lo is not None:
        print(f"  Liftover: {n_lifted:,} positions lifted successfully; "
              f"{n_lift_fail:,} failed (unmapped or one-to-many, discarded)",
              flush=True)
    if n_or_zero:
        print(f"  Warning: {n_or_zero:,} rows skipped (OR <= 0, log undefined)",
              flush=True)

    if hits:
        sample = list(hits.items())[:3]
        print("  Sample loaded hits (first 3):", flush=True)
        for (chrom, pos), v in sample:
            print(f"    {chrom}:{pos}  A1={v['allele1']}  A2={v['allele2']}  "
                  f"Effect={v['effect']:.4f}  P={v['pvalue']:.2e}", flush=True)


    return hits


# ---------------------------------------------------------------------------
# VCF loading — join GWAS hits to phase sets
# ---------------------------------------------------------------------------

def load_vcf_phase_assignments(vcf_file, gwas_hits):
    """
    Parse the phased VCF and, for each GWAS hit position, determine:
      - Which phase set (PS) the SNP belongs to
      - Which haplotype (1 or 2) carries the GWAS risk allele (Allele1)

    Only phased variants (GT uses '|') with a PS field are considered.
    Unphased variants (GT uses '/') are skipped.

    VCF phasing convention used by longcallR:
        0|1 -> REF on H1 (first allele = H1), ALT on H2
        1|0 -> ALT on H1, REF on H2

    Returns:
        phase_assignments: dict
            (chrom, pos) -> {
                "phase_set":      int,
                "risk_haplotype": int (1 or 2),  # which haplotype carries A1
                "match_type":     str,            # "direct" or "flipped"
                "strand_ambiguous": bool,
                **gwas_hit fields
            }
        phase_set_to_hits: dict
            phase_set_id -> list of (chrom, pos) keys
    """
    print(f"Reading phased VCF: {vcf_file}...", flush=True)
    open_func = gzip.open if vcf_file.endswith(".gz") else open

    # Index GWAS hits by position for O(1) lookup
    gwas_pos_set = set(gwas_hits.keys())

    phase_assignments = {}
    phase_set_to_hits = defaultdict(list)
    n_matched = 0
    n_no_match = 0
    n_unphased = 0

    with open_func(vcf_file, "rt") as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 10:
                continue

            chrom  = parts[0]
            pos    = int(parts[1])
            vcf_ref = parts[3]
            vcf_alt = parts[4]
            fmt    = parts[8].split(":")
            sample = parts[9].split(":")

            # Only consider SNPs (skip indels in VCF too)
            if len(vcf_ref) != 1 or len(vcf_alt) != 1:
                continue

            # Normalise VCF chromosome to chr-prefix to match GWAS hits,
            # which are always stored with chr prefix after load_gwas().
            # Handles VCFs aligned to references without chr prefix (e.g.
            # GRCh38 without UCSC-style naming, or hg19 integer chromosomes).
            if not chrom.startswith("chr"):
                chrom = "chr" + chrom

            # Check if this position is a GWAS hit
            key = (chrom, pos)
            if key not in gwas_pos_set:
                continue

            # Parse FORMAT fields into a dict
            fmt_dict = dict(zip(fmt, sample))
            gt = fmt_dict.get("GT", "")

            # Only phased variants (pipe-separated GT)
            if "|" not in gt:
                n_unphased += 1
                continue

            # Must have a phase set tag
            if "PS" not in fmt_dict:
                n_unphased += 1
                continue

            try:
                ps = int(fmt_dict["PS"])
            except ValueError:
                continue

            # Parse haplotype alleles: h1_allele_idx | h2_allele_idx
            gt_parts = gt.split("|")
            if len(gt_parts) != 2:
                continue
            try:
                h1_idx = int(gt_parts[0])   # 0 = REF, 1 = ALT
                h2_idx = int(gt_parts[1])
            except ValueError:
                continue

            # Skip homozygous (not informative for haplotype-specific analysis)
            if h1_idx == h2_idx:
                continue

            gwas = gwas_hits[key]
            a1, a2 = gwas["allele1"], gwas["allele2"]

            # Determine which haplotype carries the GWAS risk allele (A1)
            match_type, a1_is_alt = alleles_match(a1, a2, vcf_ref, vcf_alt)

            if match_type is None:
                # Alleles in GWAS don't match REF/ALT in VCF even after complement
                n_no_match += 1
                continue

            # a1_is_alt: True if GWAS A1 corresponds to VCF ALT (allele index 1)
            # h1_idx/h2_idx: 0=REF, 1=ALT
            if a1_is_alt:
                # Risk allele = ALT. Whichever haplotype has index 1 carries it.
                risk_hap = 1 if h1_idx == 1 else 2
            else:
                # Risk allele = REF. Whichever haplotype has index 0 carries it.
                risk_hap = 1 if h1_idx == 0 else 2

            n_matched += 1
            phase_assignments[key] = {
                "phase_set":        ps,
                "risk_haplotype":   risk_hap,
                "match_type":       match_type,
                "strand_ambiguous": gwas["strand_ambiguous"],
                "allele1":          a1,
                "allele2":          a2,
                "effect":           gwas["effect"],
                "pvalue":           gwas["pvalue"],
                "rsid":             gwas["rsid"],
                "chrom":            chrom,
                "pos":              pos,
            }
            phase_set_to_hits[ps].append(key)

    print(f"  {n_matched:,} GWAS hits matched to phased VCF positions; "
          f"{n_no_match:,} allele mismatches; {n_unphased:,} unphased",
          flush=True)

    # Diagnostic: show a few matched assignments so the user can verify
    if phase_assignments:
        sample = list(phase_assignments.items())[:5]
        print("  Sample VCF matches (first 5):", flush=True)
        for (chrom, pos), v in sample:
            print(f"    {chrom}:{pos}  PS={v['phase_set']}  "
                  f"risk_haplotype=H{v['risk_haplotype']}  "
                  f"A1={v['allele1']}  match={v['match_type']}", flush=True)
        print(f"  All matched phase sets: "
              f"{sorted(set(v['phase_set'] for v in phase_assignments.values()))}",
              flush=True)
    else:
        print("  WARNING: no GWAS hits matched to phased VCF positions.", flush=True)
        print("  This likely means either:", flush=True)
        print("    1. The GWAS SNPs are not in the VCF (not heterozygous in this sample)", flush=True)
        print("    2. Chromosome naming mismatch (check 'chr' prefix consistency)", flush=True)
        print("    3. The GWAS SNPs are present but unphased (GT uses '/' not '|')", flush=True)

    return phase_assignments, phase_set_to_hits


# ---------------------------------------------------------------------------
# longcallR result loading
# ---------------------------------------------------------------------------

def detect_mode_from_header(header_line):
    """
    Attempt to auto-detect the longcallR output mode from the header line.
    Returns "ase", "asj", or "asediting".

    Detection logic:
      - ASJ:       has "Hap1_absent" or "Hap1_present" columns
      - ASediting: has "Overall_editing_rate" column
      - ASE:       default (has H1/H2 read count columns)
    """
    h = header_line.lstrip("#").lower()
    if "hap1_absent" in h or "hap1_present" in h:
        return "asj"
    if "overall_editing" in h or "h1_editing_rate" in h:
        return "asediting"
    return "ase"


def load_longcallr_results(input_file):
    """
    Load a longcallR ASE, ASJ, or ASediting result file.

    All three formats share the key columns needed for this analysis:
        Phase_set      : integer phase set ID
        LogFC_H1_vs_H2 : log2 fold change (positive = H1 > H2)
        P_value        : per-site/gene/junction p-value
        P_value_adj    : naive BH-adjusted p-value

    Returns:
        header: list of column names
        rows:   list of dicts (one per result row)
        mode:   detected mode string
    """
    print(f"Loading longcallR results from {input_file}...", flush=True)
    open_func = gzip.open if input_file.endswith(".gz") else open

    header = None
    rows = []
    mode = None

    with open_func(input_file, "rt", encoding="utf-8-sig") as f:
        # utf-8-sig encoding automatically strips the UTF-8 BOM (\xef\xbb\xbf)
        # if present, which some tools write as the first bytes of a file and
        # which would otherwise prepend \ufeff to the first column name.
        for line in f:
            if not line.strip():
                continue
            if header is None:
                # Strip leading # comment marker and any leading/trailing whitespace
                # from the header line. Also detect delimiter: tab preferred,
                # fallback to whitespace (handles files with space-separated headers).
                raw_header = line.lstrip("#").rstrip("\n")
                delim = "\t" if "\t" in raw_header else None
                header = [col.strip() for col in
                          (raw_header.split(delim) if delim else raw_header.split())]
                mode = detect_mode_from_header(line)
                result_delim = delim
                continue
            parts = line.rstrip("\n").split(result_delim) if result_delim                     else line.rstrip("\n").split()
            if len(parts) != len(header):
                # Fallback: try both delimiters if column count doesn't match
                for fallback in ("\t", None):
                    test = line.rstrip("\n").split(fallback) if fallback                            else line.rstrip("\n").split()
                    if len(test) == len(header):
                        parts = test
                        break
            rows.append(dict(zip(header, parts)))

    print(f"  {len(rows):,} rows loaded (mode: {mode})", flush=True)
    return header, rows, mode


# ---------------------------------------------------------------------------
# Core annotation logic
# ---------------------------------------------------------------------------

def get_phase_set(row):
    """Extract integer phase set from a result row, handling various column names."""
    for col in ("Phase_set", "phase_set", "PS"):
        if col in row:
            try:
                return int(row[col])
            except (ValueError, TypeError):
                return None
    return None


def get_logfc(row):
    """
    Extract or compute LogFC_H1_vs_H2 from a result row.

    Handles column name variations across longcallR output modes:
      - ASE uses "logFC" (mixed case)
      - ASJ/ASediting use "LogFC_H1_vs_H2"

    For ASJ output, there is no precomputed LogFC column. Instead, the
    junction inclusion rate per haplotype is computed from the raw counts:
        H1_inclusion_rate = Hap1_present / (Hap1_present + Hap1_absent)
        H2_inclusion_rate = Hap2_present / (Hap2_present + Hap2_absent)
        logFC = log2(H1_inclusion_rate / H2_inclusion_rate)

    A pseudocount of 0.5 is added to numerator and 1 to denominator to
    avoid log(0) when one haplotype has zero junction-spanning reads.
    This is equivalent to a Laplace prior with alpha=0.5.
    """
    # Try direct logFC columns first (ASE and ASediting)
    for col in ("LogFC_H1_vs_H2", "logFC", "logfc", "LogFC", "log2fc"):
        if col in row and row[col] not in ("NA", "", "."):
            try:
                return float(row[col])
            except (ValueError, TypeError):
                pass

    # ASJ fallback: compute from junction inclusion counts
    # Columns: Hap1_absent, Hap1_present, Hap2_absent, Hap2_present
    try:
        h1_absent  = float(row["Hap1_absent"])
        h1_present = float(row["Hap1_present"])
        h2_absent  = float(row["Hap2_absent"])
        h2_present = float(row["Hap2_present"])
        # Laplace-smoothed inclusion rates
        h1_rate = (h1_present + 0.5) / (h1_present + h1_absent + 1.0)
        h2_rate = (h2_present + 0.5) / (h2_present + h2_absent + 1.0)
        import math
        return math.log2(h1_rate / h2_rate)
    except (KeyError, ValueError, TypeError, ZeroDivisionError):
        pass

    return None


def select_representative_snp(snp_keys, phase_assignments):
    """
    Given multiple GWAS SNPs in the same phase set, select the most
    significant (lowest p-value) as the representative SNP.

    Also checks that all SNPs agree on which haplotype carries the risk
    allele. If they disagree (e.g. due to complex haplotype structure or
    allele matching ambiguity), this is flagged.

    Returns (best_key, conflict_flag).
    """
    best_key = min(snp_keys, key=lambda k: phase_assignments[k]["pvalue"])
    best_hap = phase_assignments[best_key]["risk_haplotype"]

    # Check for haplotype conflicts among non-ambiguous SNPs
    non_ambiguous = [
        k for k in snp_keys
        if not phase_assignments[k]["strand_ambiguous"]
    ]
    conflict = False
    if len(non_ambiguous) > 1:
        haps = set(phase_assignments[k]["risk_haplotype"] for k in non_ambiguous)
        if len(haps) > 1:
            conflict = True

    return best_key, conflict


def annotate_results(rows, phase_set_to_hits, phase_assignments, keep_ambiguous):
    """
    For each longcallR result row, find GWAS hits in the same phase set and
    annotate with risk haplotype information.

    The LogFC in longcallR output is always H1 vs H2 (H1/H2).
    To get risk_logfc (risk_haplotype vs non-risk):
        if risk_haplotype == 1: risk_logfc = logfc        (already H1/H2)
        if risk_haplotype == 2: risk_logfc = -logfc       (flip to H2/H1)

    Returns list of annotated rows (dicts with extra annotation keys).
    """
    annotated = []
    n_matched = 0
    n_no_gwas = 0
    n_ambiguous_skipped = 0

    # Diagnostic: check overlap between GWAS phase sets and result phase sets
    gwas_phase_sets = set(phase_set_to_hits.keys())
    result_phase_sets = set()
    for row in rows:
        ps = get_phase_set(row)
        if ps is not None:
            result_phase_sets.add(ps)
    overlap = gwas_phase_sets & result_phase_sets
    print(f"  GWAS phase sets (from VCF): {sorted(gwas_phase_sets)}", flush=True)
    print(f"  Result phase sets (first 10): {sorted(result_phase_sets)[:10]}", flush=True)
    print(f"  Overlapping phase sets: {sorted(overlap)}", flush=True)

    for row in rows:
        ps = get_phase_set(row)
        logfc = get_logfc(row)

        if ps is None or logfc is None:
            annotated.append({**row, **_empty_annotation()})
            n_no_gwas += 1
            continue

        snp_keys = phase_set_to_hits.get(ps, [])
        if not snp_keys:
            annotated.append({**row, **_empty_annotation()})
            n_no_gwas += 1
            continue

        # Select representative SNP (most significant in phase set)
        best_key, conflict = select_representative_snp(snp_keys, phase_assignments)
        hit = phase_assignments[best_key]

        # Skip strand-ambiguous SNPs unless --keep-ambiguous is set
        if hit["strand_ambiguous"] and not keep_ambiguous:
            annotated.append({**row, **_empty_annotation(reason="strand_ambiguous")})
            n_ambiguous_skipped += 1
            continue

        risk_hap = hit["risk_haplotype"]

        # Compute risk_logfc: positive means risk haplotype has higher signal
        # LogFC_H1_vs_H2 = log2(H1/H2), so:
        #   risk_hap=1 -> risk_logfc = logfc         (H1 is risk, higher = more signal)
        #   risk_hap=2 -> risk_logfc = -logfc        (H2 is risk, flip sign)
        risk_logfc = logfc if risk_hap == 1 else -logfc

        snp_id = f"{hit['chrom']}:{hit['pos']}"
        annotation = {
            "gwas_snp":               snp_id,
            "gwas_rsid":              hit["rsid"],
            "gwas_allele1":           hit["allele1"],
            "gwas_allele2":           hit["allele2"],
            "gwas_effect":            f"{hit['effect']:.6f}",
            "gwas_pvalue":            f"{hit['pvalue']:.6e}",
            "risk_haplotype":         f"H{risk_hap}",
            "risk_logfc":             f"{risk_logfc:.4f}",
            "direction_concordant":   "TRUE" if risk_logfc > 0 else "FALSE",
            "strand_ambiguous":       "TRUE" if hit["strand_ambiguous"] else "FALSE",
            "allele_match_type":      hit["match_type"],
            "haplotype_conflict":     "TRUE" if conflict else "FALSE",
            "n_gwas_snps_in_ps":      str(len(snp_keys)),
        }
        annotated.append({**row, **annotation})
        n_matched += 1

    print(f"  {n_matched:,} results annotated with GWAS risk allele; "
          f"{n_no_gwas:,} with no GWAS hit in phase set; "
          f"{n_ambiguous_skipped:,} skipped (strand ambiguous)", flush=True)
    return annotated


def _empty_annotation(reason="no_gwas_hit"):
    """Return annotation columns filled with NA for unmatched rows."""
    return {
        "gwas_snp":             "NA",
        "gwas_rsid":            "NA",
        "gwas_allele1":         "NA",
        "gwas_allele2":         "NA",
        "gwas_effect":          "NA",
        "gwas_pvalue":          "NA",
        "risk_haplotype":       "NA",
        "risk_logfc":           "NA",
        "direction_concordant": "NA",
        "strand_ambiguous":     reason if reason == "strand_ambiguous" else "NA",
        "allele_match_type":    "NA",
        "haplotype_conflict":   "NA",
        "n_gwas_snps_in_ps":    "NA",
    }


ANNOTATION_COLS = [
    "gwas_snp", "gwas_rsid", "gwas_allele1", "gwas_allele2",
    "gwas_effect", "gwas_pvalue",
    "risk_haplotype", "risk_logfc", "direction_concordant",
    "strand_ambiguous", "allele_match_type", "haplotype_conflict",
    "n_gwas_snps_in_ps",
]


# ---------------------------------------------------------------------------
# Summary table: one row per GWAS hit x phase set
# ---------------------------------------------------------------------------

def write_summary_table(phase_set_to_hits, phase_assignments, rows, output_prefix):
    """
    Write a secondary summary table with one row per GWAS-hit x phase-set
    combination, showing all longcallR results using that phase set.

    This is useful for quickly seeing all genes/junctions associated with
    each GWAS locus.
    """
    out_path = output_prefix + ".gwas_hits_summary.tsv"
    print(f"Writing GWAS hit summary to {out_path}...", flush=True)

    # Build phase_set -> list of result rows
    ps_to_rows = defaultdict(list)
    for row in rows:
        ps = get_phase_set(row)
        if ps is not None:
            ps_to_rows[ps].append(row)

    with open(out_path, "w") as f:
        f.write(
            "#Phase_set\tGWAS_SNP\tGWAS_RSID\tGWAS_Allele1\tGWAS_Allele2\t"
            "GWAS_Effect\tGWAS_Pvalue\tRisk_haplotype\tStrand_ambiguous\t"
            "N_results_in_ps\tResult_names\tResult_logFCs\tResult_pvalues\n"
        )

        # Sort by GWAS p-value (most significant first)
        sorted_ps = sorted(
            phase_set_to_hits.keys(),
            key=lambda ps: min(
                phase_assignments[k]["pvalue"]
                for k in phase_set_to_hits[ps]
            )
        )

        for ps in sorted_ps:
            snp_keys = phase_set_to_hits[ps]
            best_key, _ = select_representative_snp(snp_keys, phase_assignments)
            hit = phase_assignments[best_key]

            ps_rows = ps_to_rows.get(ps, [])

            # Extract name column — try common longcallR name columns.
            # For ASJ, Junction is the primary identifier but Gene_name is
            # also present as the last column — include both for clarity.
            def get_name(r):
                junction = r.get("Junction", "")
                gene     = r.get("Gene_name", r.get("gene_name", ""))
                if junction and gene:
                    return f"{gene}:{junction}"
                if junction:
                    return junction
                if gene:
                    return gene
                for col in ("Chrom", "#Chrom"):
                    if col in r:
                        return r[col]
                return "?"

            names  = ";".join(get_name(r) for r in ps_rows)
            # Use get_logfc() to handle column name variations across modes
            # (ASE uses "logFC", ASJ/ASediting use "LogFC_H1_vs_H2")
            logfcs = ";".join(
                str(get_logfc(r)) if get_logfc(r) is not None else "NA"
                for r in ps_rows
            )
            # Use get_pvalue() helper to handle P_value vs P_value_adj variations
            def get_pval(r):
                for col in ("P_value_adj", "P_element_adj", "P_value", "P_value_adj"):
                    if col in r and r[col] not in ("NA", ""):
                        return r[col]
                return "NA"
            pvals  = ";".join(get_pval(r) for r in ps_rows)

            f.write(
                f"{ps}\t{hit['chrom']}:{hit['pos']}\t{hit['rsid']}\t"
                f"{hit['allele1']}\t{hit['allele2']}\t"
                f"{hit['effect']:.6f}\t{hit['pvalue']:.6e}\t"
                f"H{hit['risk_haplotype']}\t"
                f"{'TRUE' if hit['strand_ambiguous'] else 'FALSE'}\t"
                f"{len(ps_rows)}\t{names}\t{logfcs}\t{pvals}\n"
            )


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

def _parse_col_overrides(args):
    """Collect --col-* arguments into a dict for load_gwas."""
    overrides = {}
    for col in ("chrom", "pos", "allele1", "allele2", "effect", "pvalue", "rsid"):
        val = getattr(args, f"col_{col}", None)
        if val is not None:
            overrides[col] = val
    return overrides or None


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Annotate longcallR ASE/ASJ/ASediting results with GWAS "
                    "risk allele directionality."
    )
    parser.add_argument("--input", required=True,
                        help="longcallR result file (ASE, ASJ, or ASediting TSV)")
    parser.add_argument("--mode", choices=["ase", "asj", "asediting", "auto"],
                        default="auto",
                        help="Input file mode. Default: auto-detect from header")
    parser.add_argument("--vcf", required=True,
                        help="longcallR phased VCF file (plain or bgzipped)")
    parser.add_argument("--gwas", required=True,
                        help="GWAS summary statistics file (plain or gzipped). "
                             "Required columns: Allele1, Allele2, Effect, PVALUE, "
                             "chr, pos")
    parser.add_argument("--out", required=True,
                        help="Output file prefix")
    parser.add_argument("--gwas-p", type=float, default=5e-8,
                        help="GWAS p-value threshold for genome-wide significance. "
                             "Default: 5e-8")
    parser.add_argument("--liftover", default=None, metavar="BUILD",
                        help="Source genome build of the GWAS file, e.g. hg19 or "
                             "GRCh37. When provided, SNP positions are lifted to "
                             "hg38 on the fly using pyliftover. Chain files are "
                             "downloaded automatically and cached in ~/.pyliftover. "
                             "Positions that fail to lift are silently dropped. "
                             "Example: --liftover hg19")
    parser.add_argument("--liftover-chain", default=None, metavar="FILE",
                        help="Path to a local liftover chain file (e.g. "
                             "hg19ToHg38.over.chain.gz). Overrides auto-download. "
                             "Use when internet access is restricted on the HPC.")
    # Column name overrides — use when auto-detection fails or is ambiguous
    _col_group = parser.add_argument_group(
        "GWAS column overrides",
        "Explicitly specify column names for any field where auto-detection "
        "fails or produces a warning. These take precedence over the synonym "
        "table. Example: --col-chrom chromosome --col-effect BETA"
    )
    for _col in ("chrom", "pos", "allele1", "allele2", "effect", "pvalue", "rsid"):
        _col_group.add_argument(
            f"--col-{_col}", default=None, metavar="NAME",
            help=f"Column name for '{_col}' in the GWAS file"
        )
    parser.add_argument("--keep-ambiguous", action="store_true",
                        help="Include strand-ambiguous SNPs (A/T and C/G pairs) "
                             "in the annotation. By default these are excluded "
                             "since allele direction cannot be confirmed.")

    args = parser.parse_args()

    # 1. Load GWAS hits
    gwas_hits = load_gwas(args.gwas, args.gwas_p,
                         col_overrides=_parse_col_overrides(args),
                         liftover_build=args.liftover,
                         liftover_chain=args.liftover_chain)
    if not gwas_hits:
        print(f"Error: no GWAS SNPs found at P < {args.gwas_p:.2e}. "
              f"Check your GWAS file and threshold.", file=sys.stderr)
        sys.exit(1)

    # 2. Load VCF and map GWAS hits to phase sets
    phase_assignments, phase_set_to_hits = load_vcf_phase_assignments(
        args.vcf, gwas_hits
    )
    if not phase_assignments:
        print("Warning: no GWAS hits could be matched to phased VCF positions.\n"
              "Check that chromosome naming is consistent (chr1 vs 1) between\n"
              "the GWAS file and the VCF.", file=sys.stderr)

    # 3. Load longcallR results
    header, rows, mode = load_longcallr_results(args.input)
    if args.mode != "auto":
        mode = args.mode

    # 4. Annotate
    print(f"Annotating {len(rows):,} results...", flush=True)
    annotated = annotate_results(
        rows, phase_set_to_hits, phase_assignments, args.keep_ambiguous
    )

    # 5. Write annotated output
    out_path = args.out + f".{mode}_gwas_annotated.tsv"
    print(f"Writing annotated results to {out_path}...", flush=True)
    out_cols = header + ANNOTATION_COLS
    with open(out_path, "w") as f:
        f.write("#" + "\t".join(out_cols) + "\n")
        for row in annotated:
            f.write("\t".join(row.get(col, "NA") for col in out_cols) + "\n")

    # 6. Write GWAS hit summary table
    write_summary_table(phase_set_to_hits, phase_assignments, annotated, args.out)

    # 7. Print concordance summary for matched results
    matched = [r for r in annotated if r.get("risk_logfc", "NA") != "NA"]
    if matched:
        n_concordant = sum(1 for r in matched
                          if r.get("direction_concordant") == "TRUE")
        n_discord    = sum(1 for r in matched
                          if r.get("direction_concordant") == "FALSE")
        print(f"\nDirectionality summary ({len(matched)} results with GWAS annotation):",
              flush=True)
        print(f"  Concordant (risk allele -> higher signal): {n_concordant} "
              f"({100*n_concordant/len(matched):.1f}%)", flush=True)
        print(f"  Discordant (risk allele -> lower signal):  {n_discord} "
              f"({100*n_discord/len(matched):.1f}%)", flush=True)
        print(f"  (Enrichment of concordance over 50% suggests risk allele "
              f"increases expression/splicing/editing)", flush=True)

    print("Done.", flush=True)
