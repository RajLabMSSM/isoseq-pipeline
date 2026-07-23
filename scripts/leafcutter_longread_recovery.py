#!/usr/bin/env python3
"""leafcutter short-read cluster recovery by v5 long-read isoforms, PER ASSEMBLER.

For every intron leafcutter tests in the SHORT reads, does each v5 assembler's
read-supported transcript models contain that junction? Reported per tool so you can
see, e.g., which cryptic events only one tool recovers. This is the v5 re-run of the
earlier tdpkd_leafcutter_longread_recovery analysis, generalized to any cohort.

Short-read side = a leafcutter-ds differential-splicing run (effect_sizes +
cluster_significance). Every tested intron is a row; `is_significant` marks the
differentially-spliced ones (cluster p.adjust < 0.05 AND |deltapsi| >= 0.1), so both the
full tested set and the DS set are covered.

Long-read / assembler side (v5 — one consolidated run, no separate long-only StringTie):
  recover_bambu      Bambu *_extended_annotations.gtf restricted to read-supported
                     transcripts (count > 0 in *_counts_transcript.txt). PURE long-read.
  recover_isoquant   IsoQuant *.transcript_models.gtf (de-novo models; hybrid --illumina).
  recover_stringtie  StringTie *.stringtie.gtf (--mix hybrid: short-read BAM + long-read).
  recover_any        bambu OR isoquant OR stringtie.
Only bambu is pure long-read in v5; isoquant (--illumina) and stringtie (--mix) both see
the short reads, so recover_stringtie/isoquant are NOT pure-long-read evidence.

Matching: (chrom, start, end), strand-agnostic (ONT cDNA strand can flip), +/- WOBBLE bp,
implemented as a dilated query (each leafcutter junction probes its 25 neighbour keys
against the exact per-tool intron sets). novel = junction not in GENCODE (+/- WOBBLE).

Brooke Friedman
"""
import argparse
import gzip
import sys

WOBBLE = 2
PADJ_SIG = 0.05
DPSI_SIG = 0.1
TOOLS = ["bambu", "isoquant", "stringtie", "any"]


def opentext(path):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path)


def attr(field9, key):
    i = field9.find(key + ' "')
    if i < 0:
        return None
    i += len(key) + 2
    return field9[i:field9.find('"', i)]


def gtf_introns(path, keep):
    tx = {}
    with opentext(path) as fh:
        for ln in fh:
            if ln.startswith("#"):
                continue
            f = ln.split("\t")
            if len(f) < 9 or f[2] != "exon":
                continue
            tid = attr(f[8], "transcript_id")
            if keep(tid):
                tx.setdefault(tid, (f[0], []))[1].append((int(f[3]), int(f[4])))
    for chrom, exons in tx.values():
        if len(exons) >= 2:
            exons.sort()
            for (_, e1e), (e2s, _) in zip(exons, exons[1:]):
                yield (chrom, e1e + 1, e2s - 1)


def build_index(specs):
    idx = set()
    for path, keep in specs:
        try:
            for t in gtf_introns(path, keep):
                idx.add(t)
        except FileNotFoundError:
            sys.exit(f"missing GTF: {path}")
    return idx


def bambu_supported(counts_path):
    sup = set()
    with open(counts_path) as fh:
        next(fh)
        for ln in fh:
            c = ln.rstrip("\n").split("\t")
            if float(c[2]) > 0:
                sup.add(c[0])
    return sup


def hit(chrom, s, e, idx):
    for ds in range(-WOBBLE, WOBBLE + 1):
        for de in range(-WOBBLE, WOBBLE + 1):
            if (chrom, s + ds, e + de) in idx:
                return True
    return False


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--cohort", required=True)
    ap.add_argument("--data-code", required=True, help="e.g. tdpkd_nanopore")
    ap.add_argument("--v5-dir", required=True, help="results/v5 folder")
    ap.add_argument("--groups", required=True, help="comma-separated group names, e.g. TDP,CTRL")
    ap.add_argument("--leafcutter-sig", required=True)
    ap.add_argument("--leafcutter-eff", required=True)
    ap.add_argument("--gencode", required=True)
    ap.add_argument("-o", "--output", required=True)
    ap.add_argument("--summary", required=True)
    args = ap.parse_args()

    dc = args.data_code
    v5 = args.v5_dir.rstrip("/")
    groups = args.groups.split(",")

    print("indexing GENCODE introns ...", file=sys.stderr)
    gencode = build_index([(args.gencode, lambda t: True)])
    print(f"  GENCODE introns: {len(gencode):,}", file=sys.stderr)

    stringtie_specs = [(f"{v5}/stringtie3/{g}/{dc}_{g}.stringtie.gtf", lambda t: True) for g in groups]
    isoquant_specs = [(f"{v5}/isoquant/{g}/{dc}_{g}/{dc}_{g}.transcript_models.gtf", lambda t: True) for g in groups]
    bambu_specs = []
    for g in groups:
        d = f"{v5}/bambu/{g}"
        sup = bambu_supported(f"{d}/{dc}_{g}_counts_transcript.txt")
        bambu_specs.append((f"{d}/{dc}_{g}_extended_annotations.gtf",
                            (lambda s: (lambda t: t in s))(sup)))

    print("indexing per-tool assembler introns ...", file=sys.stderr)
    idx = {"bambu": build_index(bambu_specs),
           "isoquant": build_index(isoquant_specs),
           "stringtie": build_index(stringtie_specs)}
    for t in ("bambu", "isoquant", "stringtie"):
        print(f"  {t}: {len(idx[t]):,} introns", file=sys.stderr)

    clu = {}
    with open(args.leafcutter_sig) as fh:
        next(fh)
        for ln in fh:
            c = ln.rstrip("\n").split("\t")
            padj = c[5]
            clu[c[0]] = (float(padj) if padj not in ("NA", "") else None, c[6])

    rows = []
    with open(args.leafcutter_eff) as fh:
        next(fh)
        for ln in fh:
            c = ln.rstrip("\n").split("\t")
            p = c[0].split(":")
            chrom, start, end, clupart = p[0], int(p[1]), int(p[2]), p[3]
            strand = clupart.rsplit("_", 1)[1]
            cluster_id = f"{chrom}:{clupart}"
            deltapsi = float(c[4])
            padj, gene = clu.get(cluster_id, (None, "NA"))

            rec = {t: (1 if hit(chrom, start, end, idx[t]) else 0) for t in idx}
            rec["any"] = 1 if (rec["bambu"] or rec["isoquant"] or rec["stringtie"]) else 0
            is_novel = 0 if hit(chrom, start, end, gencode) else 1
            is_sig = 1 if (padj is not None and padj < PADJ_SIG and abs(deltapsi) >= DPSI_SIG) else 0

            rows.append({"junction_id": f"{chrom}:{start}-{end}:{strand}", "chrom": chrom,
                         "start": start, "end": end, "strand": strand, "cluster": cluster_id,
                         "gene": gene, "cluster_p_adjust": "" if padj is None else f"{padj:.6g}",
                         "deltapsi": f"{deltapsi:.6g}", "is_significant": is_sig,
                         "is_novel": is_novel, **{f"recover_{t}": rec[t] for t in TOOLS}})

    header = (["junction_id", "chrom", "start", "end", "strand", "cluster", "gene",
               "cluster_p_adjust", "deltapsi", "is_significant", "is_novel"]
              + [f"recover_{t}" for t in TOOLS])
    with open(args.output, "w") as out:
        out.write("\t".join(header) + "\n")
        for r in rows:
            out.write("\t".join(str(r[h]) for h in header) + "\n")

    sets = [
        ("RAW_all", rows),
        ("DS_significant", [r for r in rows if r["is_significant"]]),
        ("RAW_novel", [r for r in rows if r["is_novel"]]),
        ("DS_novel", [r for r in rows if r["is_significant"] and r["is_novel"]]),
    ]
    with open(args.summary, "w") as s:
        s.write("cohort\tjunction_set\tn_junctions\t" + "\t".join(f"{t}_n\t{t}_pct" for t in TOOLS) + "\n")
        for label, subset in sets:
            n = len(subset)
            cells = []
            for t in TOOLS:
                k = sum(r[f"recover_{t}"] for r in subset)
                cells += [str(k), f"{100*k/n:.1f}" if n else ""]
            s.write(f"{args.cohort}\t{label}\t{n}\t" + "\t".join(cells) + "\n")

    print(f"\nwrote {args.output}\nwrote {args.summary}", file=sys.stderr)
    for label, subset in sets:
        n = len(subset)
        print(f"\n{label}: {n:,} junctions", file=sys.stderr)
        for t in TOOLS:
            k = sum(r[f"recover_{t}"] for r in subset)
            print(f"  {t:12s} {k:6,}  ({100*k/n:5.1f}%)" if n else f"  {t}: 0", file=sys.stderr)


if __name__ == "__main__":
    main()
