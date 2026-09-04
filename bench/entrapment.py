#!/usr/bin/env python3
"""Entrapment FDR across arms -- the measurement this project has never made.

WHY THIS EXISTS. Every "1% FDR" number quoted in this project, ours and the reference implementation's, is NOMINAL:
it is Sage's target-decoy estimate computed from ITS OWN decoys, per run. All three adversarial
reviewers (codex, vibe, kimi) independently made the same point, and `bench/fdrcheck.py` then
confirmed it on data already on disk: the min_correlation sweep's "+7.7% peptides" was the q<=0.01
SCORE THRESHOLD sliding down as decoys were deleted, not new signal. A nominal FDR cannot detect
that, because the thing that moved IS the nominal FDR's own calibration.

Entrapment measures the ERROR RATE DIRECTLY. A foreign proteome is concatenated to the human
target database; the sample is human liver, so ANY peptide identified from the foreign proteome is
a known false positive, counted rather than modelled. This is the accepted validation when the
search space or the scoring changes (Wen, Freestone, Noble, Kall, Nesvizhskii, Nat Methods 2025).

WHAT IT SETTLES. Run over the min_correlation arms it answers, without ambiguity:
  * does the observed error rate at nominal 1% RISE as the gate tightens?  -> the gain was
    recalibration, and the nominal axis everyone optimised is not 1%
  * does it stay flat while peptides rise?                                 -> the gain is real
It equally applies to the the reference implementation arm, which bounds how much of the 3,036-peptide "coverage gap"
is the reference implementation's own FDR edge rather than peptides we genuinely lose.

ESTIMATOR. A false match distributes over the PEPTIDE search space, so the correction ratio must be
the ratio of entrapment to target PEPTIDE HYPOTHESES, not proteins (fixed 2026-07-28, codex review
#12: Arabidopsis vs human proteins differ in length/composition, so the protein ratio 0.800 is the
wrong number -- the in-silico-tryptic peptide ratio is 0.6805, giving a 1/r correction of 1.469 vs
the protein 1.250, i.e. the old estimate was ~18%% too low). Estimated false TARGET IDs = n_entrap /
r_pep (each foreign hit implies n_target_space/n_entrap_space false target hits); FDR = that / n_target.
A BOOTSTRAP 95%% CI is printed alongside every point estimate because n_entrap ~ 50-200 makes the
entrapment fraction a noisy statistic -- a point FDR without an error bar cannot separate a real
calibration change from resampling noise. Both the RAW fraction and the corrected FDR (with CI) are
printed, never the corrected number alone, so the reader can check the correction's assumptions.
"""
import argparse, csv, os, subprocess, sys, json, random
from pathlib import Path

_AA = {'G':57.02146,'A':71.03711,'S':87.03203,'P':97.05276,'V':99.06841,'T':101.04768,'C':103.00919,
       'L':113.08406,'I':113.08406,'N':114.04293,'D':115.02694,'Q':128.05858,'K':128.09496,'E':129.04259,
       'M':131.04049,'H':137.05891,'F':147.06841,'R':156.10111,'Y':163.06333,'W':186.07931}

def _digest(seq, mc=2, lo=7, hi=35, mlo=500.0, mhi=5000.0):
    """Tryptic peptides matching Sage (cleave after K/R not before P, <=mc missed, len+mass bounds)."""
    cuts = [0] + [i+1 for i in range(len(seq)-1) if seq[i] in 'KR' and seq[i+1] != 'P'] + [len(seq)]
    peps, n = set(), len(cuts)-1
    for i in range(n):
        for j in range(i+1, min(i+2+mc, n+1)):
            p = seq[cuts[i]:cuts[j]]
            if lo <= len(p) <= hi:
                m = 18.010565
                ok = True
                for a in p:
                    if a not in _AA: ok = False; break
                    m += _AA[a]
                if ok and mlo <= m <= mhi: peps.add(p)
    return peps

def peptide_ratio(fasta, tag="ENTRAP_"):
    """entrapment/target PEPTIDE-hypothesis ratio: |peptides UNIQUE to foreign| / |target peptides|."""
    tgt, ent, seq, is_ent = set(), set(), [], False
    def flush():
        if seq: (ent if is_ent else tgt).update(_digest("".join(seq)))
    with open(fasta) as f:
        for line in f:
            if line.startswith(">"):
                flush(); seq = []; is_ent = line[1:].startswith(tag)
            else: seq.append(line.strip())
        flush()
    return len(ent - tgt) / max(len(tgt), 1)

def _bootstrap_fdr(flags, ratio, B=1000, seed=1):
    """95% CI on the corrected target-FDR by resampling the accepted peptides (flags: 1=entrap)."""
    rng = random.Random(seed); n = len(flags)
    if n == 0: return (float("nan"), float("nan"))
    fdrs = []
    for _ in range(B):
        ne = 0; sample = (flags[rng.randrange(n)] for _ in range(n))
        ne = sum(sample)
        nt = n - ne
        fdrs.append(100.0 * (ne / ratio) / max(nt, 1) if ratio > 0 else float("nan"))
    fdrs.sort()
    return (fdrs[int(0.025*B)], fdrs[int(0.975*B)])


def build_entrapment(target_fa, foreign_fa, out_fa, tag="ENTRAP_"):
    """Concatenate a foreign proteome, prefixing every foreign accession."""
    n_t = n_e = 0
    with open(out_fa, "w") as out:
        with open(target_fa) as f:
            for line in f:
                if line.startswith(">"):
                    n_t += 1
                out.write(line)
        with open(foreign_fa) as f:
            for line in f:
                if line.startswith(">"):
                    n_e += 1
                    out.write(">" + tag + line[1:])
                else:
                    out.write(line)
    return n_t, n_e


def score(tsv, tag, ratio, q_field="peptide_q", q_cut=0.01):
    """Observed entrapment rate among peptides accepted at NOMINAL q_cut."""
    best = {}
    with open(tsv) as f:
        for x in csv.DictReader(f, delimiter="\t"):
            try:
                s = float(x["sage_discriminant_score"])
                q = float(x[q_field])
                prots = [p for p in x["proteins"].split(";") if p]
            except (ValueError, KeyError):
                continue
            if not prots or all(p.startswith("rev_") for p in prots):
                continue                      # decoys are Sage's model, not our measurement
            pep = x["peptide"]
            if pep not in best or s > best[pep][0]:
                # entrapment iff EVERY protein is foreign; a peptide shared with a human protein
                # is not evidence of error (shared tryptic peptides exist between proteomes)
                ent = all(p.startswith(tag) or p.startswith("rev_" + tag) for p in prots)
                best[pep] = (s, q, ent)
    acc = [(s, q, e) for s, q, e in best.values() if q <= q_cut]
    n_e = sum(1 for _, _, e in acc if e)
    n_t = len(acc) - n_e
    raw = 100.0 * n_e / max(len(acc), 1)
    # corrected TARGET-list FDR: est. false targets = n_entrap / r_pep, over accepted targets.
    fdr = 100.0 * (n_e / ratio) / max(n_t, 1) if ratio > 0 else float("nan")
    lo, hi = _bootstrap_fdr([1 if e else 0 for _, _, e in acc], ratio)
    return len(acc), n_t, n_e, raw, fdr, lo, hi


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--target", help="human FASTA (targets only, no decoys)")
    ap.add_argument("--foreign", help="foreign proteome FASTA (entrapment)")
    ap.add_argument("--prebuilt", help="an ALREADY-concatenated target+entrapment FASTA whose "
                    "foreign entries carry --tag. Counts are read back OUT of the file rather "
                    "than assumed, so a mislabelled database cannot silently set the ratio.")
    ap.add_argument("--sage", required=True, help="sage binary")
    ap.add_argument("--config", required=True, help="sage config json (its fasta is overridden)")
    ap.add_argument("--out", required=True, help="working dir for the entrapment FASTA + results")
    ap.add_argument("--arm", action="append", required=True, metavar="NAME=DIR",
                    help="arm name and the directory holding its pseudo.mzML (repeatable)")
    ap.add_argument("--tag", default="ENTRAP_")
    ap.add_argument("--q-cut", type=float, default=0.01)
    a = ap.parse_args()

    out = Path(a.out); out.mkdir(parents=True, exist_ok=True)
    if a.prebuilt:
        fa = Path(a.prebuilt)
        n_e = sum(1 for l in open(fa) if l.startswith(">" + a.tag))
        n_t = sum(1 for l in open(fa) if l.startswith(">")) - n_e
        if n_e == 0:
            sys.exit("ABORT: %s contains no '>%s' entries -- wrong tag or wrong file. Refusing to "
                     "report an entrapment rate against a database with no entrapment." % (fa, a.tag))
        print("[db] prebuilt %s: %d target + %d entrapment proteins" % (fa, n_t, n_e), flush=True)
    else:
        if not (a.target and a.foreign):
            sys.exit("ABORT: need --prebuilt, or both --target and --foreign")
        fa = out / "target_entrap.fasta"
        if not fa.exists():
            n_t, n_e = build_entrapment(a.target, a.foreign, fa, a.tag)
            print("[db] %d target + %d entrapment proteins -> %s" % (n_t, n_e, fa), flush=True)
        else:
            n_t = sum(1 for l in open(a.target) if l.startswith(">"))
            n_e = sum(1 for l in open(a.foreign) if l.startswith(">"))
            print("[db] reusing %s (%d target + %d entrapment)" % (fa, n_t, n_e), flush=True)
    prot_ratio = n_e / max(n_t, 1)
    ratio = peptide_ratio(fa, a.tag)   # PEPTIDE-hypothesis ratio drives the FDR correction (codex #12)
    print("[db] protein ratio %.3f (reference only); PEPTIDE-hypothesis ratio %.4f (drives FDR)"
          % (prot_ratio, ratio), flush=True)

    cfg = json.load(open(a.config))
    cfg["database"]["fasta"] = str(fa)
    cfg_path = out / "sage_entrap.json"
    json.dump(cfg, open(cfg_path, "w"), indent=2)

    rows = []
    for spec in a.arm:
        name, d = spec.split("=", 1)
        mzml = Path(d) / "pseudo.mzML"
        if not mzml.exists():
            # the reference implementation and other reference arms are a bare mzML path, not a run directory
            mzml = Path(d)
        if not mzml.exists():
            print("[skip] %s: no mzML at %s" % (name, d)); continue
        odir = out / name; odir.mkdir(exist_ok=True)
        tsv = odir / "results.sage.tsv"
        if not tsv.exists():
            print("[run] %s ..." % name, flush=True)
            subprocess.run([a.sage, str(cfg_path), "-o", str(odir),
                            "--disable-telemetry-i-dont-want-to-improve-sage", str(mzml)],
                           stdout=open(odir / "sage.log", "w"), stderr=subprocess.STDOUT)
        if not tsv.exists():
            print("[fail] %s: sage produced no tsv" % name); continue
        rows.append((name,) + score(str(tsv), a.tag, ratio, q_cut=a.q_cut))

    print("\n=== ENTRAPMENT FDR at nominal peptide_q <= %.3f (peptide-ratio corrected + 95%% bootstrap CI) ===" % a.q_cut)
    print("%-10s %10s %10s %8s %11s %10s %16s" % ("arm", "accepted", "target", "entrap",
                                                  "raw entrap%", "est. FDR%", "FDR 95% CI"))
    for n, tot, nt, ne, raw, fdr, lo, hi in rows:
        print("%-10s %10d %10d %8d %10.2f%% %9.2f%% %7.2f-%.2f%%" % (n, tot, nt, ne, raw, fdr, lo, hi))
    print("\nREAD: nominal FDR is 1%% BY CONSTRUCTION in every row -- that is the point. If est. FDR%%")
    print("      RISES as the gate tightens, the extra peptides are errors and the 'gain' was")
    print("      recalibration. If it stays flat while accepted counts rise, the gain is real.")


if __name__ == "__main__":
    sys.exit(main())
