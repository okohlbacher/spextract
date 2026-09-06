#!/usr/bin/env python3
"""Truth-set metrics and vault reporting for the SpeXtractor harness.

Kept separate from bench.py because these run against EXISTING output and must be re-runnable
without recomputing a benchmark. Everything here obeys the same rule: no uncorrected number is
ever emitted alone (FAILURE 4 -- a 97.8% "detection ceiling" with a 91.9% chance floor had to
be publicly retracted).
"""
from __future__ import annotations
import bisect, collections, csv, json, sys, time
from pathlib import Path

import numpy as np
import pyarrow.parquet as pq
from pyteomics import mzml
import yaml

HERE = Path(__file__).resolve().parent
DECOY_SHIFT = 11.003          # Da; off the real m/z grid but inside the same isolation window
PPM_PREC = PPM_FRAG = 20.0
D_RT, D_IM = 0.30, 0.05


def rt_minutes(scan) -> float:
    """FAILURE 8: our mzML writes RT in SECONDS, the reference implementation's in MINUTES, DIA-NN reports
    minutes. Trusting the bare number matched 31 of ~300,000 spectra. Always read the CV unit."""
    v = scan["scan start time"]
    unit = getattr(v, "unit_info", "minute")
    return float(v) / (60.0 if unit == "second" else 1.0)


def precursor_of(spec) -> tuple[float, int, float]:
    """FAILURE 8b: the reference implementation omits 'selected ion m/z' and writes the PRECURSOR as its isolation
    window (+-0.01 Da); we write the real acquisition window. Reading either file's window
    naively makes the reference implementation look perfectly pure, because its co-isolation set is empty by
    construction. Fall back explicitly and never derive co-isolation from the file's window."""
    pre = spec["precursorList"]["precursor"][0]
    ion = pre["selectedIonList"]["selectedIon"][0]
    iso = pre.get("isolationWindow", {})
    mz = float(ion.get("selected ion m/z", iso.get("isolation window target m/z", 0)) or 0)
    z = int(ion.get("charge state", 0) or 0)
    sc = spec["scanList"]["scan"][0]
    im = float(ion.get("inverse reduced ion mobility",
                       sc.get("inverse reduced ion mobility", -1)) or -1)
    return mz, z, im


def load_truth(report: Path, lib: Path):
    r = pq.read_table(report, columns=["Precursor.Id", "Precursor.Mz", "Precursor.Charge",
                                       "RT", "RT.Start", "RT.Stop", "IM", "Q.Value"]).to_pandas()
    r = r[r["Q.Value"] <= 0.01].reset_index(drop=True)
    frags = collections.defaultdict(list)
    if lib and lib.exists():
        L = pq.read_table(lib, columns=["Precursor.Id", "Product.Mz", "Decoy"]).to_pandas()
        for pid, mz in zip(L[L["Decoy"] == 0]["Precursor.Id"], L[L["Decoy"] == 0]["Product.Mz"]):
            frags[pid].append(mz)
    return r, {k: np.sort(np.array(v)) for k, v in frags.items()}


def truth_metrics(mzml_path: Path, report: Path, lib: Path, windows: list[tuple[float, float]]) -> dict:
    """Recall, purity and charge agreement -- each with its own chance control.

    FAILURE 4: the co-isolation figure was ~half random matching until a SIZE-MATCHED decoy was
    added (the own-fragment decoy controls ~9 m/z, the co-isolated set ~450). Recall likewise
    needed a floor: 74.3% raw against a 35.6% floor is 60.1% corrected, not 74.3%.
    """
    r, frags = load_truth(report, lib)
    P = r.to_dict("records")
    idx = collections.defaultdict(list)
    RTB = 0.1
    for i, p in enumerate(P):
        for b in range(int((p["RT.Start"] - D_RT) / RTB), int((p["RT.Stop"] + D_RT) / RTB) + 1):
            idx[b].append(i)

    def hits(arr, mz, ppm=PPM_FRAG) -> bool:
        if len(arr) == 0:
            return False
        i = np.searchsorted(arr, mz)
        return any(0 <= j < len(arr) and abs(arr[j] - mz) / mz * 1e6 <= ppm for j in (i - 1, i))

    def window_of(mz):
        for lo, hi in windows:
            if lo <= mz <= hi:
                return (lo, hi)
        return None

    matched, matched_dec, present = set(), set(), set()
    zc = collections.Counter()
    own, oth, own_d, oth_d = [], [], [], []
    n_spec = 0

    for s in mzml.read(str(mzml_path)):
        if s.get("ms level") != 2:
            continue
        mz, z, im = precursor_of(s)
        if mz <= 0:
            continue
        rt = rt_minutes(s["scanList"]["scan"][0])
        w = window_of(mz)
        if w is None:
            continue
        n_spec += 1
        mzs, ints = s["m/z array"], s["intensity array"]
        if len(mzs) == 0:
            continue
        co = [P[i] for i in idx.get(int(rt / RTB), ())
              if w[0] <= P[i]["Precursor.Mz"] <= w[1] and (im < 0 or abs(P[i]["IM"] - im) <= D_IM)]
        dmz = mz + DECOY_SHIFT
        for p in co:
            if abs(p["Precursor.Mz"] - dmz) / dmz * 1e6 <= PPM_PREC:
                matched_dec.add(p["Precursor.Id"])
        cand = [p for p in co if abs(p["Precursor.Mz"] - mz) / mz * 1e6 <= PPM_PREC]
        if not cand:
            continue
        best = min(cand, key=lambda p: abs(p["Precursor.Mz"] - mz))
        matched.add(best["Precursor.Id"])
        if z:
            zc[(z, best["Precursor.Charge"])] += 1

        fo = frags.get(best["Precursor.Id"], np.array([]))
        ol = [frags[p["Precursor.Id"]] for p in co
              if p["Precursor.Id"] != best["Precursor.Id"] and p["Precursor.Id"] in frags]
        others = np.sort(np.concatenate(ol)) if ol else np.array([])
        k = np.argsort(ints)[-20:]                     # equal depth: the reference implementation emits 500, we ~228
        m2, i2 = mzs[k], ints[k]
        tot = i2.sum() or 1.0
        io = ib = iod = idc = 0.0
        for q, it in zip(m2, i2):
            if hits(fo, q):
                io += it
            elif hits(others, q):
                ib += it
            if hits(fo + DECOY_SHIFT, q):
                idc += it
            if len(others) and hits(others + DECOY_SHIFT, q):
                iod += it
        own.append(io / tot); oth.append(ib / tot)
        own_d.append(idc / tot); oth_d.append(iod / tot)

    for p in P:
        present.add(p["Precursor.Id"])
    o, b = np.array(own or [0]), np.array(oth or [0])
    od, bd = np.array(own_d or [0]), np.array(oth_d or [0])
    rec = len(matched & present) / max(len(present), 1)
    recd = len(matched_dec & present) / max(len(present), 1)
    ag = sum(v for (a, t), v in zc.items() if a == t)
    tot_z = sum(zc.values())
    # FAILURE 4 extension: charge agreement is meaningless without the majority-class rate --
    # always answering z=2 scores 69.6%, so 71.8% is +2.2 points, not a triumph.
    from collections import Counter
    truth_z = Counter(p["Precursor.Charge"] for p in P)
    majority = max(truth_z.values()) / max(sum(truth_z.values()), 1)
    return {
        "spectra_scored": n_spec,
        "recall": round(100 * rec, 1), "recall_decoy_floor": round(100 * recd, 1),
        "recall_corrected": round(100 * (rec - recd) / max(1 - recd, 1e-9), 1),
        "own_pct": round(100 * o.mean(), 2), "own_decoy": round(100 * od.mean(), 2),
        "own_net": round(100 * (o.mean() - od.mean()), 2),
        "coiso_pct": round(100 * b.mean(), 2), "coiso_decoy": round(100 * bd.mean(), 2),
        "coiso_net": round(100 * (b.mean() - bd.mean()), 2),
        "charge_agreement": round(100 * ag / max(tot_z, 1), 1),
        "charge_majority_baseline": round(100 * majority, 1),
        "charge_vs_majority": round(100 * ag / max(tot_z, 1) - 100 * majority, 1),
        "charge_confusions": sorted(((v, f"{a}->{t}") for (a, t), v in zc.items() if a != t),
                                    reverse=True)[:4],
    }


# the 24 diaPASEF acquisition tiles, recovered from DiaFrameMsMsWindows. FAILURE 8b: this is
# used INSTEAD of either file's recorded isolation window so both tools are scored identically.
WINDOWS_S23 = [(327.5,468.8),(467.8,500.6),(499.6,524.3),(523.3,545.8),(544.8,564.6),
               (563.6,582.9),(581.8,602.8),(601.8,622.4),(621.4,642.6),(641.6,663.3),
               (662.3,685.3),(684.3,708.3),(707.3,733.0),(732.0,759.9),(758.9,789.5),
               (788.5,822.6),(821.6,860.3),(859.3,903.9),(902.9,955.6),(954.6,1018.6),
               (1017.6,1098.9),(1097.9,1206.0),(1205.0,1364.6),(1363.6,1650.0)]


def cmd_truth(args) -> None:
    reg = yaml.safe_load((HERE / "samples.yaml").read_text())
    plan = yaml.safe_load(Path(args.plan).read_text())
    root = Path(args.out or f"/path/to/scratch/harness/{Path(args.plan).stem}")
    for d in sorted(root.glob("*__*")):
        man = json.loads((d / "manifest.json").read_text())
        entry = reg["samples"][man["sample"]]
        rep = Path(entry["references"]["diann_report"])
        lib = Path(entry["references"].get("diann_lib", ""))
        if not rep.exists():
            print(f"[truth] {d.name}: no DIA-NN reference for {man['sample']} -- skipped")
            continue
        m = truth_metrics(d / "pseudo.mzML", rep, lib, WINDOWS_S23)
        (d / "truth.json").write_text(json.dumps(m, indent=2))
        print(f"[truth] {d.name}: recall {m['recall_corrected']}% (raw {m['recall']}, "
              f"floor {m['recall_decoy_floor']})  charge {m['charge_agreement']}% "
              f"(majority {m['charge_majority_baseline']}%, {m['charge_vs_majority']:+})")


def cmd_collate(args) -> None:
    """Write one vault report per plan. Every row carries its sample, and no ratio is computed
    across samples -- FAILURE 1 was a published number that divided dataset D by dataset A's reference."""
    reg = yaml.safe_load((HERE / "samples.yaml").read_text())
    plan = yaml.safe_load(Path(args.plan).read_text())
    root = Path(args.out or f"/path/to/scratch/harness/{Path(args.plan).stem}")
    rows = []
    for d in sorted(root.glob("*__*")):
        man = json.loads((d / "manifest.json").read_text())
        met = json.loads((d / "metrics.json").read_text()) if (d / "metrics.json").exists() else {}
        tru = json.loads((d / "truth.json").read_text()) if (d / "truth.json").exists() else {}
        rows.append({**man, **met, **tru})

    vault = Path(reg["vault"]); vault.mkdir(parents=True, exist_ok=True)
    stamp = time.strftime("%Y-%m-%d")
    out = vault / f"{Path(args.plan).stem}.md"
    L = [f"# {plan['title']}", "",
         f"Generated by `harness/bench.py` v{plan.get('harness_version','1.0')} on {stamp}. "
         "Do not hand-edit; re-run `collate` instead.", "",
         "All peptide counts are **`peptide_q <= 0.01`**. `spectrum_q` is shown only to expose "
         "`fdr_loss`, which is itself the redundancy signal.", ""]

    by_sample = collections.defaultdict(list)
    for r in rows:
        by_sample[r["sample"]].append(r)

    for sample, rs in sorted(by_sample.items()):
        L += [f"## {sample}", "",
              "| arm | peptide_q | spectrum_q | fdr_loss | recall (corr) | charge | vs majority | wall | RSS |",
              "|---|---|---|---|---|---|---|---|---|"]
        for r in sorted(rs, key=lambda x: -(x.get("peptide_q_peptides") or 0)):
            L.append("| {arm} | **{pq}** | {sq} | {fl}% | {rc} | {ch} | {cm} | {w}s | {m} MB |".format(
                arm=r["arm"], pq=r.get("peptide_q_peptides", "-"),
                sq=r.get("spectrum_q_peptides", "-"), fl=r.get("fdr_loss_pct", "-"),
                rc=f"{r['recall_corrected']}%" if "recall_corrected" in r else "-",
                ch=f"{r['charge_agreement']}%" if "charge_agreement" in r else "-",
                cm=f"{r['charge_vs_majority']:+}" if "charge_vs_majority" in r else "-",
                w=r.get("wall_s", "-"), m=r.get("peak_rss_mb", "-")))
        L.append("")
        b = rs[0]
        L += [f"Binary `{b['binary']['sha256']}` built {b['binary']['mtime']}; "
              f"git `{b.get('git_ref') or 'n/a'}`; host `{b['host']}`.", ""]

    L += ["## Provenance", "",
          "Each row's full config, binary hash and sample fingerprint are in "
          "`manifest.json` beside its outputs. The sample was verified by reading its "
          "fingerprint back out of the produced mzML, not taken from the invocation.", ""]
    out.write_text("\n".join(L))
    print(f"[collate] wrote {out}")
