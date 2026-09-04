# Real DDA ground truth (PXD027359 matched pair)

## Why this is different from every other number in this repo

Every previous benchmark is anchored on a DIA-NN library that is tryptic / z2–4 / ≤2 mods, and
therefore **structurally cannot contain the non-canonical species this tool exists to find**.
PXD027359 provides the same HeLa digest on the same 6-min Evosep gradient acquired **twice** —
once diaPASEF, once ddaPASEF. Real DDA has no library restriction.

One engine (MSFragger, which reads Bruker `.d` natively so the real-DDA side needs no
conversion), one shared target+decoy FASTA, one FDR procedure. The only difference is where the
spectra came from.

| | peptides @1% FDR |
|---|---|
| **real DDA** | **3,067** |
| **our pseudo-DDA** (from the matched DIA run) | **2,065** |
| both | 1,477 |
| real DDA only — we miss | **1,590 (51.8% of DDA)** |
| ours only | **588 (28.5% of ours)** |

**Recall of real-DDA peptides: 48.2%.**

## Reading this honestly

**We miss half of what real DDA finds.** That is a harder number than the 86.9% the reference implementation ratio,
and it is measured against a reference with no structural blind spot. It should be the headline
figure for "how good are these pseudo-spectra", not the DIA-NN-anchored ones.

**The 588 "ours only" peptides are NOT automatically a win.** DDA targets a limited number of
precursors per cycle while DIA fragments everything in the window, so some of the 588 are genuine
DIA depth — peptides DDA never selected — and some are false positives. **The direction cannot be
read off this number**, and characterising that set is the obvious next analysis.

## Caveats that bound the comparison

* **Different acquisitions, not different processing of one acquisition.** DDA and DIA sample the
  peptide population differently by design, so neither set is a superset of the other. 48.2% is
  not a defect rate — part of the 51.8% is DDA targeting precursors our DIA cycle never isolated
  cleanly, and part is our extraction failing on precursors it did.
* **6-min gradient, 361 MB file.** Both counts are small (3,067 / 2,065) relative to the 31-min
  samples used elsewhere. Whether recall holds at 10,000-peptide scale is untested.
* Single sample, single replicate. No confidence interval.

## Operational notes

* MSFragger cannot read Bruker `.d` without the Bruker native libraries on `LD_LIBRARY_PATH`; it
  reports `Scans = 0` and calls the file corrupt. Same JNA trap that made the reference implementation look like a
  12-hour hang. The libraries ship with the reference implementation (`ext/bruker`).
* MSFragger also cannot read OpenMS mzML without `fix_mzml_cvparams.py` — OpenMS writes valueless
  cvParams (legal mzML) and batmass-io calls `.trim()` on the null.
