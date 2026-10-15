"""Compact record of the exploratory (diagnostic) MOPAC campaigns cited in Sections 2.3-2.4 of the manuscript.

Three exploratory campaigns, run before the production protocol was fixed, document three statements of the
Methods section:
  v1  Sparkle/PM7 as distributed, COSMO (EPS=78.4), Cartesian coordinates (no O-H constraint)
      -> "COSMO ... gave artificially short Eu-O contacts and optimizer failures"
  v2  Sparkle/PM7 as distributed, gas phase, O-H bonds frozen
      -> "hydroxyl and alkoxide-like oxygens collapsed onto the metal at 1.3-1.45 A"
  v5  Sparkle/PM7* (Eu-O correction), gas phase, Cartesian optimisation (XYZ EF, no O-H constraint); the campaign
      was stopped after 70 Sparkle/PM7* jobs; the keyword line actually run is recorded for every job
      -> "diagnostic runs without constraints showed exactly this" (proton transfer from water)
The raw output folders of these campaigns (~1.7 GB, mostly .aux trajectories) were not kept. This script was run
once on them (2026-10-13) to write:
  diagnostic_runs_summary.csv           one row per job that produced output
  diagnostic_runs_final_geometries.xyz  last geometry of each of those jobs (multi-structure XYZ)
The original input packages (eu_pectin_campaign_v1/v2/v5.tar.gz) are kept next to these files.
Descriptors are those of analyze.py (cutoff 3.0 A); d_min_any = shortest Eu-O distance (any oxygen);
collapsed = an oxygen within 2.0 A of Eu; n_hydroxide / n_oxonium = water oxygens that lost / gained a proton.

Usage (historical): python summarize_diagnostic_runs.py <campaign_v1_dir> <campaign2_dir> <campaign5_dir> <out_dir>
"""
import json
import os
import re
import sys

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import analyze                      # noqa: E402
from mopio import read_aux          # noqa: E402

CAMPAIGNS = [
    ("v1", "Sparkle/PM7 as distributed; COSMO EPS=78.4 (one control without COSMO); Cartesian (XYZ), no O-H constraint", "PM7"),
    ("v2", "Sparkle/PM7 as distributed; gas phase; O-H bonds frozen", "PM7"),
    ("v5", "Sparkle/PM7* (EXTERNAL Eu-O correction); gas phase; Cartesian optimisation (XYZ EF), no O-H constraint", None),
]
COLS = ["campaign", "protocol", "job", "model", "keywords_used", "status", "ncycles", "gnorm", "heat", "CN", "n_pectin", "n_water",
        "d_min_any", "collapsed", "n_hydroxide", "n_oxonium", "n_water_defects", "n_alkoxide", "n_acidOH", "n_carboxylate",
        "n_hydroxyl", "n_ring_or_glycosidic", "motif"]


def keywords_used(outfile):
    """Keyword line echoed by MOPAC in the .out file (the input actually run)."""
    for line in open(outfile, errors="ignore"):
        s = line.strip()
        if re.match(r"^(PM7|RM1|PM6|AM1)\s.*CHARGE=", s):
            return s
    return ""


def main(dirs, out):
    os.makedirs(out, exist_ok=True)
    rows = []
    xyz = open(os.path.join(out, "diagnostic_runs_final_geometries.xyz"), "w")
    for (tag, protocol, model), camp in zip(CAMPAIGNS, dirs):
        man = json.load(open(os.path.join(camp, "manifest.json")))
        for rec in man:
            rec = dict(rec)
            rec.setdefault("model", model or rec["job"].split("_")[0])
            if isinstance(rec.get("residues"), dict):          # v1/v2 manifests store residues as a dict
                rec["residues"] = [rec["residues"][k] for k in sorted(rec["residues"], key=int)]
            jd = os.path.join(camp, rec["job"])
            name, stage = analyze.last_stage(jd, rec["job"])
            if not os.path.exists(os.path.join(jd, name + ".aux")):
                continue
            r = analyze.analyze_job(camp, rec)
            if r.get("status") in ("missing", "noxyz"):
                continue
            r.update(campaign=tag, protocol=protocol, keywords_used=keywords_used(os.path.join(jd, name + ".out")))
            rows.append({c: r.get(c) for c in COLS})
            a = read_aux(os.path.join(jd, name + ".aux"))
            xyz.write(f"{len(a['el'])}\n{tag} {rec['job']} stage={stage} status={r.get('status')} d_min_EuO={r.get('d_min_any', float('nan')):.3f}\n")
            for e, p in zip(a["el"], a["xyz"]):
                xyz.write(f"{e:2s} {p[0]:12.6f} {p[1]:12.6f} {p[2]:12.6f}\n")
    xyz.close()
    df = pd.DataFrame(rows, columns=COLS)
    df.to_csv(os.path.join(out, "diagnostic_runs_summary.csv"), index=False)
    for tag, d in df.groupby("campaign", sort=False):
        c = d[d.CN.notna()]
        print(f"{tag}: {len(d)} jobs with output; status {d.status.value_counts().to_dict()}; Eu-O < 2.0 A in {int(c.collapsed.sum())} "
              f"(shortest {c.d_min_any.min():.2f} A); hydroxide formed in {int((c.n_hydroxide > 0).sum())} jobs")


if __name__ == "__main__":
    main(sys.argv[1:4], sys.argv[4])
