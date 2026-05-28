#!/usr/bin/env python3
"""Parse strong-scaling SLURM logs into a tidy CSV.

Reads apps/logs/scaling/scaling_*.log produced by run_strong_scaling.sh and
emits apps/output/data/scaling_data.csv with one row per (ranks, case,
trial) measurement.
"""
from __future__ import annotations

import re
import sys
from pathlib import Path

import pandas as pd

_HERE     = Path(__file__).parent              # apps/scripts/
APPS_DIR  = _HERE.parent                       # apps/
LOG_DIR   = APPS_DIR / "logs" / "scaling"
OUT_CSV   = APPS_DIR / "output" / "data" / "scaling_data.csv"

_RE_HEADER     = re.compile(r"── (\S+)\s+trial (\d+)/\d+\s+ranks=(\d+) ──")
_RE_MLUPS_TOT  = re.compile(r"\[\s*0\s*\]\[RESULT\s*\].*?Total MLUPS:\s*([\d.]+)")
_RE_MLUPS_PROC = re.compile(r"\[\s*0\s*\]\[RESULT\s*\].*?MLUPS per process:\s*([\d.]+)")
_RE_WALL       = re.compile(r"\[\s*0\s*\]\[RESULT\s*\].*?Total simulation time:\s*([\d.]+)\s+seconds")


def parse_log(path: Path) -> list[dict]:
    rows: list[dict] = []
    text = path.read_text(errors="replace")
    # Split on the trial header so each segment contains exactly one run's RESULT lines.
    segments = re.split(_RE_HEADER, text)
    # re.split with capture groups yields: [pre, case1, trial1, ranks1, body1, case2, trial2, ranks2, body2, ...]
    if len(segments) < 5:
        return rows
    for case, trial, ranks, body in zip(segments[1::4], segments[2::4], segments[3::4], segments[4::4]):
        m_tot   = _RE_MLUPS_TOT.search(body)
        m_proc  = _RE_MLUPS_PROC.search(body)
        m_wall  = _RE_WALL.search(body)
        if not (m_tot and m_wall):
            continue
        rows.append({
            "case":          case.strip(),
            "ranks":         int(ranks),
            "trial":         int(trial),
            "mlups_total":   float(m_tot.group(1)),
            "mlups_per_proc": float(m_proc.group(1)) if m_proc else float("nan"),
            "wall_time_s":   float(m_wall.group(1)),
        })
    return rows


def main() -> int:
    logs = sorted(LOG_DIR.glob("scaling_*.log"))
    if not logs:
        print(f"No scaling logs found in {LOG_DIR}", file=sys.stderr)
        return 1
    all_rows: list[dict] = []
    for log in logs:
        rs = parse_log(log)
        print(f"  {log.name}: {len(rs)} measurements")
        all_rows.extend(rs)
    if not all_rows:
        print("No measurements parsed (likely jobs still queued or failed).", file=sys.stderr)
        return 1
    df = pd.DataFrame(all_rows).sort_values(["case", "ranks", "trial"]).reset_index(drop=True)
    OUT_CSV.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(OUT_CSV, index=False, float_format="%.4f")
    print(f"\nWrote {OUT_CSV}  ({len(df)} rows)")
    print(df.to_string(index=False))
    return 0


if __name__ == "__main__":
    sys.exit(main())
