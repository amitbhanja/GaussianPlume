#!/usr/bin/env python3
"""Summarise Output/forward_model_data.csv: measured vs. modelled per stability class.

Usage (from the repository root):
    python3 scripts/summarize_output.py [path/to/forward_model_data.csv]

Only the Python standard library is used.
"""
import csv
import math
import sys

CLASSES = "ABCDEF"


def main(path="Output/forward_model_data.csv"):
    rows = list(csv.DictReader(open(path, newline="")))
    print(f"{path}: {len(rows)} rows ({len(rows) // 6} sensor readings x 6 classes)\n")
    print("Class  non-zero   mean meas.   mean model   max model    RMSE      median model/meas")
    print("                  (g/m3)       (g/m3)       (g/m3)       (g/m3)    (non-zero rows)")
    for k in range(6):
        pairs = [
            (float(r["Sensor Measurement"]), float(r["Forward-Model-Concentration"]))
            for r in rows
            if r["Stability"] == str(k)
        ]
        meas = [m for m, _ in pairs]
        model = [p for _, p in pairs]
        ratios = sorted(p / m for m, p in pairs if p > 0 and m > 0)
        rmse = math.sqrt(sum((m - p) ** 2 for m, p in pairs) / len(pairs))
        median = f"{ratios[len(ratios) // 2]:.3g}" if ratios else "-"
        print(
            f"  {CLASSES[k]}    {sum(p > 0 for p in model):4d}/{len(pairs)}"
            f"   {sum(meas) / len(meas):.3e}    {sum(model) / len(model):.3e}"
            f"    {max(model):.3e}    {rmse:.3e}   {median}"
        )


if __name__ == "__main__":
    main(*sys.argv[1:2])
