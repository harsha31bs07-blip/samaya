#!/usr/bin/env python3
"""HiGHS presolve sizes with one presolve rule switched off at a time.

Shows which HiGHS rule removes the columns samaya keeps on an instance (sizes only, no solve):

    python3 bench/highs_rule_sizes.py bench/instances/miplib_small app1-1 b1c1s1

Rules whose absence changes nothing are not printed. The bits are HiGHS's presolve_rule_off
bits, as in bench/attribute.py.
"""
import sys

import highspy

from attribute import HIGHS_PRESOLVE_RULE_BITS


def presolved_size(path: str, off: int) -> tuple[int, int]:
    h = highspy.Highs()
    h.setOptionValue("output_flag", False)
    h.setOptionValue("threads", 1)
    if off:
        h.setOptionValue("presolve_rule_off", off)
    h.readModel(path)
    h.presolve()
    lp = h.getPresolvedLp()
    return lp.num_row_, lp.num_col_


def main() -> int:
    if len(sys.argv) < 3:
        print(__doc__)
        return 2
    for name in sys.argv[2:]:
        path = f"{sys.argv[1]}/{name}.mps"
        rows, cols = presolved_size(path, 0)
        print(f"{name}: default rows {rows} cols {cols}")
        for rule, bit in HIGHS_PRESOLVE_RULE_BITS.items():
            r, c = presolved_size(path, 1 << bit)
            if (r, c) != (rows, cols):
                print(f"  {rule:24s} off: rows {r} cols {c} ({c - cols:+d} cols)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
