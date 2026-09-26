#!/usr/bin/env python3
"""Rate of the monotonic clock against real time, for bench/harness.py --samaya-clock-factor.

samaya times its limit on the monotonic clock (std::chrono::steady_clock), while HiGHS, SCIP, CBC
and GLPK stop on the real-time clock. On WSL2 on some laptops the monotonic clock runs a few
percent fast, and the real-time clock is corrected every half minute, so at the same nominal limit
samaya would get less real time. This measures monotonic seconds per real second over a window
(Windows' own clock through powershell.exe when it is there, else the real-time clock) and prints
the factor, rounded to 3 decimals.

    python3 bench/clock_factor.py [SECONDS]
"""
import shutil
import subprocess
import sys
import time


def windows_seconds():
    out = subprocess.run(["powershell.exe", "-NoProfile", "-Command", "[DateTime]::UtcNow.Ticks"],
                         capture_output=True, text=True, timeout=60).stdout.strip()
    return int(out) / 1e7


def main():
    window = float(sys.argv[1]) if len(sys.argv) > 1 else 120.0
    real = windows_seconds if shutil.which("powershell.exe") else time.time
    r0, m0 = real(), time.monotonic()
    time.sleep(window)
    r1, m1 = real(), time.monotonic()
    factor = (m1 - m0) / (r1 - r0)
    print(f"{factor:.3f}")
    print(f"monotonic {m1 - m0:.2f} s over {r1 - r0:.2f} s of real time "
          f"({'Windows clock' if real is windows_seconds else 'real-time clock'})", file=sys.stderr)


if __name__ == "__main__":
    main()
