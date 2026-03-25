#!/usr/bin/env python3
"""
Automated performance measurement for LINPACK benchmark.
"""

import subprocess
import sys
import re
import time
import csv
import statistics
from pathlib import Path
Path("results").mkdir(exist_ok=True)
SCRIPT_DIR = Path(__file__).parent.resolve()

# Configuration
N_RUNS = 20
# Binary files
ORIG_BIN =  SCRIPT_DIR / "linpack"
LLVM_BIN =  SCRIPT_DIR / "llvm_linpack"
TIGRESS_BIN =  SCRIPT_DIR / "tigress_linpack"
VMPROTECT_BIN =  SCRIPT_DIR / "vmprotect_linpack"
THEMIDA_BIN =  SCRIPT_DIR / "themida_linpack"

def run_benchmark(exe, N):
    """Run the benchmark and return (time_seconds, mflops)."""
    result = subprocess.run([exe, str(N)], capture_output=True, text=True)
    output = result.stdout + result.stderr
    # Parse time and MFLOPS from output
    time_match = re.search(r"Time in seconds\s*=\s*([0-9.eE+-]+)", output)
    mflops_match = re.search(r"MegaFLOPS\s*=\s*([0-9.eE+-]+)", output)
    if time_match and mflops_match:
        t = float(time_match.group(1))
        mf = float(mflops_match.group(1))
        return t, mf
    else:
        print("Failed to parse output:")
        print(output)
        return None, None

def measure(exe, runs, N):
    """Run benchmark multiple times, return list of (t, mf)."""
    results = []
    for i in range(runs):
        print(f"  Run {i+1}/{runs}...")
        t, mf = run_benchmark(exe, N)
        if t is not None:
            results.append((t, mf))
        else:
            print("Skipping this run due to error.")
        # Small pause
        time.sleep(0.5)
    return results

def save_csv(data, filename, labels):
    """Save results to CSV file."""
    with open(filename, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["Configuration", "Run", "Time(s)", "MFLOPS"])
        for label, runs in zip(labels, data):
            for i, (t, mf) in enumerate(runs):
                writer.writerow([label, i+1, t, mf])

def print_stats(data, labels):
    """Print statistics for each configuration."""
    for label, runs in zip(labels, data):
        times = [t for t,_ in runs]
        mflops = [mf for _,mf in runs]
        if times:
            print(f"\n{label}:")
            print(f"  Time: mean={statistics.mean(times):.3f}s, std={statistics.stdev(times):.3f}s")
            print(f"  MFLOPS: mean={statistics.mean(mflops):.1f}, std={statistics.stdev(mflops):.1f}")

def main():
    user_input = input("Enter matrix size: ")
    N = int(user_input)
    LDA = N + 1
    # run measurements
    print("\nMeasuring original version...")
    orig_results = measure(ORIG_BIN, N_RUNS, N)
    print("\nMeasuring Tigress version...")
    tigress_results = measure(TIGRESS_BIN, N_RUNS, N)
    print("\nMeasuring VMProtect version...")
    vmprotect_results = measure(VMPROTECT_BIN, N_RUNS, N)
    # themida binary file crashes with Segfault
    #print("\nMeasuring Themida version...")
    #themida_results = measure(THEMIDA_BIN, N_RUNS, N)
    print("\nMeasuring my version...")
    llvm_results = measure(LLVM_BIN, N_RUNS, N)

    # Save results
    all_results = [orig_results, llvm_results, tigress_results, vmprotect_results]
    labels = ["Original", "Softcom LLVM obfuscator", "Tigress", "VMProtect"]
    save_csv(all_results, "results/linpack_bench_results_" + str(N) + ".csv", labels)

    # Print statistics
    print_stats(all_results, labels)

if __name__ == "__main__":
    main()
