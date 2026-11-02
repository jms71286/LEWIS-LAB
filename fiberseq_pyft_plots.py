#!/usr/bin/env python3
"""
Fiber-seq visualization using pyft — Neurospora crassa
=======================================================
Covers all major plot types available through pyft plus
custom plots for nucleosome spacing heterogeneity analysis.

Install dependencies:
    pip install pyft altair vegafusion pandas numpy matplotlib seaborn

Usage:
    python fiberseq_pyft_plots.py \
        --bam  /path/to/sample.nucs.bam \
        --tss  /path/to/neurospora.bed \
        --out  /path/to/output_dir \
        --sample WT_Eddie

    # Compare two samples (WT vs mutant):
    python fiberseq_pyft_plots.py \
        --bam  /path/to/WT.nucs.bam \
        --bam2 /path/to/cac1.nucs.bam \
        --tss  /path/to/neurospora.bed \
        --out  /path/to/output_dir \
        --sample WT_Eddie \
        --sample2 cac-1
"""

import argparse
import os
import sys
import warnings
warnings.filterwarnings("ignore")

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.colors import LinearSegmentedColormap
import seaborn as sns

try:
    import pyft
    import altair as alt
    alt.data_transformers.enable("vegafusion")
except ImportError:
    sys.exit("ERROR: pyft not installed. Run: pip install pyft altair vegafusion")


# ── Color scheme ──────────────────────────────────────────────────────────────
COLORS = {
    "WT_Eddie":       "#2166AC",
    "WT_Rochelle":    "#4DAF4A",
    "cac-1":          "#D6604D",
    "cac-2":          "#F4A582",
    "rtt109":         "#762A83",
    "rtt109FLAG":     "#9970AB",
    "default_wt":     "#2166AC",
    "default_mut":    "#D6604D",
}

def get_color(sample):
    for key in COLORS:
        if key in sample:
            return COLORS[key]
    return "#888888"


# =============================================================================
# PLOT 1 — centered_chart (pyft native)
# =============================================================================
# The pyft centered_chart shows nucleosome (nuc) and MSP (accessible span)
# positions on each fiber, centered on a genomic locus of interest.
# Each row = one fiber. Blocks = nucleosomes. Gaps = accessible DNA (MSPs).
# This is the primary single-fiber visualization.
# =============================================================================

def plot_centered_chart(fiberbam, region, strand, sample_label, outdir,
                        max_flank=1000, max_fibers=200):
    """
    pyft centered_chart around a single locus.
    region = (chrom, start, end)
    """
    chrom, start, end = region
    print(f"  Fetching fibers at {chrom}:{start}-{end}...")

    df = pyft.utils.region_to_centered_df(
        fiberbam, region, strand=strand, max_flank=max_flank
    )

    if df is None or len(df) == 0:
        print(f"  ⚠️  No fibers found at {chrom}:{start}-{end}")
        return

    n_fibers = df["fiber_name"].nunique()
    print(f"  {n_fibers} fibers found")

    chart = pyft.plot.centered_chart(df, width=700, height=min(n_fibers * 3, 400))

    locus_str = f"{chrom}_{start}_{end}"
    out_html = os.path.join(outdir, f"{sample_label}.centered_{locus_str}.html")
    chart.save(out_html)
    print(f"  ✅  Saved: {out_html}  (open in browser)")
    return df


# =============================================================================
# PLOT 2 — TSS-centered aggregate: nucleosome + MSP occupancy profiles
# =============================================================================
# For each position relative to TSS, calculates:
#   - fraction of fibers with a nucleosome overlapping that position
#   - fraction of fibers with an MSP (accessible) overlapping that position
# This is the single-fiber equivalent of a deepTools profile plot.
# =============================================================================

def plot_tss_occupancy_profile(fiberbam, tss_df, sample_label, outdir,
                                flank=2000, bin_size=25, max_loci=500):
    """
    Aggregate nucleosome and MSP occupancy around all TSSs.
    """
    print(f"  Computing TSS occupancy profile (flank={flank}bp, {len(tss_df)} loci)...")

    bins = np.arange(-flank, flank + bin_size, bin_size)
    bin_mids = (bins[:-1] + bins[1:]) / 2
    nuc_counts = np.zeros(len(bin_mids))
    msp_counts = np.zeros(len(bin_mids))
    total_fibers = 0

    loci_used = 0
    for _, locus in tss_df.iterrows():
        if loci_used >= max_loci:
            break
        chrom = locus["chrom"]
        tss   = int(locus["tss"])
        strand = locus.get("strand", "+")

        region = (chrom, max(0, tss - flank), tss + flank)
        try:
            df = pyft.utils.region_to_centered_df(
                fiberbam, region, strand=strand, max_flank=flank
            )
        except Exception:
            continue

        if df is None or len(df) == 0:
            continue

        for fiber_name, fgrp in df.groupby("fiber_name"):
            total_fibers += 1
            nucs = fgrp[fgrp["type"] == "nuc"]
            msps = fgrp[fgrp["type"] == "msp"]

            for _, nuc in nucs.iterrows():
                mask = (bin_mids >= nuc["start"]) & (bin_mids <= nuc["end"])
                nuc_counts[mask] += 1

            for _, msp in msps.iterrows():
                mask = (bin_mids >= msp["start"]) & (bin_mids <= msp["end"])
                msp_counts[mask] += 1

        loci_used += 1

    if total_fibers == 0:
        print("  ⚠️  No fibers found for TSS profile")
        return

    nuc_occ = nuc_counts / total_fibers
    msp_occ = msp_counts / total_fibers

    fig, axes = plt.subplots(2, 1, figsize=(10, 7), sharex=True)
    fig.suptitle(f"{sample_label} — nucleosome & accessibility occupancy around TSS\n"
                 f"({loci_used} loci, {total_fibers} fiber-locus observations)",
                 fontsize=11, fontweight="bold")

    color = get_color(sample_label)

    axes[0].fill_between(bin_mids, nuc_occ, alpha=0.4, color=color)
    axes[0].plot(bin_mids, nuc_occ, color=color, lw=1.5)
    axes[0].axvline(0, color="black", lw=1, linestyle="--", alpha=0.6)
    axes[0].set_ylabel("Fraction of fibers\nwith nucleosome")
    axes[0].spines[["top", "right"]].set_visible(False)

    axes[1].fill_between(bin_mids, msp_occ, alpha=0.4, color="#2CA02C")
    axes[1].plot(bin_mids, msp_occ, color="#2CA02C", lw=1.5)
    axes[1].axvline(0, color="black", lw=1, linestyle="--", alpha=0.6)
    axes[1].set_ylabel("Fraction of fibers\nwith MSP (accessible)")
    axes[1].set_xlabel("Distance from TSS (bp)")
    axes[1].spines[["top", "right"]].set_visible(False)

    plt.tight_layout()
    out = os.path.join(outdir, f"{sample_label}.TSS_occupancy_profile.png")
    plt.savefig(out, dpi=180, bbox_inches="tight")
    plt.close()
    print(f"  ✅  Saved: {out}")


# =============================================================================
# PLOT 3 — Single-fiber heatmap sorted by +1 nucleosome position
# =============================================================================
# Each row = one fiber. Each column = position bin relative to TSS.
# Color = nucleosome present (dark) vs accessible (light).
# Fibers sorted by the position of the +1 nucleosome.
# Heterogeneous samples will show a "fuzzy" band; well-phased will be sharp.
# =============================================================================

def plot_fiber_heatmap_tss(fiberbam, tss_df, sample_label, outdir,
                            flank=1500, bin_size=25, max_fibers=300,
                            max_loci=100):
    """
    Single-fiber nucleosome heatmap centered on TSS, sorted by +1 nuc position.
    """
    print(f"  Building single-fiber heatmap ({max_fibers} fibers max)...")

    bins     = np.arange(-flank, flank + bin_size, bin_size)
    bin_mids = (bins[:-1] + bins[1:]) / 2
    rows     = []
    plus1_positions = []

    for _, locus in tss_df.iterrows():
        if len(rows) >= max_fibers:
            break
        chrom  = locus["chrom"]
        tss    = int(locus["tss"])
        strand = locus.get("strand", "+")
        region = (chrom, max(0, tss - flank), tss + flank)

        try:
            df = pyft.utils.region_to_centered_df(
                fiberbam, region, strand=strand, max_flank=flank
            )
        except Exception:
            continue
        if df is None or len(df) == 0:
            continue

        for fiber_name, fgrp in df.groupby("fiber_name"):
            if len(rows) >= max_fibers:
                break
            row = np.zeros(len(bin_mids))
            nucs = fgrp[fgrp["type"] == "nuc"]
            for _, nuc in nucs.iterrows():
                mask = (bin_mids >= nuc["start"]) & (bin_mids <= nuc["end"])
                row[mask] = 1.0

            rows.append(row)

            # +1 nuc = closest downstream nucleosome
            downstream_nucs = nucs[nucs["start"] > 0]
            if len(downstream_nucs) > 0:
                plus1_positions.append(downstream_nucs["start"].min())
            else:
                plus1_positions.append(np.nan)

    if len(rows) == 0:
        print("  ⚠️  No fibers for heatmap")
        return

    matrix = np.array(rows)
    # sort by +1 nuc position
    order  = np.argsort([p if not np.isnan(p) else 9999 for p in plus1_positions])
    matrix = matrix[order]

    color  = get_color(sample_label)
    cmap   = LinearSegmentedColormap.from_list("nuc", ["#f7f7f7", color])

    fig, ax = plt.subplots(figsize=(10, max(4, len(rows) * 0.025)))
    ax.imshow(matrix, aspect="auto", cmap=cmap,
              extent=[bin_mids[0], bin_mids[-1], 0, len(rows)],
              interpolation="nearest", vmin=0, vmax=1)
    ax.axvline(0, color="black", lw=1.2, linestyle="--", label="TSS")
    ax.set_xlabel("Distance from TSS (bp)")
    ax.set_ylabel(f"Individual fibers (n={len(rows)}, sorted by +1 nuc)")
    ax.set_title(f"{sample_label} — single-fiber nucleosome map\nsorted by +1 nucleosome position",
                 fontsize=11)
    ax.legend(loc="upper right", fontsize=9)

    plt.tight_layout()
    out = os.path.join(outdir, f"{sample_label}.fiber_heatmap_TSS.png")
    plt.savefig(out, dpi=180, bbox_inches="tight")
    plt.close()
    print(f"  ✅  Saved: {out}")


# =============================================================================
# PLOT 4 — MSP length distribution
# =============================================================================
# MSP = methylated span = accessible element (linker DNA or NDR).
# Distribution of MSP lengths tells you:
#   - peak ~38-80 bp = short linker/TF-bound elements
#   - peak ~150-200 bp = nucleosome-sized gaps (NDRs)
#   - compare WT vs mutant to detect global changes in accessibility
# =============================================================================

def plot_msp_length_distribution(fiberbam, sample_label, outdir,
                                  max_fibers=5000, color=None):
    """MSP length histogram from all fibers genome-wide."""
    print(f"  Collecting MSP lengths (up to {max_fibers} fibers)...")

    msp_lengths = []
    nuc_lengths = []
    count = 0

    for fiber in fiberbam:
        if count >= max_fibers:
            break
        msp_lengths.extend(fiber.msp.reference_lengths.tolist())
        nuc_lengths.extend(fiber.nuc.reference_lengths.tolist())
        count += 1

    if not msp_lengths:
        print("  ⚠️  No MSP data found")
        return

    color = color or get_color(sample_label)

    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    fig.suptitle(f"{sample_label} — MSP and nucleosome length distributions\n"
                 f"({count} fibers)", fontsize=11, fontweight="bold")

    # MSP lengths
    axes[0].hist(msp_lengths, bins=np.arange(0, 600, 5), color=color,
                 alpha=0.7, density=True)
    axes[0].axvline(38,  color="gray", lw=1, linestyle=":", label="38 bp (min MSP)")
    axes[0].axvline(150, color="gray", lw=1, linestyle="--", label="150 bp (NDR)")
    axes[0].set_xlabel("MSP length (bp)")
    axes[0].set_ylabel("Density")
    axes[0].set_title("Accessible span (MSP) lengths")
    axes[0].set_xlim(0, 600)
    axes[0].legend(fontsize=8)
    axes[0].spines[["top", "right"]].set_visible(False)

    # Nucleosome lengths
    axes[1].hist(nuc_lengths, bins=np.arange(50, 350, 5), color="#D6604D",
                 alpha=0.7, density=True)
    axes[1].axvline(147, color="gray", lw=1, linestyle="--", label="147 bp (canonical)")
    axes[1].set_xlabel("Nucleosome footprint length (bp)")
    axes[1].set_ylabel("Density")
    axes[1].set_title("Nucleosome footprint lengths")
    axes[1].set_xlim(50, 350)
    axes[1].legend(fontsize=8)
    axes[1].spines[["top", "right"]].set_visible(False)

    plt.tight_layout()
    out = os.path.join(outdir, f"{sample_label}.MSP_nuc_length_distributions.png")
    plt.savefig(out, dpi=180, bbox_inches="tight")
    plt.close()
    print(f"  ✅  Saved: {out}")

    return msp_lengths, nuc_lengths


# =============================================================================
# PLOT 5 — WT vs Mutant MSP length overlay
# =============================================================================

def plot_msp_comparison(msp_wt, msp_mut, nuc_wt, nuc_mut,
                         label_wt, label_mut, outdir):
    """Overlay MSP and nucleosome length distributions for two samples."""
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    fig.suptitle(f"MSP and nucleosome lengths: {label_wt} vs {label_mut}",
                 fontsize=11, fontweight="bold")

    color_wt  = get_color(label_wt)
    color_mut = get_color(label_mut)

    for ax, wt_data, mut_data, xmax, xlabel, title in [
        (axes[0], msp_wt, msp_mut, 600, "MSP length (bp)", "Accessible span lengths"),
        (axes[1], nuc_wt, nuc_mut, 350, "Nucleosome length (bp)", "Nucleosome footprint lengths"),
    ]:
        bins = np.arange(0, xmax, 5)
        ax.hist(wt_data,  bins=bins, color=color_wt,  alpha=0.5,
                density=True, label=label_wt)
        ax.hist(mut_data, bins=bins, color=color_mut, alpha=0.5,
                density=True, label=label_mut)
        ax.set_xlabel(xlabel)
        ax.set_ylabel("Density")
        ax.set_title(title)
        ax.set_xlim(0, xmax)
        ax.legend(fontsize=9)
        ax.spines[["top", "right"]].set_visible(False)

    plt.tight_layout()
    out = os.path.join(outdir, f"{label_wt}_vs_{label_mut}.MSP_comparison.png")
    plt.savefig(out, dpi=180, bbox_inches="tight")
    plt.close()
    print(f"  ✅  Saved: {out}")


# =============================================================================
# PLOT 6 — Per-fiber m6A fraction distribution
# =============================================================================
# Checks Hia5 labeling quality. Each fiber gets a score = m6A / total_A.
# Should peak around 5-6%. Bimodal distribution or high zero-peak = problem.
# =============================================================================

def plot_m6a_fraction(fiberbam, sample_label, outdir,
                       max_fibers=5000, color=None):
    """Per-fiber m6A fraction histogram — QC metric."""
    print(f"  Computing per-fiber m6A fractions...")

    fracs = []
    counts = []
    count = 0

    for fiber in fiberbam:
        if count >= max_fibers:
            break
        seq_len = fiber.get_seq_length()
        if seq_len == 0:
            continue
        n_m6a = len(fiber.m6a.starts)
        fracs.append(n_m6a / seq_len * 100)
        counts.append(n_m6a)
        count += 1

    if not fracs:
        print("  ⚠️  No m6A data found")
        return

    color = color or get_color(sample_label)
    zero_frac = sum(1 for f in fracs if f == 0) / len(fracs) * 100

    fig, axes = plt.subplots(1, 2, figsize=(12, 4))
    fig.suptitle(f"{sample_label} — m6A labeling QC\n"
                 f"({count} fibers, {zero_frac:.1f}% with zero m6A)",
                 fontsize=11, fontweight="bold")

    axes[0].hist(fracs, bins=50, color=color, alpha=0.7, density=True)
    axes[0].axvline(5, color="gray", lw=1, linestyle="--", label="Expected ~5%")
    axes[0].axvline(6, color="gray", lw=1, linestyle="--")
    axes[0].set_xlabel("m6A / total bases (%)")
    axes[0].set_ylabel("Density")
    axes[0].set_title("Per-fiber m6A fraction\n(should peak at ~5-6%)")
    axes[0].legend(fontsize=8)
    axes[0].spines[["top", "right"]].set_visible(False)

    axes[1].hist(counts, bins=50, color=color, alpha=0.7)
    axes[1].set_xlabel("m6A count per fiber")
    axes[1].set_ylabel("Number of fibers")
    axes[1].set_title("Per-fiber m6A count")
    axes[1].spines[["top", "right"]].set_visible(False)

    # QC annotation
    median_frac = np.median(fracs)
    axes[0].text(0.95, 0.95, f"Median: {median_frac:.2f}%\n0 m6A: {zero_frac:.1f}%",
                 transform=axes[0].transAxes, ha="right", va="top",
                 fontsize=9, bbox=dict(boxstyle="round", facecolor="white", alpha=0.7))

    plt.tight_layout()
    out = os.path.join(outdir, f"{sample_label}.m6A_QC.png")
    plt.savefig(out, dpi=180, bbox_inches="tight")
    plt.close()
    print(f"  ✅  Saved: {out}")
    print(f"      Median m6A fraction: {median_frac:.2f}%  |  Zero-m6A fibers: {zero_frac:.1f}%")


# =============================================================================
# PLOT 7 — Nucleosomes-per-fiber distribution
# =============================================================================

def plot_nucs_per_fiber(fiberbam, sample_label, outdir,
                         max_fibers=5000, color=None):
    """Distribution of nucleosome count per fiber — reflects read length + occupancy."""
    print(f"  Computing nucleosomes per fiber...")

    n_nucs = []
    read_lengths = []
    count = 0

    for fiber in fiberbam:
        if count >= max_fibers:
            break
        n_nucs.append(len(fiber.nuc.starts))
        read_lengths.append(fiber.get_seq_length())
        count += 1

    if not n_nucs:
        return

    color = color or get_color(sample_label)

    fig, axes = plt.subplots(1, 2, figsize=(11, 4))
    fig.suptitle(f"{sample_label} — nucleosomes per fiber ({count} fibers)",
                 fontsize=11, fontweight="bold")

    axes[0].hist(n_nucs, bins=40, color=color, alpha=0.7)
    axes[0].axvline(np.median(n_nucs), color="black", lw=1.5, linestyle="--",
                    label=f"Median: {np.median(n_nucs):.0f}")
    axes[0].set_xlabel("Nucleosomes per fiber")
    axes[0].set_ylabel("Number of fibers")
    axes[0].set_title("Nucleosomes per fiber")
    axes[0].legend(fontsize=9)
    axes[0].spines[["top", "right"]].set_visible(False)

    # nucleosome density vs read length scatter
    axes[1].scatter(read_lengths, n_nucs, alpha=0.1, s=3, color=color)
    axes[1].set_xlabel("Read length (bp)")
    axes[1].set_ylabel("Nucleosome count")
    axes[1].set_title("Nucleosome count vs read length")
    axes[1].spines[["top", "right"]].set_visible(False)

    plt.tight_layout()
    out = os.path.join(outdir, f"{sample_label}.nucs_per_fiber.png")
    plt.savefig(out, dpi=180, bbox_inches="tight")
    plt.close()
    print(f"  ✅  Saved: {out}")


# =============================================================================
# Utility: load TSS BED
# =============================================================================

def load_tss_bed(path, max_loci=None):
    cols = ["chrom", "start", "end", "name", "score", "strand"]
    n = sum(1 for _ in open(path))
    ncols = len(open(path).readline().split("\t"))
    df = pd.read_csv(path, sep="\t", header=None,
                     names=cols[:ncols], dtype={"chrom": str})
    df["tss"] = df.apply(
        lambda r: r["end"] if r.get("strand") == "-" else r["start"], axis=1
    )
    if max_loci:
        df = df.head(max_loci)
    print(f"  Loaded {len(df)} loci from {path}")
    return df


# =============================================================================
# Main
# =============================================================================

def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--bam",     required=True, help="*.nucs.bam for sample 1 (or WT)")
    parser.add_argument("--bam2",    default=None,  help="*.nucs.bam for sample 2 (mutant) — optional")
    parser.add_argument("--tss",     required=True, help="TSS BED file (neurospora.bed)")
    parser.add_argument("--out",     required=True, help="Output directory")
    parser.add_argument("--sample",  default="sample1", help="Name for sample 1")
    parser.add_argument("--sample2", default="sample2", help="Name for sample 2")
    parser.add_argument("--locus",   default=None,
                        help="Single locus for centered chart, e.g. CM002237.1:1866000-1874000")
    parser.add_argument("--max-fibers", default=3000, type=int,
                        help="Max fibers for genome-wide plots (default 3000)")
    parser.add_argument("--max-loci", default=200, type=int,
                        help="Max TSS loci for aggregate plots (default 200)")
    args = parser.parse_args()

    os.makedirs(args.out, exist_ok=True)

    print(f"\n=== Loading {args.sample} ===")
    fiberbam1 = pyft.Fiberbam(args.bam, threads=8)
    tss_df    = load_tss_bed(args.tss, max_loci=args.max_loci)

    fiberbam2 = None
    if args.bam2:
        print(f"\n=== Loading {args.sample2} ===")
        fiberbam2 = pyft.Fiberbam(args.bam2, threads=8)

    # ── Plot 1: centered_chart at a specific locus ───────────────
    if args.locus:
        print(f"\n=== Plot 1: Centered chart at {args.locus} ===")
        chrom, coords = args.locus.split(":")
        start, end    = [int(x) for x in coords.replace(",", "").split("-")]
        region = (chrom, start, end)
        plot_centered_chart(fiberbam1, region, "+", args.sample, args.out)
        if fiberbam2:
            plot_centered_chart(fiberbam2, region, "+", args.sample2, args.out)

    # ── Plot 2: TSS occupancy profiles ───────────────────────────
    print(f"\n=== Plot 2: TSS occupancy profile ===")
    plot_tss_occupancy_profile(fiberbam1, tss_df, args.sample, args.out,
                                max_loci=args.max_loci)
    if fiberbam2:
        plot_tss_occupancy_profile(fiberbam2, tss_df, args.sample2, args.out,
                                    max_loci=args.max_loci)

    # ── Plot 3: Single-fiber heatmap ─────────────────────────────
    print(f"\n=== Plot 3: Single-fiber heatmap ===")
    plot_fiber_heatmap_tss(fiberbam1, tss_df, args.sample, args.out,
                            max_loci=args.max_loci // 5)
    if fiberbam2:
        plot_fiber_heatmap_tss(fiberbam2, tss_df, args.sample2, args.out,
                                max_loci=args.max_loci // 5)

    # ── Plot 4+5: MSP / nucleosome length distributions ──────────
    print(f"\n=== Plot 4: MSP + nucleosome length distributions ===")
    result1 = plot_msp_length_distribution(fiberbam1, args.sample, args.out,
                                            max_fibers=args.max_fibers)
    result2 = None
    if fiberbam2:
        result2 = plot_msp_length_distribution(fiberbam2, args.sample2, args.out,
                                                max_fibers=args.max_fibers)

    if result1 and result2:
        print(f"\n=== Plot 5: WT vs Mutant MSP comparison ===")
        plot_msp_comparison(result1[0], result2[0],
                             result1[1], result2[1],
                             args.sample, args.sample2, args.out)

    # ── Plot 6: m6A QC ───────────────────────────────────────────
    print(f"\n=== Plot 6: m6A labeling QC ===")
    # Re-open to reset iterator
    fiberbam1 = pyft.Fiberbam(args.bam, threads=8)
    plot_m6a_fraction(fiberbam1, args.sample, args.out,
                       max_fibers=args.max_fibers)
    if args.bam2:
        fiberbam2 = pyft.Fiberbam(args.bam2, threads=8)
        plot_m6a_fraction(fiberbam2, args.sample2, args.out,
                           max_fibers=args.max_fibers)

    # ── Plot 7: Nucleosomes per fiber ────────────────────────────
    print(f"\n=== Plot 7: Nucleosomes per fiber ===")
    fiberbam1 = pyft.Fiberbam(args.bam, threads=8)
    plot_nucs_per_fiber(fiberbam1, args.sample, args.out,
                         max_fibers=args.max_fibers)
    if args.bam2:
        fiberbam2 = pyft.Fiberbam(args.bam2, threads=8)
        plot_nucs_per_fiber(fiberbam2, args.sample2, args.out,
                             max_fibers=args.max_fibers)

    print(f"\n✅  All done. Outputs in: {args.out}/")
    print("""
Summary of outputs:
  *.centered_*.html              — interactive single-fiber view (open in browser)
  *.TSS_occupancy_profile.png   — nucleosome + MSP occupancy around TSS
  *.fiber_heatmap_TSS.png       — single-fiber heatmap sorted by +1 nuc position
  *.MSP_nuc_length_distributions.png — MSP and nucleosome length histograms
  *_vs_*.MSP_comparison.png     — WT vs mutant MSP overlay
  *.m6A_QC.png                  — Hia5 labeling QC (m6A fraction per fiber)
  *.nucs_per_fiber.png          — nucleosome count per fiber
""")


if __name__ == "__main__":
    main()
