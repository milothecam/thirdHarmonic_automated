"""
Standalone analysis script for harmonic generation experiment.

Reads CSV files produced by harmonic_measurement.py and generates:
  - One figure per fundamental frequency with one curve per harmonic
    X-axis: P_TKIP_in [dBm], Y-axis: P_TKIP_out [dBm]
  - Single temperature time-series plot (fixed width)

Power corrections (optional, set APPLY_CORRECTIONS = True):

  All S-parameter and attenuation values are negative magnitudes.

  loss(f)       = 0.5 * [S2P_CRYO_IN_OUT(f) - ATT_in]

  P_TKIP_IN(f)  = P_GEN + S2P_PUMP_LINE(f) + S2P_PUMP_IN(f)
                        + 0.5*ATT_in - 0.5*S2P_CRYO_IN_OUT(f)

  P_TKIP_OUT(f,h) = P_out - loss(f*h)
                  = P_out - 0.5*S2P_CRYO_IN_OUT(f*h) + 0.5*ATT_in
"""

import pathlib
import csv
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import skrf as rf

############################
# USER CONFIGURATION
# Must match harmonic_measurement.py
############################

frequency_list = [2, 3, 4, 5]            # GHz — fundamental tone frequencies
power_list     = np.linspace(15, 25, 51).tolist()  # dBm — signal generator output powers
harmonic_list  = [1, 3]            # harmonic index (multiplier of fundamental)

SAVE_DIR      = pathlib.Path(r"C:\data\Camilo\Wafer UVA microstrip\microstrip E\harmonicsAutomated\microstripE_2345GHz_3rdH_10mK_fineStep")
FILE_PREFIX   = "harmAuto"

PLOT_SAVE_DIR = SAVE_DIR / "plotsCorr"
PLOT_SAVE_DIR.mkdir(parents=True, exist_ok=True)

############################
# POWER CORRECTION CONFIGURATION
#
# Set APPLY_CORRECTIONS = False to plot raw (uncorrected) data.
#
# All S-parameter values and attenuations are negative magnitudes [dB].
#
# Three s2p files:
#   S2P_PUMP_LINE  : A→B  outside cryostat, signal generator side
#   S2P_PUMP_IN    : B→D  inside cryostat, pump input path
#   S2P_CRYO_IN_OUT: C→D  VNA in to VNA out, full cryostat path
#
# ATT_p  : attenuator on pump line        (negative dB, e.g. -20)
# ATT_in : attenuator on VNA/input line   (negative dB, e.g. -20)
#
# Correction formulas:
#   loss(f)        = 0.5 * [S2P_CRYO_IN_OUT(f) - ATT_in]
#   P_TKIP_IN(f)   = P_GEN + S2P_PUMP_LINE(f) + S2P_PUMP_IN(f)
#                          + 0.5*ATT_in - 0.5*S2P_CRYO_IN_OUT(f)
#   P_TKIP_OUT(f,h)= P_out - loss(f*h)
#                  = P_out - 0.5*S2P_CRYO_IN_OUT(f*h) + 0.5*ATT_in
############################

APPLY_CORRECTIONS = True

S2P_PUMP_LINE    = pathlib.Path(r"C:\data\Camilo\Wafer UVA microstrip\microstrip E\harmonicsAutomated\pumpPath.s2p")
S2P_PUMP_IN      = pathlib.Path(r"C:\data\Camilo\Wafer UVA microstrip\microstrip E\harmonicsAutomated\S21_pumpOut_oneBiasTee_20001pts.s2p")
S2P_CRYO_IN_OUT  = pathlib.Path(r"C:\data\Camilo\Wafer UVA microstrip\microstrip E\harmonicsAutomated\S21_inOut_wBiasTees_20001pts.s2p")

ATT_p   = -20.0   # dB — attenuator on pump line   (negative magnitude)
ATT_in  = -50.0   # dB — attenuator on input line  (negative magnitude)

# Smoothing window for s2p curves (number of points, must be odd; set to 1 to disable)
SMOOTH_WINDOW = 501

############################
# LOAD S2P FILES + INTERPOLATORS
############################

def load_and_smooth(s2p_path, port_i=0, port_j=1, window=SMOOTH_WINDOW):
    """
    Load an s2p file, smooth the S-parameter magnitude using a moving average,
    and return a callable that interpolates magnitude [dB] at any frequency [Hz].
    port_i, port_j: 0-indexed skrf port indices for the desired S-parameter.
    window: number of points for the moving average (must be odd; 1 = no smoothing).
    """
    ntwk  = rf.Network(str(s2p_path))
    freqs = ntwk.f
    raw   = 20.0 * np.log10(np.abs(ntwk.s[:, port_i, port_j]) + 1e-300)

    if window > 1:
        kernel = np.ones(window) / window
        s_db   = np.convolve(raw, kernel, mode='same')
        # Restore edges where convolution boundary effects occur
        half         = window // 2
        s_db[:half]  = raw[:half]
        s_db[-half:] = raw[-half:]
    else:
        s_db = raw

    def interpolate(freq_hz):
        return float(np.interp(freq_hz, freqs, s_db))

    return interpolate, freqs, s_db


if APPLY_CORRECTIONS:
    print("Loading and smoothing s2p correction files...")
    get_pump_line,   f_pl,  s_pl  = load_and_smooth(S2P_PUMP_LINE)
    get_pump_in,     f_pi,  s_pi  = load_and_smooth(S2P_PUMP_IN)
    get_cryo_in_out, f_cio, s_cio = load_and_smooth(S2P_CRYO_IN_OUT)
    print(f"  S2P_PUMP_LINE   : {S2P_PUMP_LINE.name}")
    print(f"  S2P_PUMP_IN     : {S2P_PUMP_IN.name}")
    print(f"  S2P_CRYO_IN_OUT : {S2P_CRYO_IN_OUT.name}")
    print(f"  ATT_p  = {ATT_p} dB")
    print(f"  ATT_in = {ATT_in} dB")
    print(f"  Smoothing window: {SMOOTH_WINDOW} points")

    def loss(freq_hz):
        """One-way cryostat loss at freq_hz [Hz]. Returns negative dB value."""
        return 0.5 * (get_cryo_in_out(freq_hz) - ATT_in)

    def corrected_x(power_dbm, freq_ghz):
        """P_TKIP_IN: correct input power for all line losses up to TKIP input."""
        f = freq_ghz * 1e9
        return (power_dbm
                + get_pump_line(f)
                + get_pump_in(f)
                + 0.5 * ATT_in
                - 0.5 * get_cryo_in_out(f))

    def corrected_y(peak_db, freq_ghz, harmonic):
        """P_TKIP_OUT: correct measured output for cryostat loss at harmonic freq."""
        f = freq_ghz * harmonic * 1e9
        return peak_db - loss(f)

    x_label      = "P_TKIP_in [dBm]"
    y_label      = "P_TKIP_out [dBm]"
    fname_suffix = "_corrected"

else:
    def corrected_x(power_dbm, freq_ghz):
        return power_dbm

    def corrected_y(peak_db, freq_ghz, harmonic):
        return peak_db

    x_label      = "Input Power [dBm]"
    y_label      = "Peak Magnitude [dB]"
    fname_suffix = "_raw"

############################
# CORRECTION SANITY CHECK — CONSOLE
############################

print("\n--- Input correction sanity check (X-axis) ---")
print(f"{'Freq':>10}  {'PumpLine':>10}  {'PumpIn':>10}  {'CryoInOut':>11}  {'loss':>10}  {'X corr':>10}")
print("-" * 72)
for freq in frequency_list:
    f = freq * 1e9
    if APPLY_CORRECTIONS:
        pl  = get_pump_line(f)
        pi  = get_pump_in(f)
        cio = get_cryo_in_out(f)
        l   = loss(f)
        xc  = pl + pi + 0.5 * ATT_in - 0.5 * cio
    else:
        pl = pi = cio = l = xc = 0.0
    print(f"{freq:>8.3f} GHz  {pl:>+10.3f}  {pi:>+10.3f}  {cio:>+11.3f}  {l:>+10.3f}  {xc:>+10.3f}")

print()
print("--- Output correction sanity check (Y-axis) ---")
print(f"{'Freq':>8}  {'Harmonic':>10}  {'Harm. freq':>12}  {'loss(f,h)':>11}  {'Y corr':>10}")
print("-" * 60)
for freq in frequency_list:
    for harm in harmonic_list:
        f_h = freq * harm * 1e9
        if APPLY_CORRECTIONS:
            l   = loss(f_h)
            yc  = -l
        else:
            l = yc = 0.0
        print(f"{freq:>6.3f} GHz  {harm:>10}  {freq*harm:>8.3f} GHz  {l:>+11.3f}  {yc:>+10.3f}")

print("--- End of sanity check ---\n")

############################
# PLOT LOADED S2P CURVES
############################

if APPLY_CORRECTIONS:
    fig_s2p, ax_s2p = plt.subplots(figsize=(9, 5))
    ax_s2p.plot(f_pl  / 1e9, s_pl,  label=f"Pump line    ({S2P_PUMP_LINE.name})",   linewidth=1.5)
    ax_s2p.plot(f_pi  / 1e9, s_pi,  label=f"Pump in      ({S2P_PUMP_IN.name})",     linewidth=1.5)
    ax_s2p.plot(f_cio / 1e9, s_cio, label=f"Cryo in-out  ({S2P_CRYO_IN_OUT.name})", linewidth=1.5)

    ax_s2p.set_xlabel("Frequency [GHz]", fontsize=12)
    ax_s2p.set_ylabel("S Magnitude [dB]", fontsize=12)
    ax_s2p.set_title("Correction curves — loaded & smoothed S-parameters", fontsize=13)
    ax_s2p.legend(fontsize=8)
    ax_s2p.grid(True, linestyle="--", alpha=0.5)
    fig_s2p.tight_layout()

    plot_path_s2p = PLOT_SAVE_DIR / "correction_S_curves.png"
    fig_s2p.savefig(plot_path_s2p, dpi=150)
    print(f"Saved: {plot_path_s2p.name}")

############################
# LOAD MASTER LOG
############################

master_csv = SAVE_DIR / f"{FILE_PREFIX}_log.csv"
log_lookup = {}

if master_csv.exists():
    with master_csv.open("r") as f:
        reader = csv.DictReader(f)
        for row in reader:
            key = (float(row["frequency_GHz"]),
                   float(row["power_dBm"]),
                   int(row["harmonic"]))
            log_lookup[key] = {
                "peak_dB":        float(row["peak_dB"]),
                "temperature_mK": float(row["temperature_mK"]),
            }
    print(f"Loaded master log: {len(log_lookup)} entries from {master_csv.name}")
else:
    print(f"Master log not found ({master_csv}). Falling back to reading individual CSVs.")

############################
# CSV READER (fallback)
############################

def read_csv_peak(filepath):
    mags = []
    with open(filepath, "r") as f:
        reader = csv.DictReader(f)
        for row in reader:
            mags.append(float(row["magnitude_dB"]))
    return float(np.max(mags)) if mags else np.nan

############################
# BUILD DATA ARRAYS
############################

peak_data = {freq: {harm: {} for harm in harmonic_list} for freq in frequency_list}
temp_data = {freq: {harm: {} for harm in harmonic_list} for freq in frequency_list}
missing   = []

for freq in frequency_list:
    for pwr in power_list:
        for harm in harmonic_list:
            key = (freq, pwr, harm)
            if key in log_lookup:
                peak_data[freq][harm][pwr] = log_lookup[key]["peak_dB"]
                temp_data[freq][harm][pwr] = log_lookup[key]["temperature_mK"]
            else:
                fname = f"{FILE_PREFIX}_{freq:.3f}GHz_{pwr:.1f}dBm_{harm}H.csv"
                fpath = SAVE_DIR / fname
                if not fpath.exists():
                    missing.append(fname)
                    peak_data[freq][harm][pwr] = np.nan
                    temp_data[freq][harm][pwr] = np.nan
                else:
                    try:
                        peak_data[freq][harm][pwr] = read_csv_peak(fpath)
                        temp_data[freq][harm][pwr] = np.nan
                    except Exception as e:
                        print(f"  WARNING: could not read {fname}: {e}")
                        peak_data[freq][harm][pwr] = np.nan
                        temp_data[freq][harm][pwr] = np.nan

if missing:
    print(f"\nMissing files ({len(missing)}):")
    for m in missing:
        print(f"  {m}")

############################
# PLOTTING
############################

cmap   = cm.get_cmap("tab10", len(harmonic_list))
colors = [cmap(i) for i in range(len(harmonic_list))]

# ── Per-frequency figures ─────────────────────────────────────────────────────

for freq in frequency_list:
    fig1, ax1 = plt.subplots(figsize=(7, 5))

    for idx, harm in enumerate(harmonic_list):
        x_vals = [corrected_x(pwr, freq) for pwr in power_list]
        y_vals = [corrected_y(peak_data[freq][harm][pwr], freq, harm)
                  for pwr in power_list]

        ax1.plot(x_vals, y_vals,
                 marker="o", color=colors[idx], linewidth=1.5,
                 label=f"Harmonic {harm}  ({freq * harm:.2f} GHz)")

    ax1.set_xlabel(x_label, fontsize=12)
    ax1.set_ylabel(y_label, fontsize=12)
    ax1.set_title(f"Harmonic Generation — Fundamental: {freq:.3f} GHz", fontsize=13)
    ax1.legend(fontsize=9, loc="best")
    ax1.grid(True, linestyle="--", alpha=0.5)
    fig1.tight_layout()

    plot_path1 = PLOT_SAVE_DIR / f"harmonic_{freq:.3f}GHz_magnitude{fname_suffix}.png"
    fig1.savefig(plot_path1, dpi=150)
    print(f"Saved: {plot_path1.name}")

# ── Temperature time-series (fixed width) ────────────────────────────────────

ordered_temps = []
for freq in frequency_list:
    for pwr in power_list:
        for harm in harmonic_list:
            t = temp_data[freq][harm].get(pwr, np.nan)
            ordered_temps.append((freq, pwr, harm, t))

has_temp = any(not np.isnan(t) for _, _, _, t in ordered_temps)

if has_temp:
    fig2, ax2 = plt.subplots(figsize=(12, 4))   # fixed width

    x_vals = list(range(len(ordered_temps)))
    t_vals = [t for _, _, _, t in ordered_temps]

    ax2.plot(x_vals, t_vals, color="steelblue", linewidth=1.5, marker=".", markersize=4)

    step  = len(power_list) * len(harmonic_list)
    y_top = ax2.get_ylim()[1]
    for block_idx, freq in enumerate(frequency_list):
        x_line = block_idx * step
        if x_line > 0:
            ax2.axvline(x=x_line, color="tomato", linestyle="--", linewidth=1.2)
        ax2.text(x_line + 0.3, y_top,
                 f"{freq:.3f} GHz",
                 fontsize=8, color="tomato", va="top", ha="left",
                 bbox=dict(boxstyle="round,pad=0.2", fc="white", ec="tomato", alpha=0.7))

    ax2.set_xlabel("Measurement index  (freq \u2192 power \u2192 harmonic)", fontsize=12)
    ax2.set_ylabel("Temperature [mK]", fontsize=12)
    ax2.set_title("Temperature — full time series", fontsize=13)
    ax2.grid(True, linestyle="--", alpha=0.4)
    fig2.tight_layout()

    plot_path2 = PLOT_SAVE_DIR / "temperature_timeseries.png"
    fig2.savefig(plot_path2, dpi=150)
    print(f"Saved: {plot_path2.name}")
else:
    print("  (No temperature data found — skipping temperature time-series plot)")

plt.show()
print("\nAnalysis complete.")