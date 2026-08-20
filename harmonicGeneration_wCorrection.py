import math
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import skrf as rf  # scikit-rf for s2p files

def plot_first_and_third_curves(csv_lists, x_values,
                                s2p_file, s2p_pump, s2p_longLine,
                                target_freq, attenuation,
                                labels=None,
                                colors=None,
                                markers=None,
                                linestyles=None,
                                sizes=None,
                                alphas=None,
                                fit_points=3,
                                title="First & Third Harmonic vs Input",
                                xlabel="Input power to TKIP [dBm]",
                                ylabel="Output power from TKIP [dBm]"):
    """
    Plot first and third harmonic curves using adjusted x-axis corrected by two S2P files.

    New correction:
      x_adj = np.array(x_values) + S12_pump(target_freq) - S12_longLine(target_freq)

    Parameters
    ----------
    s2p_file : str
        Main device S2P file for harmonic loss extraction
    s2p_pump : str
        S2P file for pump line
    s2p_longLine : str
        S2P file for long line
    """

    # --- Load s2p file and extract losses ---
    ntwk = rf.Network(s2p_file)
    idx = np.argmin(np.abs(ntwk.f - target_freq))
    loss = 20 * math.log10(np.abs(ntwk.s[idx, 0, 1]))

    idx3 = np.argmin(np.abs(ntwk.f - 3 * target_freq))
    loss_3h = 20 * math.log10(np.abs(ntwk.s[idx3, 0, 1]))

    # --- Load correction S2Ps ---
    pump = rf.Network(s2p_pump)
    longLine = rf.Network(s2p_longLine)

    idx_pump = np.argmin(np.abs(pump.f - target_freq))
    s12_pump = 20 * math.log10(np.abs(pump.s[idx_pump, 1, 0]))  # S12 magnitude in dB

    idx_long = np.argmin(np.abs(longLine.f - target_freq))
    s12_long = 20 * math.log10(np.abs(longLine.s[idx_long, 1, 0]))  # S12 magnitude in dB

    # --- Apply correction ---
    correction = s12_pump - s12_long + loss + 3
    x_adj = np.array(x_values) + correction
    print("Adjusted x-axis values:", x_adj)

    print(f"S12_cryostat({target_freq/1e9:.3f} GHz) = {loss:.3f} dB")
    print(f"S21_cryostat @ {3*target_freq/1e9:.3f} GHz = {loss_3h:.3f} dB")
    print(f"S12_pump({target_freq/1e9:.3f} GHz) = {s12_pump:.3f} dB")
    print(f"S12_longLine({target_freq/1e9:.3f} GHz) = {s12_long:.3f} dB")
    print(f"Applied correction to x-axis: Δ = {correction:.3f} dB")

    plt.figure(figsize=(10, 6))

    # defaults
    if labels is None:
        labels = ["Tone at input", "3rd Harmonic at output"]
    if colors is None:
        colors = ["blue", "magenta"]
    if markers is None:
        markers = ["o", "d"]
    if linestyles is None:
        linestyles = ["-", "-"]
    if sizes is None:
        sizes = [60, 60]
    if alphas is None:
        alphas = [0.8, 0.8]

    # Loop only over curve 0 and curve 2
    for idx_curve, i in enumerate([0, 2]):
        file_group = csv_lists[i]

        # max of second column
        y_values = []
        for file in file_group:
            data = pd.read_csv(file)
            col2 = pd.to_numeric(data.iloc[:, 1], errors="coerce")
            y_values.append(col2.dropna().max())
        y_values = np.array(y_values)

        # y-value corrections
        if i == 0:  # fundamental
            y_values = y_values - loss + attenuation + 3
            print("Adjusted tone values:", y_values)
        elif i == 2:  # 3rd harmonic
            y_values = y_values - loss_3h + attenuation + 3
            print("Adjusted 3rd harmonic values:", y_values)

        # plot data
        plt.plot(x_adj, y_values,
                 marker=markers[idx_curve],
                 linestyle=linestyles[idx_curve],
                 color=colors[idx_curve],
                 markersize=np.sqrt(sizes[idx_curve]),
                 alpha=alphas[idx_curve])

        # linear fit (first N points)
        if fit_points > len(x_adj):
            raise ValueError("fit_points is greater than the number of data points.")
        coeffs = np.polyfit(x_adj[:fit_points], y_values[:fit_points], 1)
        slope, intercept = coeffs
        fit_line = np.polyval(coeffs, x_adj)

        plt.plot(x_adj, fit_line,
                 linestyle="--",
                 color=colors[idx_curve],
                 alpha=0.7)

        # put fit equation in legend
        plt.scatter([], [], marker=markers[idx_curve],
                    color=colors[idx_curve],
                    label=f"{labels[idx_curve]} (fit: y={slope:.2f}x+{intercept:.2f})")

    plt.title(title)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.grid(True)
    plt.legend()
    plt.tight_layout()
    plt.show()

    


# Example usage
csv_lists = [
    ["5GHz_T_15dBm.csv",  "5GHz_T_17dBm.csv",  "5GHz_T_19dBm.csv",  "5GHz_T_21dBm.csv",  "5GHz_T_23dBm.csv",  "5GHz_T_25dBm.csv"],
    ["5GHz_2H_15dBm.csv", "5GHz_2H_17dBm.csv", "5GHz_2H_19dBm.csv", "5GHz_2H_21dBm.csv", "5GHz_2H_23dBm.csv", "5GHz_2H_25dBm.csv"],
    ["5GHz_3H_15dBm.csv", "5GHz_3H_17dBm.csv", "5GHz_3H_19dBm.csv", "5GHz_3H_21dBm.csv", "5GHz_3H_23dBm.csv", "5GHz_3H_25dBm.csv"]
]

s2p_file = "microstripE_noLNA_nomKcoupler_10mK_0to20GHz.s2p"
s2p_pump = "pumpPath_withVeryLongCable_S21.s2p"
s2p_longLine = "veryLongCable_S21.s2p"

target_freq = 5e9
x_values = [15, 17, 19, 21, 23, 25]
attenuation = -20

plot_first_and_third_curves(csv_lists, x_values,
                             s2p_file, s2p_pump, s2p_longLine,
                             target_freq, attenuation,
                             labels=["Tone at TKIP input", "3rd Harmonic at TKIP output"],
                             colors=[(.3, 0, .3), "magenta"],
                             markers=["o", "d"],
                             linestyles=["-", "-"],
                             sizes=[60, 60],
                             alphas=[0.8, 0.8],
                             fit_points=6,
                             title="Output per Harmonic Component vs Input Power @ 5 GHz")
