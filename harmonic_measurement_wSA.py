import time
import pathlib
import pyvisa
import csv
import numpy as np

############################
# USER CONFIGURATION
############################

# Input lists — non-integer values are fine
# frequency_list = [2.16, 2.213, 2.626, 2.668]            # GHz — fundamental tone frequencies
frequency_list = [3.5] # np.linspace(2, 5, 7)            # GHz — fundamental tone frequencies
power_list     = np.linspace(5, 25, 21).tolist()  # dBm — signal generator output powers
harmonic_list  = [1, 3, 5]            # harmonic index (multiplier of fundamental)


# Spectrum analyzer settings
SA_SPAN_HZ    = 2e4    # Hz — fixed span around each harmonic centre frequency
SA_RBW_HZ     = 100    # Hz — resolution bandwidth (equivalent to VNA IFBW)
SA_VBW_HZ     = 100    # Hz — video bandwidth (set equal to RBW; reduce to smooth noise floor)
SA_NUM_POINTS = 501    # number of sweep points
# N9010A reference level note:
#   - Preamp OFF (recommended for signals > -40 dBm): range is roughly -50 to +30 dBm
#   - Preamp ON:  range is restricted to approx. -100 to -35 dBm
#   Set ~10 dB above your expected peak signal. Run preamp diagnostic to confirm your range.
SA_REF_LEVEL  = -40    # dBm — safe default; adjust after confirming preamp state
SA_PREAMP_ON  = False  # set True only if measuring signals below ~ -50 dBm
SA_ATTEN_DB   = 10     # dB  — input attenuation
SA_DETECTOR   = "NORM" # detector mode: NORM | AVER | PEAK | SAMP

# Timing
PAUSE_TIME = 31        # seconds — wait between measurements

# Output directory for CSV files
SAVE_DIR = pathlib.Path(r"C:\data\Camilo\Wafer Delft 1st run 2023\1b\harmonicsAutomated_1b\twoTones_3.5GHz")
SAVE_DIR.mkdir(parents=True, exist_ok=True)

# File prefix
FILE_PREFIX = "Delft1stRun_1b_10mK_SA_twoTones"

# Lakeshore temperature log
LAKESHORE_LOG_ROOT = pathlib.Path("C:/Users/bluefors/Documents/logging/temperature")
LAKESHORE_CHANNEL  = 6

############################
# OPEN VISA INSTRUMENTS
############################

rm = pyvisa.ResourceManager()

SG = rm.open_resource("TCPIP0::10.12.96.50::inst0::INSTR")   # Agilent E8257D
SA = rm.open_resource("USB0::0x2A8D::0x1B0B::MY62060378::0::INSTR")   # Keysight N9010B  ← update IP

SG.read_termination  = "\n"
SG.write_termination = "\n"
SA.read_termination  = "\n"
SA.write_termination = "\n"
SA.timeout = 60_000   # ms — generous timeout for slow sweeps

print("SG ID:", SG.query("*IDN?").strip())
print("SA ID:", SA.query("*IDN?").strip())

############################
# INSTRUMENT FUNCTIONS
############################

def read_lakeshore_temperature():
    """Read last line of latest Lakeshore log folder. Returns temperature in mK."""
    latest_folder = sorted([x for x in LAKESHORE_LOG_ROOT.iterdir()])[-1]
    log_file = latest_folder / f"CH{LAKESHORE_CHANNEL} T {latest_folder.stem}.log"
    with log_file.open() as f:
        last = f.readlines()[-1]
        temp_K = float(last.split(",")[-1])
    return temp_K * 1e3   # convert to mK


def sg_set_tone(freq_ghz, power_dbm):
    """Configure and enable the Agilent E8257D CW tone."""
    SG.write("*RST")
    SG.write(f"FREQ {freq_ghz}GHz")
    SG.write(f"POW {power_dbm}DBM")
    SG.write("OUTP ON")
    time.sleep(0.3)


def sg_off():
    """Turn off the signal generator RF output."""
    SG.write("OUTP OFF")


def configure_sa():
    """
    One-time spectrum analyzer setup before the measurement loop.

    Critical sequence for N9010B:
      1. Reset and wait.
      2. Explicitly select SA application — *RST can leave the instrument
         in a state where SCPI commands go to a background app and the
         front panel does not update.
      3. Select the swept-SA measurement (:CONF:SAN).
      4. Configure all parameters while continuous sweep is ON so the
         display refreshes and you can confirm settings visually.
      5. Switch to single-sweep mode only at the very end.
    """
    SA.write("*RST")
    SA.query("*OPC?")                               # wait for reset

    # ── Force SA application and measurement ──────────────────
    SA.write(":INST:SEL SA")                        # select Spectrum Analyzer app
    SA.query("*OPC?")
    SA.write(":CONF:SAN")                           # select swept-power SA measurement
    SA.query("*OPC?")

    # ── Leave continuous sweep ON while configuring so the ────
    # ── display updates and settings are visible on screen  ────
    SA.write(":INIT:CONT ON")

    SA.write(":UNIT:POW DBM")                       # amplitude units → dBm
    # ── Preamp — must be set BEFORE ref level ─────────────────
    # On the N9010A the preamp restricts the ref level ceiling to ~ -35 dBm.
    # Turn it off unless you specifically need it for very weak signals.
    SA.write(f":SENS:POW:GAIN:STAT {'ON' if SA_PREAMP_ON else 'OFF'}")
    SA.query("*OPC?")
    SA.write(f":SENS:DET:FUNC {SA_DETECTOR}")       # detector mode
    SA.write(f":SENS:BWID:RES {SA_RBW_HZ}")        # resolution bandwidth
    SA.write(f":SENS:BWID:VID {SA_VBW_HZ}")        # video bandwidth
    SA.write(f":SENS:POW:RF:ATT {SA_ATTEN_DB}")    # input attenuation (set before ref level)
    SA.write(f":DISP:WIND:TRAC:Y:RLEV {SA_REF_LEVEL}")  # reference level
    SA.write(f":SENS:SWE:POIN {SA_NUM_POINTS}")    # number of sweep points
    SA.write(":FORM:DATA ASC")                      # ASCII trace transfer
    SA.query("*OPC?")

    # ── Now freeze into single-sweep mode ─────────────────────
    SA.write(":INIT:CONT OFF")
    SA.query("*OPC?")

    # ── Drain error queue so we catch any bad commands above ──
    errors = []
    for _ in range(10):
        err = SA.query(":SYST:ERR?").strip()
        if err.startswith("+0") or err.startswith("0,"):
            break
        errors.append(err)
    if errors:
        print(f"  WARNING — SA reported SCPI errors during configure: {errors}")

    print(f"SA configured: RBW={SA_RBW_HZ/1e3:.0f} kHz | VBW={SA_VBW_HZ/1e3:.0f} kHz | "
          f"Atten={SA_ATTEN_DB} dB | Det={SA_DETECTOR} | Pts={SA_NUM_POINTS}")


def measure_and_save_csv(center_freq_hz, filepath):
    """
    Update SA centre frequency, trigger a single sweep, retrieve the
    frequency axis (reconstructed) and magnitude trace, and save to a
    two-column CSV.
    Returns (freqs_hz, mag_dbm) as numpy arrays.

    NOTE — The N9010B does not return a frequency axis with the trace;
    we reconstruct it linearly from start/stop readback, which matches
    the instrument's internal grid exactly.
    """
    # ── Update frequency axis ──────────────────────────────────
    SA.write(f":SENS:FREQ:CENT {center_freq_hz}")
    SA.write(f":SENS:FREQ:SPAN {SA_SPAN_HZ}")
    SA.query("*OPC?")

    # Readback start/stop so the reconstructed axis is exact
    f_start = float(SA.query(":SENS:FREQ:STAR?").strip())
    f_stop  = float(SA.query(":SENS:FREQ:STOP?").strip())

    # ── Trigger one sweep and wait ─────────────────────────────
    SA.write(":INIT:IMM")     # trigger single sweep
    SA.query("*OPC?")         # block until sweep complete

    # ── Fetch trace ────────────────────────────────────────────
    raw = SA.query(":TRAC:DATA? TRACE1").strip()
    mag_dbm = np.array([float(x) for x in raw.split(",")])

    # Reconstruct frequency axis
    freqs_hz = np.linspace(f_start, f_stop, len(mag_dbm))

    # ── Save CSV ───────────────────────────────────────────────
    with open(filepath, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["frequency_Hz", "magnitude_dBm"])
        for freq, mag in zip(freqs_hz, mag_dbm):
            writer.writerow([freq, mag])

    return freqs_hz, mag_dbm


############################
# CONFIGURE SA ONCE
############################

configure_sa()

############################
# MAIN MEASUREMENT LOOP
# Order: frequency → power → harmonic
############################

master_log = []

total = len(frequency_list) * len(power_list) * len(harmonic_list)
count = 0

for freq in frequency_list:
    for pwr in power_list:
        for harm in harmonic_list:
            count += 1
            center_hz = freq * harm * 1e9

            print(f"\n[{count}/{total}]  freq={freq} GHz | power={pwr} dBm | "
                  f"harmonic={harm}  →  centre={center_hz/1e9:.4f} GHz")

            # Read temperature BEFORE turning on the signal generator
            T = read_lakeshore_temperature()
            print(f"  Temperature : {T:.2f} mK")

            # Apply tone at the fundamental frequency
            sg_set_tone(freq, pwr)

            # Measure and save CSV
            fname = f"{FILE_PREFIX}_{freq:.3f}GHz_{pwr:.1f}dBm_{harm}H.csv"
            fpath = SAVE_DIR / fname
            freqs, mags = measure_and_save_csv(center_hz, fpath)
            peak = float(mags.max())
            print(f"  Saved       : {fname}  (peak = {peak:.2f} dBm)")

            # Accumulate metadata
            master_log.append({
                "filename":       fname,
                "frequency_GHz":  freq,
                "power_dBm":      pwr,
                "harmonic":       harm,
                "center_GHz":     center_hz / 1e9,
                "temperature_mK": T,
                "peak_dBm":       peak,
            })

            # Turn off signal generator between measurements
            sg_off()

            # Pause before next iteration
            time.sleep(PAUSE_TIME)

############################
# SAVE MASTER METADATA CSV
############################

master_csv = SAVE_DIR / f"{FILE_PREFIX}_log.csv"
fieldnames = ["filename", "frequency_GHz", "power_dBm", "harmonic",
              "center_GHz", "temperature_mK", "peak_dBm"]

with master_csv.open("w", newline="") as f:
    writer = csv.DictWriter(f, fieldnames=fieldnames)
    writer.writeheader()
    writer.writerows(master_log)

print(f"\nMaster log saved to: {master_csv}")

############################
# CLEANUP
############################

print("Closing instrument connections...")
sg_off()
SG.close()
SA.close()
rm.close()
print("Done.")