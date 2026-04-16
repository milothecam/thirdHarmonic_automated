import time
import pathlib
import pyvisa
import csv
import numpy as np

############################
# USER CONFIGURATION
############################

# Input lists — non-integer values are fine
frequency_list = [2, 3, 4, 5]            # GHz — fundamental tone frequencies
power_list     = np.linspace(15, 25, 51).tolist()  # dBm — signal generator output powers
harmonic_list  = [1, 3]            # harmonic index (multiplier of fundamental)

# VNA test receiver settings
VNA_RECEIVER   = "C"       # test receiver parameter — maps to T3(1) on S4VNA
VNA_SPAN_HZ    = 1e6       # Hz — fixed span around each harmonic centre frequency
VNA_IFBW_HZ    = 10e3        # Hz — IF (resolution) bandwidth
VNA_NUM_POINTS = 1001       # number of frequency points

# Timing
PAUSE_TIME = 60              # seconds — wait between measurements

# Output directory for CSV files
SAVE_DIR = pathlib.Path(r"C:\data\Camilo\Wafer UVA microstrip\microstrip E\harmonicsAutomated\microstripE_2345GHz_3rdH_10mK_fineStep")
SAVE_DIR.mkdir(parents=True, exist_ok=True)

# File prefix
FILE_PREFIX = "harmAuto"

# Lakeshore temperature log
LAKESHORE_LOG_ROOT = pathlib.Path("C:/Users/bluefors/Documents/logging/temperature")
LAKESHORE_CHANNEL  = 6

############################
# OPEN VISA INSTRUMENTS
############################

rm = pyvisa.ResourceManager()

SG  = rm.open_resource("TCPIP0::10.12.96.50::inst0::INSTR")   # Agilent E8257D
VNA = rm.open_resource("TCPIP0::127.0.0.1::5025::SOCKET")   # Copper Mountain VNA

SG.read_termination  = "\n"
SG.write_termination = "\n"
VNA.read_termination = "\n"
VNA.timeout = 60000   # ms

print("SG  ID:", SG.query("*IDN?").strip())
print("VNA ID:", VNA.query("*IDN?").strip())

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


def configure_vna_receiver():
    """
    Set up VNA in test receiver mode once before the measurement loop.
    Frequency/span/IFBW are updated per iteration inside measure_and_save_csv().
    """
    VNA.write("CALC:PAR:COUN 1")
    VNA.write(f"CALC:PAR1:DEF {VNA_RECEIVER}")
    VNA.write("CALC:PAR1:SEL")
    VNA.write("CALC:FORM MLOG")
    VNA.query("*OPC?")
    print(f"VNA configured: test receiver '{VNA_RECEIVER}'")


def measure_and_save_csv(center_freq_hz, filepath):
    """
    Update VNA centre frequency, trigger a single sweep, retrieve the
    frequency axis and magnitude trace, and save to a two-column CSV.
    Returns (freqs_hz, mag_db) as numpy arrays.
    """
    # Update frequency axis and sweep parameters
    VNA.write(f"SENS:FREQ:CENT {center_freq_hz}")
    VNA.write(f"SENS:FREQ:SPAN {VNA_SPAN_HZ}")
    VNA.write(f"SENS:SWE:POIN {VNA_NUM_POINTS}")
    VNA.write(f"SENS:BWID {VNA_IFBW_HZ}")
    VNA.query("*OPC?")

    # 
    # # Abort any running sweep and force hold before triggering
    # VNA.write("ABOR")
    # VNA.write("INIT:CONT OFF")      # stop continuous sweep — hold mode
    # VNA.query("*OPC?")
    # 

    # Single triggered sweep
    VNA.write("TRIG:SOUR BUS")
    VNA.write("TRIG:SING")
    VNA.query("*OPC?")

    # Retrieve frequency axis
    raw_freq = VNA.query("CALC:DATA:XAX?")
    freqs_hz = np.array([float(x) for x in raw_freq.strip().split(",")])

    # Retrieve trace data — for a MLOG-formatted test receiver trace, FDAT returns
    # interleaved (formatted_value, 0) pairs where the real part is already in dBm.
    # No log conversion needed — unlike S-parameter traces which return raw re/im.
    raw_data = VNA.query("CALC:DATA:FDAT?")
    pairs    = np.array([float(x) for x in raw_data.strip().split(",")])
    mag_db   = pairs[0::2]   # real part = dBm value directly

    # Save CSV
    with open(filepath, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["frequency_Hz", "magnitude_dB"])
        for freq, mag in zip(freqs_hz, mag_db):
            writer.writerow([freq, mag])

    return freqs_hz, mag_db


############################
# CONFIGURE VNA ONCE
############################

configure_vna_receiver()

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
            print(f"  Saved       : {fname}  (peak = {peak:.2f} dB)")

            # Accumulate metadata
            master_log.append({
                "filename":       fname,
                "frequency_GHz":  freq,
                "power_dBm":      pwr,
                "harmonic":       harm,
                "center_GHz":     center_hz / 1e9,
                "temperature_mK": T,
                "peak_dB":        peak,
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
              "center_GHz", "temperature_mK", "peak_dB"]

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
VNA.close()
rm.close()
print("Done.")
