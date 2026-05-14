import pyvisa
import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import find_peaks, spectrogram, hilbert, correlate
import csv
import time
import re
from pathlib import Path

# =============================
# USER SETTINGS
# =============================

import os
from dotenv import load_dotenv

load_dotenv()

SCOPE_IP = os.environ["SCOPE_IP"]
CHANNEL = "C3"
SAVE_DIR = Path(".")
MAX_HARMONICS = 30
FUNDAMENTAL_SEARCH_MIN_HZ = 40
FUNDAMENTAL_SEARCH_MAX_HZ = 500
PLOT_MAX_HZ = 4000
SAVE_RAW_WAVEFORM_CSV = True
SAVE_HARMONIC_SUMMARY_CSV = True

# =============================
# HELPERS
# =============================

def extract_number(txt: str) -> float:
    m = re.search(r"[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?", txt)
    return float(m.group(0)) if m else 0.0

def save_waveform_csv(path: Path, time_axis: np.ndarray, volts: np.ndarray) -> None:
    with open(path, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["time_s", "voltage_v"])
        for t, v in zip(time_axis, volts):
            writer.writerow([t, v])

def save_harmonic_csv(path: Path, rows) -> None:
    with open(path, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["harmonic", "frequency_hz", "ratio", "amplitude", "amplitude_ratio"])
        writer.writerows(rows)

def robust_peak_pick(freq: np.ndarray, amp: np.ndarray, fmin: float, fmax: float):
    """
    Pick candidate spectral peaks in a range.
    """
    band = (freq >= fmin) & (freq <= fmax)
    f_band = freq[band]
    a_band = amp[band]

    if len(f_band) < 10:
        return np.array([], dtype=int), {}

    noise_floor = np.median(a_band)
    threshold = max(noise_floor * 8.0, np.max(a_band) * 0.05)

    peaks, props = find_peaks(
        a_band,
        height=threshold,
        distance=5
    )

    # Convert back to full-spectrum indices
    full_idx = np.where(band)[0][peaks]
    return full_idx, props

def nearest_bin(freq: np.ndarray, target: float) -> int:
    return int(np.argmin(np.abs(freq - target)))

# =============================
# CONNECT
# =============================

rm = pyvisa.ResourceManager("@py")
scope = rm.open_resource(f"TCPIP::{SCOPE_IP}::INSTR")
scope.timeout = 20000

print("Connected to:")
print(scope.query("*IDN?").strip())

# =============================
# ACQUISITION SETUP
# =============================

# Stop acquisition so the data is stable while reading
scope.write(":STOP")

# Try common Siglent waveform settings
# RAW mode gives full memory when available
try:
    scope.write("WFSU SP,0,NP,0,FP,0")
except Exception:
    pass

try:
    scope.write("WAV:MODE RAW")
except Exception:
    pass

try:
    scope.write(f"WAV:SOURCE {CHANNEL}")
except Exception:
    try:
        scope.write(f":WAV:SOUR {CHANNEL}")
    except Exception:
        pass

# =============================
# SCALE PARAMETERS
# =============================

# Different Siglent models use slightly different query names
sample_rate = 0.0
for q in ["SARA?", ":ACQ:SRAT?"]:
    try:
        sample_rate = extract_number(scope.query(q))
        if sample_rate > 0:
            break
    except Exception:
        pass

vscale = 0.0
for q in [f"{CHANNEL}:VDIV?", f":{CHANNEL}:VDIV?"]:
    try:
        vscale = extract_number(scope.query(q))
        if vscale > 0:
            break
    except Exception:
        pass

offset = 0.0
for q in [f"{CHANNEL}:OFST?", f":{CHANNEL}:OFST?"]:
    try:
        offset = extract_number(scope.query(q))
        break
    except Exception:
        pass

print(f"Sample rate: {sample_rate}")
print(f"Vertical scale: {vscale}")
print(f"Offset: {offset}")

if sample_rate <= 0:
    raise RuntimeError("Could not read sample rate from scope.")

# =============================
# READ RAW DATA
# =============================

raw = None

# Method 1: raw block read
try:
    scope.write("WAV:DATA?")
    raw = scope.read_raw()
except Exception:
    raw = None

voltage = None

if raw is not None and len(raw) > 10 and raw[:1] == b"#":
    # SCPI definite-length block: #<n><len><data>
    header_digits = int(raw[1:2].decode(errors="ignore"))
    data_len = int(raw[2:2 + header_digits].decode(errors="ignore"))
    data_start = 2 + header_digits
    data_end = data_start + data_len
    data = raw[data_start:data_end]

    wave = np.frombuffer(data, dtype=np.int8)

    if len(wave) == 0:
        raise RuntimeError("Scope returned empty waveform block.")

    # This scale was already working in your earlier script
    voltage = wave * vscale / 25.0 - offset

else:
    # Method 2: byte query fallback
    try:
        scope.write(":WAV:FORM BYTE")
        data = scope.query_binary_values(":WAV:DATA?", datatype="b")
        wave = np.array(data, dtype=np.int8)

        if len(wave) == 0:
            raise RuntimeError("Scope returned empty waveform data.")

        voltage = wave * vscale / 25.0 - offset
    except Exception as e:
        raise RuntimeError(f"Could not retrieve waveform data from scope: {e}")

# =============================
# PREP SIGNAL
# =============================

samples = len(voltage)
time_axis = np.arange(samples) / sample_rate

# Remove DC
voltage = voltage - np.mean(voltage)

print(f"Samples: {samples}")
print(f"Duration: {samples / sample_rate:.6f} s")
print(f"FFT bin spacing: {sample_rate / samples:.6f} Hz")

# =============================
# FFT
# =============================

fft_vals = np.fft.rfft(voltage)
freq = np.fft.rfftfreq(samples, d=1.0 / sample_rate)
amp = np.abs(fft_vals)

# Ignore DC for peak search
search_peaks, props = robust_peak_pick(
    freq, amp,
    FUNDAMENTAL_SEARCH_MIN_HZ,
    FUNDAMENTAL_SEARCH_MAX_HZ
)

if len(search_peaks) == 0:
    print("No strong peaks found in the search band.")
    f0 = None
    harmonic_rows = []
else:
    # Fundamental = largest peak in the search band
    best_local = search_peaks[np.argmax(amp[search_peaks])]
    f0 = freq[best_local]
    a0 = amp[best_local]

    print("\nDetected fundamental candidate:")
    print(f"f0 = {f0:.6f} Hz, A0 = {a0:.6f}")

    harmonic_rows = []
    print("\nHarmonics:")
    for n in range(1, MAX_HARMONICS + 1):
        target = n * f0
        if target > freq[-1]:
            break

        idx = nearest_bin(freq, target)
        fn = freq[idx]
        an = amp[idx]
        ratio = fn / f0 if f0 != 0 else np.nan
        amp_ratio = an / a0 if a0 != 0 else np.nan

        harmonic_rows.append([n, fn, ratio, an, amp_ratio])
        print(
            f"{n:2d}: {fn:10.6f} Hz | ratio {ratio:8.4f} | "
            f"A={an:12.6f} | A/A0={amp_ratio:10.6f}"
        )

# =============================
# EXTRA ANALYSIS
# =============================

# SNR estimate around the fundamental
snr = np.nan
if f0 is not None:
    signal_idx = nearest_bin(freq, f0)
    noise_band = (freq > 2000) & (freq < min(10000, freq[-1]))
    if np.any(noise_band):
        noise_floor = np.mean(amp[noise_band])
        snr = amp[signal_idx] / noise_floor if noise_floor > 0 else np.nan
        print(f"\nEstimated SNR: {snr:.3f}")

# Instantaneous phase
analytic = hilbert(voltage)
phase = np.unwrap(np.angle(analytic))

# Spectrogram
f_spec, t_spec, Sxx = spectrogram(voltage, fs=sample_rate)

# =============================
# SAVE FILES
# =============================

timestamp = int(time.time())

if SAVE_RAW_WAVEFORM_CSV:
    wf_path = SAVE_DIR / f"waveform_run_{timestamp}.csv"
    save_waveform_csv(wf_path, time_axis, voltage)
    print(f"Saved waveform: {wf_path}")

if SAVE_HARMONIC_SUMMARY_CSV and len(harmonic_rows) > 0:
    harm_path = SAVE_DIR / f"harmonic_summary_{timestamp}.csv"
    save_harmonic_csv(harm_path, harmonic_rows)
    print(f"Saved harmonic summary: {harm_path}")

# =============================
# PLOTS
# =============================

plt.figure(figsize=(10, 4))
plt.plot(time_axis, voltage, lw=1)
plt.title("Oscilloscope Waveform (DC removed)")
plt.xlabel("Time (s)")
plt.ylabel("Voltage (V)")
plt.grid(True)
plt.tight_layout()

plt.figure(figsize=(10, 5))
plt.plot(freq, amp, lw=1)
if f0 is not None:
    plt.axvline(f0, linestyle="--")
plt.xlim(0, PLOT_MAX_HZ)
plt.title("FFT Spectrum")
plt.xlabel("Frequency (Hz)")
plt.ylabel("Amplitude")
plt.grid(True)
plt.tight_layout()

if len(harmonic_rows) > 0:
    harm_nums = [r[0] for r in harmonic_rows]
    harm_amps = [r[3] for r in harmonic_rows]

    plt.figure(figsize=(8, 4))
    plt.plot(harm_nums, harm_amps, marker="o")
    plt.title("Harmonic Amplitudes")
    plt.xlabel("Harmonic Number")
    plt.ylabel("Amplitude")
    plt.grid(True)
    plt.tight_layout()

plt.figure(figsize=(10, 5))
plt.pcolormesh(t_spec, f_spec, Sxx, shading="gouraud")
plt.ylim(0, min(2000, f_spec[-1]))
plt.title("Spectrogram")
plt.xlabel("Time (s)")
plt.ylabel("Frequency (Hz)")
plt.colorbar(label="Power")
plt.tight_layout()

plt.figure(figsize=(10, 4))
plt.plot(phase, lw=1)
plt.title("Instantaneous Phase")
plt.xlabel("Sample")
plt.ylabel("Phase (rad)")
plt.grid(True)
plt.tight_layout()

plt.show()