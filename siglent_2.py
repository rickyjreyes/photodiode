import pyvisa
import numpy as np
import cupy as cp
import matplotlib.pyplot as plt
from scipy.signal import find_peaks
import csv
import re
import time

import os
from dotenv import load_dotenv

load_dotenv()

SCOPE_IP = os.environ["SCOPE_IP"]
CHANNEL = "C3"

# -----------------------------
# Connect to oscilloscope
# -----------------------------

rm = pyvisa.ResourceManager('@py')
scope = rm.open_resource(f"TCPIP::{SCOPE_IP}::INSTR")

scope.timeout = 20000

print(scope.query("*IDN?"))

# -----------------------------
# Utility parser
# -----------------------------

def parse_value(s):
    val = re.findall(r"[-+]?\d*\.\d+E[-+]?\d+|[-+]?\d*\.\d+|[-+]?\d+", s)
    return float(val[0])

# -----------------------------
# Stop acquisition
# -----------------------------

scope.write(":STOP")

# -----------------------------
# Channel scaling
# -----------------------------

vdiv = parse_value(scope.query(f"{CHANNEL}:VDIV?"))
offset = parse_value(scope.query(f"{CHANNEL}:OFST?"))

sample_rate = parse_value(scope.query(":ACQ:SRAT?"))

print("Vertical scale:", vdiv)
print("Offset:", offset)
print("Sample rate:", sample_rate)

# -----------------------------
# Configure waveform transfer
# -----------------------------

scope.write(f":WAV:SOUR {CHANNEL}")
scope.write(":WAV:FORM BYTE")

# -----------------------------
# Retrieve waveform
# -----------------------------

data = scope.query_binary_values(":WAV:DATA?", datatype='B')
wave = np.array(data)

print("Samples:", len(wave))

# -----------------------------
# Convert ADC → voltage
# -----------------------------

adc_mid = 128
volts = (wave - adc_mid) * (vdiv / 25.0) - offset

# remove DC
volts = volts - np.mean(volts)

# time axis
dt = 1.0 / sample_rate
time_axis = np.arange(len(volts)) * dt

# -----------------------------
# GPU FFT
# -----------------------------

volts_gpu = cp.asarray(volts)

fft_gpu = cp.fft.rfft(volts_gpu)

spectrum = cp.abs(fft_gpu)

freq = cp.fft.rfftfreq(len(volts_gpu), d=dt)

spectrum = cp.asnumpy(spectrum)
freq = cp.asnumpy(freq)

# -----------------------------
# Peak detection
# -----------------------------

threshold = np.max(spectrum) * 0.1
peaks, props = find_peaks(spectrum, height=threshold)

peak_freqs = freq[peaks]
peak_amps = props['peak_heights']

# sort by frequency
order = np.argsort(peak_freqs)
peak_freqs = peak_freqs[order]
peak_amps = peak_amps[order]

# -----------------------------
# Compute harmonic ratios
# -----------------------------

results = []

if len(peak_freqs) > 0:

    f1 = peak_freqs[0]
    A1 = peak_amps[0]

    print("\nDetected peaks:")

    for i in range(len(peak_freqs)):
        f = peak_freqs[i]
        A = peak_amps[i]

        ratio = f / f1
        amp_ratio = A / A1

        print(
            f"{i+1}: {f:.3f} Hz | "
            f"ratio {ratio:.3f} | "
            f"A={A:.3f} | "
            f"A/A1={amp_ratio:.3f}"
        )

        results.append([i+1, f, ratio, A, amp_ratio])

# -----------------------------
# Save results
# -----------------------------

filename = f"spectral_run_{int(time.time())}.csv"

with open(filename, "w", newline="") as f:
    writer = csv.writer(f)

    writer.writerow([
        "mode",
        "frequency_Hz",
        "harmonic_ratio",
        "amplitude",
        "amplitude_ratio"
    ])

    for row in results:
        writer.writerow(row)

print("\nSaved:", filename)

# -----------------------------
# Plot waveform
# -----------------------------

plt.figure(figsize=(8,4))
plt.plot(time_axis, volts)
plt.title("Oscilloscope Waveform (DC removed)")
plt.xlabel("Time (s)")
plt.ylabel("Voltage (V)")
plt.grid()

# -----------------------------
# Plot spectrum
# -----------------------------

plt.figure(figsize=(8,4))
plt.plot(freq, spectrum)

if len(peaks) > 0:
    plt.scatter(freq[peaks], spectrum[peaks], color="red")

plt.xlim(0, 500)
plt.title("FFT Spectrum")
plt.xlabel("Frequency (Hz)")
plt.ylabel("Amplitude")
plt.grid()

plt.show()