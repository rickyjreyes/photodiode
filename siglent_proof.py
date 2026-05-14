import pyvisa
import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import find_peaks
import csv
import time
import re

# -----------------------------
# SCOPE SETTINGS
# -----------------------------

import os
from dotenv import load_dotenv

load_dotenv()

SCOPE_IP = os.environ["SCOPE_IP"]
CHANNEL = "C3"

# -----------------------------
# CONNECT
# -----------------------------

rm = pyvisa.ResourceManager('@py')

scope = rm.open_resource(f"TCPIP::{SCOPE_IP}::INSTR")

scope.timeout = 20000

print(scope.query("*IDN?"))

# -----------------------------
# HELPER
# -----------------------------

def extract_number(txt):
    m = re.search(r'[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?', txt)
    if m:
        return float(m.group(0))
    return 0.0

# -----------------------------
# STOP ACQUISITION
# -----------------------------

scope.write(":STOP")

# -----------------------------
# CONFIGURE WAVEFORM
# -----------------------------

scope.write("WFSU SP,0,NP,0,FP,0")
scope.write("WAV:MODE RAW")
scope.write("WAV:SOURCE " + CHANNEL)

# -----------------------------
# SCALE PARAMETERS
# -----------------------------

vscale = extract_number(scope.query(f"{CHANNEL}:VDIV?"))
offset = extract_number(scope.query(f"{CHANNEL}:OFST?"))

sample_rate = extract_number(scope.query("SARA?"))

print("Vertical scale:", vscale)
print("Offset:", offset)
print("Sample rate:", sample_rate)

# -----------------------------
# READ DATA
# -----------------------------

scope.write("WAV:DATA?")
raw = scope.read_raw()

# strip header
header_len = int(raw[1:2])
data = raw[2 + header_len:]

wave = np.frombuffer(data, dtype=np.int8)

samples = len(wave)

print("Samples:", samples)

# -----------------------------
# VOLTAGE CONVERSION
# -----------------------------

voltage = wave * vscale / 25 - offset

# remove DC
voltage -= np.mean(voltage)

# -----------------------------
# TIME AXIS
# -----------------------------

t = np.arange(samples) / sample_rate

# -----------------------------
# FFT
# -----------------------------

fft = np.fft.rfft(voltage)
freq = np.fft.rfftfreq(samples, 1/sample_rate)

amp = np.abs(fft)

# -----------------------------
# PEAK DETECTION
# -----------------------------

noise_floor = np.median(amp)

peaks, _ = find_peaks(
    amp,
    height=noise_floor * 10,
    distance=20
)

peak_freqs = freq[peaks]
peak_amps = amp[peaks]

order = np.argsort(peak_freqs)

peak_freqs = peak_freqs[order]
peak_amps = peak_amps[order]

# -----------------------------
# HARMONIC ANALYSIS
# -----------------------------

print("\nDetected peaks:")

if len(peak_freqs) > 0:

    f0 = peak_freqs[0]

    for i,(f,a) in enumerate(zip(peak_freqs,peak_amps)):

        ratio = f/f0

        print(
            f"{i+1}: {f:.3f} Hz | ratio {ratio:.3f} | A={a:.3f}"
        )

# -----------------------------
# SAVE CSV
# -----------------------------

timestamp = int(time.time())

csv_name = f"spectral_proof_{timestamp}.csv"

with open(csv_name,"w",newline="") as f:

    writer = csv.writer(f)

    writer.writerow(["harmonic","frequency","ratio","amplitude"])

    for i,(f,a) in enumerate(zip(peak_freqs,peak_amps)):

        writer.writerow([
            i+1,
            f,
            f/peak_freqs[0] if len(peak_freqs) else 0,
            a
        ])

print("\nSaved:",csv_name)

# -----------------------------
# PLOTS
# -----------------------------

plt.figure()
plt.plot(t, voltage)
plt.title("Oscilloscope Waveform (DC removed)")
plt.xlabel("Time (s)")
plt.ylabel("Voltage (V)")
plt.grid()

plt.figure()
plt.plot(freq, amp)

if len(peaks) > 0:
    plt.scatter(peak_freqs, peak_amps, color='red')

plt.title("FFT Spectrum")
plt.xlabel("Frequency (Hz)")
plt.ylabel("Amplitude")
plt.grid()

plt.show()