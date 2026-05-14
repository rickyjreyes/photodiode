import pyvisa
import numpy as np
import cupy as cp
import matplotlib.pyplot as plt
from scipy.signal import find_peaks

import os
from dotenv import load_dotenv

load_dotenv()

SCOPE_IP = os.environ["SCOPE_IP"]
CHANNEL = "C3"

# -----------------------------
# connect to oscilloscope
# -----------------------------

rm = pyvisa.ResourceManager('@py')
scope = rm.open_resource(f"TCPIP::{SCOPE_IP}::INSTR")

scope.timeout = 20000

print(scope.query("*IDN?"))

# -----------------------------
# helper for parsing siglent responses
# -----------------------------

def parse_value(cmd):
    resp = scope.query(cmd)
    val = resp.split()[1]
    val = val.replace("V","")
    return float(val)

# -----------------------------
# stop acquisition
# -----------------------------

scope.write(":STOP")

# -----------------------------
# read channel scaling
# -----------------------------

vscale = parse_value("C3:VDIV?")
voffset = parse_value("C3:OFST?")
sample_rate = float(scope.query(":ACQ:SRAT?"))

print("Vertical scale:", vscale)
print("Offset:", voffset)
print("Sample rate:", sample_rate)

# -----------------------------
# waveform capture
# -----------------------------

scope.write(":WAV:SOUR C3")
scope.write(":WAV:FORM BYTE")

data = scope.query_binary_values(":WAV:DATA?", datatype='B')

wave = np.array(data)

print("Samples:", len(wave))

# -----------------------------
# convert ADC -> voltage
# -----------------------------

adc_center = 128
volts = (wave - adc_center) * (vscale / 25) + voffset

# remove DC component
volts = volts - np.mean(volts)

# -----------------------------
# build time axis
# -----------------------------

dt = 1 / sample_rate
time = np.arange(len(volts)) * dt

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
# detect peaks
# -----------------------------

threshold = np.max(spectrum) * 0.1
peaks,_ = find_peaks(spectrum, height=threshold)

peak_freqs = freq[peaks]

print("\nDetected peaks:")

for f in peak_freqs:
    print(f"{f:.3f} Hz")

# -----------------------------
# harmonic ratios
# -----------------------------

if len(peak_freqs) > 1:

    f0 = peak_freqs[0]
    ratios = peak_freqs / f0

    print("\nHarmonic ratios:")

    for r in ratios:
        print(f"{r:.3f}")

# -----------------------------
# waveform plot
# -----------------------------

plt.figure()
plt.plot(time, volts)
plt.title("Oscilloscope Waveform (DC removed)")
plt.xlabel("Time (s)")
plt.ylabel("Voltage (V)")
plt.grid()

# -----------------------------
# FFT plot
# -----------------------------

mask = freq < 500   # focus on low frequency region

plt.figure()
plt.plot(freq[mask], spectrum[mask])
plt.scatter(freq[peaks], spectrum[peaks], color='red')
plt.title("FFT Spectrum")
plt.xlabel("Frequency (Hz)")
plt.ylabel("Amplitude")
plt.grid()

plt.show()