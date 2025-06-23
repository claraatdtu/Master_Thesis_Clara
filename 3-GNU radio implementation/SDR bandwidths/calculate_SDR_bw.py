# This script opens the sdr input frames and read/display them. It computes the FFT power spectrum. This is done for several data files for each modulation scheme.
#Author: Clara SORRE

import numpy as np
import matplotlib.pyplot as plt
fft_size = 1024
sample_rate = 1e6  # 1 MHz
center_freq = 868.1e6  # Hz
# FFT log power spectrum in dB from gnu
data = np.fromfile(r"C:\Users\clsor\OneDrive\Documents\MATLAB\Master_Thesis_Clara\Master_Thesis_Clara\3-GNU radio implementation\SDR bandwidths\BPSKbw", dtype=np.float32)
frames = data.reshape(-1, fft_size) # reshape if multiple frames
freqs = np.linspace(center_freq - sample_rate / 2, center_freq + sample_rate / 2, fft_size) / 1e6  # frequency axis in MHz
bin_width = sample_rate / fft_size # width of each frequency bin

# compute 99% BW for each frame and average
bw_list = []
for frame in frames:
    vec_lin = 10 ** (frame / 10)
    total_power = np.sum(vec_lin)
    cumsum = np.cumsum(vec_lin)
    low = np.searchsorted(cumsum, 0.005* total_power)
    high = np.searchsorted(cumsum, 0.995 * total_power)
    bw = (high - low) * bin_width
    bw_list.append(bw)
bw_99_avg = np.mean(bw_list)
print(f"Average 99% Occupied Bandwidth: {bw_99_avg:.2f} Hz")
# plot last frame for visualization
vec_db = frames[-1]
plt.plot(freqs, vec_db)
plt.title("FFT Power Spectrum for BPSK (last frame, dB)")
plt.xlabel("Frequency (MHz)")
plt.ylabel("Power (dB)")
plt.grid(True)
plt.show()

# vec_db = frames[-1] #pick a frame
# vec_lin = 10 ** (vec_db / 10)
# total_power = np.sum(vec_lin) #normalize the energy
# cumsum = np.cumsum(vec_lin)
# # indices that contain 99% of energy
# low = np.searchsorted(cumsum, 0.005 * total_power)
# high = np.searchsorted(cumsum, 0.995 * total_power)
#
# bw_99 = (high - low) * bin_width # occupied bandwidth
# print(f"The 99% occupied bandwidth is {bw_99:.2f} Hz")
#
#
# # plot
# plt.plot(freqs, vec_db)  #frequency axis
# plt.title("FFT Power Spectrum for 8FSK(dB)")
# plt.xlabel("Frequency (MHz)")
# plt.ylabel("Power (dB)")
# plt.grid(True)
# plt.show()
