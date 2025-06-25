# Master_Thesis_Clara

## Explanation of the content of the scripts
This folder contains MATLAB and Python files to simulate and compare modulation techniques (PSK, M-FSK) and LoRa regarding their sensitivity and spectrum efficiency. Each one of the scripts is commented.

**1-Mathematical Model**: This folder contains 4 MATLAB scripts used to set up the modeling of the following modulation schemes: LoRa CSS, MFSK, GFSK and PSK. The objective is to compute and visualize:
- the BER (Bit Error Rate: P_b) vs Eb/N0
- the Q-function interpretation
- for LoRa: the SNR versus the spreading factor and coding rate (parity bits values) 
- the receiver sensitivity as a function of SNR and BER
- the sensitivity tables for different SF and parity set ups


**2-MATLAB simulations**: This folder contains 4 BER simulations in MATLAB to plot BER vs Eb/N0 for the following schemes:
- BPSK and QPSK (`test_psk_only.m`)
- M-FSK with M = 2, 4, 8, and 16 (`test_mfsk_only.m`)
- LoRa CSS modulation (`test_lora_only.m`)
- a combined simulation (`test_total_modulationsandLora.m`) that gathers all of them
For each, a random bitstream is generated, modulated, noise (AWGN) is added, and demodulated. The resulting BER is then calculated for several values of Eb/N0.


**3-GNU radio implementation**: This folder contains several subfolders that are used to perform step by step the real life SDR experiments with the 2 Hack RF Ones and GNU Radio described in the Thesis report. It includes .grc files and Python scripts for transmitting and receiving bitstreams, do some bandwidth measurements and sensitivity testing. Here is a small explanation for each subfolder:
- *SDR files of bits*: It contains the input bitsream used for every SDR test (`sdrinput`), the output folder containing the output bitstreams obtained for each scheme for specific Eb/N0, the Python script `BERcomparison_input_output` computing the BER for each output bitstream and the MATLAB script `plot_ber_ebn0.m` plotting the final graph of all BER obtained with the SDR experiment.
- *loopbacks working*:
- *TX working*:
- *RX working*:
- *SDR bandwidths*:
- *SDR evaluation of NF*:
- *old*: This folder is just a backup of earlier test files  and flowgraphs that were no selected.

**4-pictures report**: This file gathers the pictures (pdf and png) used in the Thesis report.

## The dependencies I used (software and hardware)
- MATLAB R2024a 
- GNU Radio Companion 3.10.10.0 
- Python 3.11.9 (for GNU Radio use and also for post-processing and plotting)
- 2 HackRF One from GREAT SCOTT GADGETS (for the SDR measurements)
