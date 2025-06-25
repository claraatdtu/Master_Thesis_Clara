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

**3-GNU radio implementation**: This folder contains several subfolders that are used to perform step by step the real life SDR experiments with the 2 Hack RF Ones and GNU Radio described in the Thesis report. It includes .grc files and Python scripts for transmitting and receiving bitstreams, do some bandwidth measurements and sensitivity testing. When creating a .grc file, it creates automatically the main application block in the Python script (using GNU Radio and PyQt). Here is a small explanation for each subfolder:
- *SDR files of bits*: It contains the input bitstream used for every SDR test (`sdrinput`), the output folder containing the output bitstreams obtained for each scheme for specific Eb/N0, the Python script `BERcomparison_input_output.py` computing the BER for each output bitstream and the MATLAB script `plot_ber_ebn0.m` plotting the final graph of all BER obtained with the SDR experiment.
- *loopbacks working*: It contains the Python scripts and .grc files for each modulation scheme loopback flowgraph. The flowgraphs contain the modulation and demodulation of the input bitstream and the insertion of noise via the "Add" and "Noise Source" blocks. The Python scripts used are `LoRaloopback.py`, `BFSKloopback_real.py`, `QuatreFSKloopback.py`, `EightFSKloopback.py`, `SixteenFSKloopback.py`, `PSKloopback.py` and `PSKloopback_wform.py`. The used loopback flowgraphs in GNU Radio Companion are: `LoRaloopback.grc`, `BFSKloopback_real.grc`, `4FSKloopback.grc`, `8FSKloopback.grc`, `16FSKloopback.grc` and`PSKloopback.grc` (changing the constellation Object and Points for QPSK or BPSK). 
- *TX working*: It contains the scripts used to transmit the input bitstream for each scheme via the HackRF One number 1 connected to the computer. The TX scripts are obtained from the *loopbacks working* scripts by removing the receiving part and adding a Soapy Hackrf Sink: `LoRaTX.grc`, `BFSKTX.grc`, `4FSKTX.grc`, `8FSKTX.grc`, `16FSKTX.grc` and`PSKTX.grc`
- *RX working*: It contains the scripts used to receive the input bitstream for each scheme via the HackRF One number 2 connected to the computer. The TX scripts are obtained from the *loopbacks working* scripts by removing the transmitting part and adding a Soapy Hackrf Source/ `LoRaRX.grc`, `BFSKRX.grc`, `4FSKRX.grc`, `8FSKRX.grc`, `16FSKRX.grc` and`PSKRX.grc`
- *SDR bandwidths*: It contains the Python script `calculate_SDR_bw.py` used to compute the 99% bandwidth for each modulation scheme, running at different bitrates. It also gathers 2 folders containing the output frequency bins used by each scheme, for bit rates similar to LoRa CSS SF=7 and SF=12.
- *SDR evaluation of NF*: It contains 2 Python scripts/ GNU Radio: `TX_evaluation_NF.py` and `RX_evaluation_NF.py`. Both need a connection to a HackRF One. It evaluates the Noise Figure of the receiving HackRF. The objective is to introduce this measured value in the sensitivity calculation.
- *old*: This folder is just a backup of earlier test files  and flowgraphs that were no selected.

**4-pictures report**: This file gathers the pictures (pdf and png) used in the Thesis report.


## The dependencies I used (software and hardware)
- MATLAB R2024a 
- GNU Radio Companion 3.10.10.0 
- Python 3.11.9 (for GNU Radio use and also for post-processing and plotting)
- 2 HackRF One from GREAT SCOTT GADGETS (for the SDR measurements)


## Problem Formulation of the Thesis report
This project focuses on evaluating the spectrum efficiency and sensitivity of LoRa modulation
(CSS) compared to X-SK modulation. The LoRa SX1276 chip claims a sensitivity of -137 dBm
due to its FEC, for specific parameters (125kHz Bandwidth, Spreading Factor of 12 and data rate
< 300 bps) as summarized in Figure 2.4. As a comparison, a very good FSK transceiver can have
a sensitivity of -127 dBm (12kHz Bandwidth, for a similar data rate). Through mathematical
analysis and testing, this project aims to understand the contributions of FEC and modulation
to LoRa’s sensitivity and spectrum efficiency. This project’s goal is also to comprehend the
redundant data and collision avoidance contributions to sensitivity.
The key objectives of this study are to:
- Calculate the sensitivity of LoRa without FEC and evaluate FEC’s contribution to
sensitivity.
- Quantify the redundant information transmitted in a LoRa transmission.
- Develop a mathematical approach to compare fairly the spectrum efficiency and sensitivity
of LoRa and X-SK.
- Analyze the role of Chirp Spread Spectrum in collision avoidance and its impact on
sensitivity.
- Simulate the different modulation schemes transmissions on MATLAB
- Prototype the set up on SDR and perform tests to validate the theoretical findings.