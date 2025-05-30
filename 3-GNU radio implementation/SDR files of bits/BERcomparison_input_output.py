# This script opens the sdr input bitstream of 100 000 bits and read/display it. It does the same with the sdr output bitstream after reception with the hackrf. Then, it compares input and output to determine the BER. This is done for several output data files for each different SNR and modulation scheme.
#Author: Clara SORRE

import os

## READ THE INPUT DATA
#path to the input file
file_path = r"C:\Users\clsor\OneDrive\Documents\MATLAB\Master_Thesis_Clara\Master_Thesis_Clara\3-GNU radio implementation\SDR files of bits\sdrinput"

with open(file_path, 'rb') as f:
    byte_data = f.read()

   # convert the bytes to bits
try:
    input_bits = ''.join(str(b) for b in byte_data)  # No need to format to 8 bits
except ValueError:
    print("File contains unexpected byte values (not 0 or 1).")
else:
    print("Bits:")
    print(input_bits)
    print(f"\nTotal length at the input: {len(bits)} bits")


## READ THE OUTPUT DATA
#path to the input file
file_path = r"C:\Users\clsor\OneDrive\Documents\MATLAB\Master_Thesis_Clara\Master_Thesis_Clara\3-GNU radio implementation\SDR files of bits\BFSKsdroutput"

with open(file_path, 'rb') as f:
    byte_data = f.read()

   # convert the bytes to bits
try:
    output_bits = ''.join(str(b) for b in byte_data)  # No need to format to 8 bits
except ValueError:
    print("File contains unexpected byte values (not 0 or 1).")
else:
    print("Bits:")
    print(output_bits)
    print(f"\nTotal length at the output: {len(bits)} bits")



## COMPARE INPUT AND OUTPUT FOR BER
if input_bits is not None and output_bits is not None:
    #length for comparison: should be 100 000
    min_len = min(len(input_bits), len(output_bits))
    # compare the bits
    errors = sum(1 for i in range(min_len) if input_bits[i] != output_bits[i])
    ber = errors / min_len #if min_len > 0 else None

    print("\n BER Computation") #ber results
    print(f"Bit errors: {errors}")
    if ber is not None:
        print(f"BER: {ber:.6f}")
    else:
        print("BER issue of computation")



