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
    print(f"\nTotal length at the input: {len(input_bits)} bits")


## bit vector
# Convert bitstring to list of integers
input_bits_list = [int(bit) for bit in input_bits]

print("Bits as list:")
print(input_bits_list)
print(f"\nTotal length at the input: {len(input_bits_list)} bits")


## READ THE OUTPUT DATA
#path to the output file
file_path = r"C:\Users\clsor\OneDrive\Documents\MATLAB\Master_Thesis_Clara\Master_Thesis_Clara\3-GNU radio implementation\SDR files of bits\LoRa12sdroutput-5"

with open(file_path, 'rb') as f:
    byte_data = f.read()

   # convert the bytes to bits
try:
    output_bits = ''.join(str(byte) for byte in byte_data)  # No need to format to 8 bits
except ValueError:
    print("File contains unexpected byte values (not 0 or 1).")
else:
    print("Bits:")
    print(output_bits)
    print(f"\nTotal length at the output: {len(output_bits)} bits")




## bit vector
# Convert bitstring to list of integers
output_bits_list = [int(bit) for bit in output_bits]

print("Bits as list:")
print(output_bits_list)
print(f"\nTotal length at the output: {len(output_bits_list)} bits")




## COMPARE INPUT AND OUTPUT FOR BER
delay=0 #for lora
#delay=32 #for qpsk
#delay=16 #for bpsk
#delay= #for 2Fsk
#delay= #for 4Fsk
#delay= #for 8Fsk
#delay= #for 16Fsk

if input_bits_list is not None and output_bits_list is not None:
    #length for comparison: should be 100 000
    min_len = min(len(input_bits_list), len(output_bits_list)-delay)
    # compare the bits
    errors = sum(1 for i in range(min_len) if input_bits_list[i] != output_bits_list[i+delay])
    ber = errors / (min_len) #if min_len > 0 else None

    print("\nBER Computation") #ber results to display
    print(f"Bit errors: {errors}")
    if ber is not None:
        print(f"BER: {ber:.6f}")
    else:
        print("BER issue of computation")


print(min_len)



##
L_output = []
L_input = []
for i in range (0,80000):
    L_output.append(int(output_bits_list[i]))
    L_input.append(int(input_bits_list[i]))


print(L_input)
print(L_output)

if L_input is not None and L_output is not None:
    #length for comparison: should be 100 000
    min_len = min(len(L_input), len(L_output))
    # compare the bits
    errors = sum(1 for i in range(min_len) if L_input[i] != L_output[i])
    ber = errors / min_len #if min_len > 0 else None

    print("\nBER Computation") #ber results to display
    print(f"Bit errors: {errors}")
    if ber is not None:
        print(f"BER: {ber:.6f}")
    else:
        print("BER issue of computation")



## preamble reco
preamble = '00110011011000111'
positions = []

start = 0
while True:
    index = output_bits.find(preamble, start)
    if index == -1:
        break
    positions.append(index)
    start = index + 1  # Move forward to find next occurrence

print("Preamble found at indices:", positions)
           #Preamble found at indices: [107199, 143127]