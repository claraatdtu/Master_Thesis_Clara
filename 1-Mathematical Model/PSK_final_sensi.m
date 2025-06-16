%% PSK sensitivity/spectrum efficiency study
% AUTHOR: Clara SORRE
%This MATLAB code analyzes the performance of PSK modulation defined in
% the report, in terms of BER, sensitivity and spectral efficiency.
clc; clear; close all;

%% 1- Assumptions: change of BW but rest fixed
marker_list = {'o', 's', 'd', '^', 'v'};
BW_LoRa = 125e3;
% BW_FSK_1=12e3;
% BW_FSK_2=1.2e3;
NF= 6;

%% 2- PSK bit rate = same of lora for comparing
SF_values = 7:12;
p_values = 0:4;
% initialize 
spec_eff = zeros(length(SF_values), length(p_values));
bit_rate = zeros(length(SF_values), length(p_values));
% spectral efficiency values
for i = 1:length(SF_values)
    SF = SF_values(i);
    for j = 1:length(p_values)
        p = p_values(j);
        spec_eff(i, j) = (SF / (2^SF)) * (4 / (4 + p));
        bit_rate(i, j)=spec_eff(i, j).*BW_LoRa;
    end
end

%% 3- BPSK and QPSK P_b as a function of Eb/N0
EbN0_dB = -5:1:20; % Eb/N0 in dB
EbN0_lin = 10.^(EbN0_dB / 10); % convert to linear scale
Pb_BPSK = 0.5 * erfc(sqrt(EbN0_lin)); % compute Pb for BPSK and QPSK using the same formula
figure;
semilogy(EbN0_dB, Pb_BPSK, 'b-', 'DisplayName', 'BPSK/QPSK', 'LineWidth', 1.5,'Marker', 'o');
grid on;
xlabel('E_b/N_0 (dB)');
ylabel('Bit Error Rate (P_b)');
title('BPSK and QPSK P_b as a function of Eb/N0');
legend show;
ylim([1e-7 1]);
xlim([-5 20]);


%% 4- sensitivity plot: BPSK and QPSK with conditions similar at LoRa
Pb_values = linspace(10e-10, 1.00e-3  , 30); 
p_values = 0:4; % parity bits
SF_values = 7:12;
figure; 
hold on;
for i = 1:length(SF_values)%recreate the sensitivity of LoRa
    SF = SF_values(i);
    for j = 1:length(p_values)
        p = p_values(j);
        z_primeprime = sqrt(2) * erfcinv(2*Pb_values);
        EbN0_linear =   ((z_primeprime).^2) / 2;%  linear scale
        Sensitivity_dBm = -174 + 10 * log10(bit_rate(i, 1) * EbN0_linear) + NF;
    end
    plot(Pb_values, Sensitivity_dBm, 'LineWidth', 1.5,'DisplayName', sprintf('BPSK and QPSK with Bit rate equal to LoRa, for SF = %d', SF),  'Marker', marker_list{j});
end
xlabel('P_{b}');
ylabel('Sensitivity (dBm)');
title('Sensitivity_{dBm} of BPSK and QPSK as a function of P_{b} for similar conditions as LoRa CSS different parity bits');
grid on;
legend show;
hold off;

%% 5- SNR function of z' for different p for BPSK and QPSK (non plotted in the analysis)
z_primeprime = linspace(0, 111, 1000); % start from 0
p_values = 0:4; % parity bits
k_values = 1:2; % bits per symbol
figure; 
hold on;
for p_idx = 1:length(p_values)
    p = p_values(p_idx);
    SNR_linear = spec_eff(1, p_idx) .* ((z_primeprime).^2) / 2;
    SNR_dB = 10 * log10(SNR_linear);
    plot(z_primeprime, SNR_dB, 'LineWidth', 1.5,'DisplayName', sprintf('QPSK and BPSK for conditions similar at LoRa, for p = %d', p));
end
xlabel(' z_{doubleprime} ');
ylabel('SNR (dB)');
title('SNR as a function of z_{doubleprime} for BPSK and QPSK');
grid on;
legend show;
hold off;

%% 6- z' function of Pb (non plotted in the analysis)
Pb_values = linspace(0, 10e-6  , 1000000); 
k_values = 1:2; % bits per symbol, fixed to 12 for similar study
p_values = 0:4;
figure; 
hold on;
for idx = 1:length(k_values)
    k = k_values(idx);
    z_primeprime = sqrt(2) * erfcinv(2*Pb_values);
    if k ==1
        plot(Pb_values, z_primeprime, 'LineWidth', 1.5,'DisplayName', sprintf('BPSK'));
    end
    if k ==2
        plot(Pb_values, z_primeprime, 'LineWidth', 1.5,'DisplayName', sprintf('QPSK'));
    end
end
xlabel('P_b');
ylabel('z_{primeprime}');
title('z_{primeprime} as a function of P_b');
grid on;
legend show;
hold off;

%% 7- snr function of pb
Pb_values = linspace(0, 10e-6 , 1000000); 
k_values = 1:2; % k % bits per symbol, fixed to 12 for similar study
p_values = 0:4; % parity bits
figure; 
hold on;
for idx = 1:length(p_values)
    p = p_values(idx);
    z_primeprime = sqrt(2) * erfcinv(2*Pb_values);
    SNR_linear = spec_eff(1, idx) .* ((z_primeprime).^2) / 2;
    SNR_dB = 10 * log10(SNR_linear);
    plot(Pb_values, SNR_dB, 'LineWidth', 1.5,'DisplayName', sprintf('BPSK and QPSK for conditions similar of LoRa, for p = %d', p));
end
xlabel('P_{b}');
ylabel('SNR (dB)');
title('SNR_{dB} as a function of P_{b} for cases similar to LoRa CSS for different parity bits');
grid on;
legend show;
hold off;