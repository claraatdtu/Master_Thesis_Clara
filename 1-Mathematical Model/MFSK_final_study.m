%% MFSK sensitivity/spectrum efficiency study
% AUTHOR: Clara SORRE
% This MATLAB code analyzes the performance of MFSK modulation defined in
% the report, in terms of BER, sensitivity and spectral efficiency.
clc; clear; close all;

%% 1- compute Pb in function of EbN0: function reused in the following
function Pb = compute_Pb(EbN0_vec, M)
    Pb = zeros(size(EbN0_vec)); %preallocate the vector
    for k = 1:length(EbN0_vec)
        EbN0 = EbN0_vec(k);
        Pb_k = 0;
        for n = 1:M-1
            term = ((M / 2) / (M - 1)) * ((-1)^(n+1) / (n + 1)) * nchoosek(M - 1, n) * exp((-n * log2(M) * EbN0) / (n + 1)); %log2(M) because Es=m*Eb
            Pb_k = Pb_k + term; %each term of the approximation formula
        end
        Pb(k) = Pb_k;
    end
end

%% 2- Main script: MFSK P_b as a function of E_b/N_0
target_Pb = 1e-3; %target ber
m_values = 1:5; %this corresponds to M from 2 to 32
marker_list = {'o', 's', 'd', '^', 'v'}; % the markers 
EbN0_dB_range = -5:1:15;
EbN0_lin_range = 10.^(EbN0_dB_range / 10); %converting to linear scale
EbN0_estimated = zeros(1, length(m_values));
figure;
for idx = 1:length(m_values)
    m = m_values(idx);
    M = 2^m;
    Pb_values = compute_Pb(EbN0_lin_range, M); %interpolate to find the required Eb/N0 for the target ber
    EbN0_estimated(idx) = interp1(Pb_values, EbN0_dB_range, target_Pb, 'linear', 'extrap'); 
    semilogy(EbN0_dB_range, Pb_values, 'LineWidth', 1.5, 'Marker', marker_list{idx},'DisplayName', sprintf('%d-FSK', M)); hold on;
end
title('MFSK P_b as a function of E_b/N_0');
xlabel('E_b/N_0 (dB)');
ylabel('P_b (Average Bit Error Probability)');
legend show;
grid on;
ylim([1e-7, 1]);

%% 3- LoRa spectral efficiency & bit rate (SF fixed)
SF_values = 7:12;
p_values = 0:4;
BW_LoRa = 125e3;
spec_eff = zeros(length(SF_values), length(p_values));
bit_rate = zeros(length(SF_values), length(p_values));
for i = 1:length(SF_values)
    SF = SF_values(i);
    for j = 1:length(p_values)
        p = p_values(j);
        spec_eff(i,j) = (SF / (2^SF)) * (4 / (4 + p));
        bit_rate(i,j) = spec_eff(i,j) * BW_LoRa;
    end
end
disp(bit_rate);

%% 4- Sensitivity fixed like LoRa CSS: Bit rate reachable and bandwidth needed for each modulation scheme at fixed Pb
BW_LoRa=125e3;
NF= 6; %the receiver noise figure in dB
Pb_target = 1e-3; % target probability of bit error
z = sqrt(2) * erfcinv(4 * Pb_target); % corresponding z value
SF_values = 7:12;
p_values = 0:4;
Sensitivity_LoRa = zeros(length(SF_values), length(p_values)); % initialize result matrix
for i = 1:length(SF_values)%recreate the sensitivity of LoRa
    SF = SF_values(i);
    for j = 1:length(p_values)
        p = p_values(j);
        SNR_linear = ((z * sqrt(2) * SF) / (log10(SF)/log10(12))).* (4 / (4 + p)) ./ (2.^SF);% compute SNR in linear scale
        SNR_dB = 10 * log10(SNR_linear);
        Sensitivity_dBm = -174 + 10 * log10(BW_LoRa) + NF + SNR_dB;
        Sensitivity_LoRa(i, j) = Sensitivity_dBm;  % store in table
    end
end
Rb_table = zeros(length(SF_values), length(p_values), length(m_values));
BW_table = zeros(length(SF_values), length(p_values), length(m_values));
for idx1 = 1:length(m_values)
    m = m_values(idx1);
    M = 2^m;
    EbN0_lin = 10^(EbN0_estimated(idx1)/10);
    for idx2 = 1:length(SF_values)
        for idx3 = 1:length(p_values)
            p = p_values(idx3);
            S = Sensitivity_LoRa(idx2,idx3); % fixed sensitivity like LoRa
            Rb = ((10^((S + 174 - NF)/10)) / EbN0_lin)*(4/(p+4)); % Rb from noise power and sensitivity
            BW = (M * Rb) / m;  % required BW for M-FSK
            Rb_table(idx2, idx3, idx1) = Rb;
            BW_table(idx2, idx3, idx1) = BW;
        end
    end
    fprintf('Achievable Bit Rate for %d-FSK (bps): \n', M);
    disp(Rb_table(:,:,idx1));
    fprintf('and Needed Bandwidth for %d-FSK (Hz): \n', M);
    disp(BW_table(:,:,idx1));
end
% heatmaps for Bit Rate and Bandwidth at Sensitivity Fixed (p = 0)
SF_values = 7:12;
m_values = [1 2 3 4 5]; % → M = 2, 4, 8, 16, 32
M_labels = {'2-FSK', '4-FSK', '8-FSK', '16-FSK', '32-FSK'};
p_index = 1; % p = 0
Rb_p0 = squeeze(Rb_table(:, p_index, :));  % [SF x MFSK]for p = 0
BW_p0 = squeeze(BW_table(:, p_index, :));  % [SF x MFSK]
bluemap = [linspace(0.8,0,64)', linspace(0.9,0.4,64)', ones(64,1)];
greenmap = [linspace(0.9,0.2,64)', ones(64,1), linspace(0.9,0.4,64)'];
figure; % Heatmap: Bit Rate
imagesc(Rb_p0);
colormap(greenmap);
colorbar;
title('Bit Rate (bps) for Parity = 0');
xlabel('Modulation');
ylabel('Spreading Factor (SF)');
xticks(1:length(M_labels));
xticklabels(M_labels);
yticks(1:length(SF_values));
yticklabels(arrayfun(@(x) sprintf('SF%d', x), SF_values, 'UniformOutput', false));
for i = 1:size(Rb_p0,1)
    for j = 1:size(Rb_p0,2)
        text(j, i, sprintf('%.0f', Rb_p0(i,j)), 'HorizontalAlignment', 'center', 'Color', 'k');
    end
end
% heatmap: Bandwidth
figure;
imagesc(BW_p0);
colormap(bluemap);
colorbar;
title('Bandwidth (Hz) for Parity = 0');
xlabel('Modulation');
ylabel('Spreading Factor (SF)');
xticks(1:length(M_labels));
xticklabels(M_labels);
yticks(1:length(SF_values));
yticklabels(arrayfun(@(x) sprintf('SF%d', x), SF_values, 'UniformOutput', false));
for i = 1:size(BW_p0,1)
    for j = 1:size(BW_p0,2)
        text(j, i, sprintf('%.0f', BW_p0(i,j)), 'HorizontalAlignment', 'center', 'Color', 'k');
    end
end

%% 5- Bit rate fixed like LoRa CSS: Sensitivity reachable and bandwidth needed for each modulation scheme at fixed Pb
NF = 6; % noise figure
S = zeros(length(SF_values), length(p_values), length(m_values));
BW = zeros(length(SF_values), length(p_values), length(m_values));
for idx1 = 1:length(m_values) %5
    m = m_values(idx1);
    fprintf('Estimated Sensitivity for %d-FSK (dB): \n', 2^m);
    for idx2 = 1:length(SF_values) %7
        for idx3 = 1:length(p_values) %5
            EbN0_lin = 10^(EbN0_estimated(idx1)/10);
            Rb = bit_rate(idx2, idx3);
            S(idx2, idx3, idx1) = -174 + NF + 10 * log10(Rb * EbN0_lin);
            BW(idx2, idx3, idx1) = ((2^m) * Rb) / m;
        end
    end
    disp(S(:,:,idx1));
    fprintf('and Estimated Bandwidth for %d-FSK (dB): \n', 2^m);
    disp(BW(:,:,idx1));
end
SF_values = 7:12;
m_values = [1 2 3 4 5]; % → M = 2, 4, 8, 16, 32
M_labels = {'2-FSK', '4-FSK', '8-FSK', '16-FSK', '32-FSK'}; % display for bit rate fixed
p_index = 1; % Parity = 0
S_p0 = squeeze(S(:, p_index, :));  % [SF x MFSK] % extract sensitivity and bandwidth for parity = 0
BW_p0 = squeeze(BW(:, p_index, :)); % [SF x MFSK]
% heatmap:sensitivity
figure;
imagesc(S_p0);
redmap = [ones(64,1), linspace(0.8,0,64)', linspace(0.8,0,64)'];
colormap(flipud(redmap));
colorbar;
title('Sensitivity Achievable (dBm) for Parity = 0');
xlabel('Modulation');
ylabel('Spreading Factor (SF)');
xticks(1:length(M_labels));
xticklabels(M_labels);
yticks(1:length(SF_values));
yticklabels(arrayfun(@(x) sprintf('SF%d', x), SF_values, 'UniformOutput', false));
for i = 1:size(S_p0,1)
    for j = 1:size(S_p0,2)
        text(j, i, sprintf('%.1f', S_p0(i,j)), 'HorizontalAlignment', 'center', 'Color', 'k');
    end
end
bluemap = [linspace(0.8,0,64)', linspace(0.9,0.4,64)', ones(64,1)];
figure; % heatmap: Bandwidth
imagesc(BW_p0);
colormap(bluemap);
colorbar;
title('Bandwidth Needed (Hz) for Parity = 0');
xlabel('Modulation');
ylabel('Spreading Factor (SF)');
xticks(1:length(M_labels));
xticklabels(M_labels);
yticks(1:length(SF_values));
yticklabels(arrayfun(@(x) sprintf('SF%d', x), SF_values, 'UniformOutput', false));
for i = 1:size(BW_p0,1)
    for j = 1:size(BW_p0,2)
        text(j, i, sprintf('%.0f', BW_p0(i,j)), 'HorizontalAlignment', 'center', 'Color', 'k');
    end
end


%% 6- BW fixed to 125e3 Hz: compute achievable bit rate and sensi for fixed bw
BW=125e3;
m=0;
NF= 6;
Rb = zeros(1, length(m_values));
sensi = zeros(1, length(m_values));
for idx2 = 1:length(m_values)
    m = m_values(idx2);
    Rb(idx2)=(m*BW)/(2^m);
    sensi(idx2) = -174 + NF + 10*log10(Rb(idx2) * 10^(EbN0_estimated(idx2)/10));

end
display(Rb)
display(sensi)




%% 7- OFDM-MFSK STUDY for sensi fixed like LoRa
% use bit rate from LoRa and bit rate obtained from MFSK study when
% sensitivity fixed as lora
% formula: N>= M* int(bit_rate_lora/bit_rate_mfsk)
subcarriers_spacing = zeros(length(SF_values), length(p_values), length(m_values));
bandwidth_needed_persub= zeros(length(SF_values), length(p_values), length(m_values));
can_achieve_lora_bit_rate = false(length(SF_values), length(p_values), length(m_values));
BW_used=zeros(length(SF_values), length(p_values), length(m_values));
subcarriers_needed=zeros(length(SF_values), length(p_values), length(m_values));
g=0;
for x1 = 1:length(m_values)
    m = m_values(x1);
    for x2 = 1:length(p_values)
        p = p_values(x2);
        for x3 = 1:length(SF_values)
            SF = SF_values(x3);
            subcarriers_needed(x3, x2, x1)=(2^m)*ceil(bit_rate(x3, x2)/Rb_table(x3, x2, x1));
            if subcarriers_needed(x3, x2, x1)==0
                BW_used(x3, x2, x1)= BW_table(x3, x2, x1);
            end

            if subcarriers_needed(x3, x2, x1) ~= 0
                subcarriers_spacing(x3, x2, x1) = BW_LoRa/subcarriers_needed(x3, x2, x1);
                bandwidth_needed_persub(x3, x2, x1)=(Rb_table(x3, x2, x1) *2^m)/m;
                if subcarriers_spacing(x3, x2, x1)>=bandwidth_needed_persub(x3, x2, x1)
                    can_achieve_lora_bit_rate(x3, x2, x1)= true;
                    if can_achieve_lora_bit_rate(x3, x2, x1) == true
                        BW_used(x3, x2, x1)=subcarriers_needed(x3, x2, x1)*bandwidth_needed_persub(x3, x2, x1);
                    else
                        BW_used(x3, x2, x1)=subcarriers_needed(x3, x2, x1)*bandwidth_needed_persub(x3, x2, x1);
                    end

                else
                    BW_used(x3, x2, x1)=subcarriers_needed(x3, x2, x1)*bandwidth_needed_persub(x3, x2, x1);
                end    
            
            end
        
        end

    end
end
disp(BW_used)
disp(subcarriers_needed)

