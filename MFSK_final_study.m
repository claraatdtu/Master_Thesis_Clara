clc; clear; close all;

%% Vectorized!!!!!!!!!!!!!!!!!
function Pb = compute_Pb(EbN0_vec, M)
    Pb = zeros(size(EbN0_vec));
    for k = 1:length(EbN0_vec)
        EbN0 = EbN0_vec(k);
        Pb_k = 0;
        for n = 1:M-1
            term = ((M / 2) / (M - 1)) * ((-1)^(n+1) / (n + 1)) * nchoosek(M - 1, n) * exp((-n * log2(M) * EbN0) / (n + 1)); %log2(M) because Es=m*Eb
            Pb_k = Pb_k + term;
        end
        Pb(k) = Pb_k;
    end
end

%% Main script
target_Pb = 1e-5;
m_values = 1:5;
EbN0_dB_range = 0:0.1:20;
EbN0_lin_range = 10.^(EbN0_dB_range / 10);
EbN0_estimated = zeros(1, length(m_values));

figure;
for idx = 1:length(m_values)
    m = m_values(idx);
    M = 2^m;
    Pb_values = compute_Pb(EbN0_lin_range, M);
    
    % Estimate required Eb/N0 for target Pb
    EbN0_estimated(idx) = interp1(Pb_values, EbN0_dB_range, target_Pb, 'linear', 'extrap');
    
    semilogy(EbN0_dB_range, Pb_values, 'DisplayName', sprintf('%d-FSK', M)); hold on;
end

title('MFSK P_b as a function of E_b/N_0');
xlabel('E_b/N_0 (dB)');
ylabel('P_b (Average Bit Error Probability)');
legend show;
grid on;
ylim([1e-7, 1]);

%% LoRa spectral efficiency & bit rate (SF fixed)
%SF = 7;  % can change to 7 if needed
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
%% Sensitivity and Bandwidth computation: bit rate fixed
NF = 6; % Noise figure
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


%% Display for bit rate fixed
% Settings
SF_values = 7:12;
m_values = [1 2 3 4 5]; % → M = 2, 4, 8, 16, 32
M_labels = {'2-FSK', '4-FSK', '8-FSK', '16-FSK', '32-FSK'};
p_index = 1; % Parity = 0

% Extract sensitivity and bandwidth for parity = 0
S_p0 = squeeze(S(:, p_index, :));  % [SF x MFSK]
BW_p0 = squeeze(BW(:, p_index, :)); % [SF x MFSK]

% --- Heatmap: Sensitivity ---
figure;
imagesc(S_p0);
colormap(flipud(hot));
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

% --- Heatmap: Bandwidth ---
figure;
imagesc(BW_p0);
colormap(parula);
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




%% Sensitivity fixed like LoRa SF12

%% Final Sensitivity Table: Sensitivity (dBm) for each (SF, p) at fixed Pb
BW_LoRa=125e3;
NF= 6;
Pb_target = 1e-5; % Target probability of bit error
z = sqrt(2) * erfcinv(4 * Pb_target); % Corresponding z value

SF_values = 7:12;
p_values = 0:4;

Sensitivity_LoRa = zeros(length(SF_values), length(p_values)); % Initialize result matrix

for i = 1:length(SF_values)
    SF = SF_values(i);
    for j = 1:length(p_values)
        p = p_values(j);

        % Compute SNR in linear scale
        SNR_linear = ((z + 1.28 * sqrt(SF) - 0.4) / 1.28)^2 * (4 / (4 + p)) / (2^SF);

        % Convert to dB
        SNR_dB = 10 * log10(SNR_linear);

        % Compute Sensitivity
        Sensitivity_dBm = -174 + 10 * log10(BW_LoRa) + NF + SNR_dB;

        % Store in table
        Sensitivity_LoRa(i, j) = Sensitivity_dBm;
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
            S = Sensitivity_LoRa(idx2,idx3); % Fixed sensitivity like LoRa
        % Rb from noise power and sensitivity
            Rb = ((10^((S + 174 - NF)/10)) / EbN0_lin)*(4/(p+4));
    %     Rb(idx2) = (10^((S+174-NF)/10))./EbN0(idx2);
% 
%     BW(idx2) = ((2^m)*Rb(idx2))/m;
            % Required BW for M-FSK
            BW = (M * Rb) / m;
    
            Rb_table(idx2, idx3, idx1) = Rb;
            BW_table(idx2, idx3, idx1) = BW;
        end
    end
    fprintf('Achievable Bit Rate for %d-FSK (bps): \n', M);
    disp(Rb_table(:,:,idx1));
    fprintf('and Needed Bandwidth for %d-FSK (Hz): \n', M);
    disp(BW_table(:,:,idx1));
end



% 
% Rb = zeros(1, length(m_values));
% BW = zeros(1, length(m_values));
% for idx2 = 1:length(m_values)
%     m = m_values(idx2);
%     Rb(idx2) = (10^((S+174-NF)/10))./EbN0(idx2);
% 
%     BW(idx2) = ((2^m)*Rb(idx2))/m;
% 
% end
% display(Rb)
% display(BW)
%% display heat map

% Plotting heatmaps for Rb and BW
% Setup
SF_values = 7:12;
m_values = [1 2 3 4 5]; % → 2^m = 2, 4, 8, 16, 32
M_labels = {'2-FSK', '4-FSK', '8-FSK', '16-FSK', '32-FSK'};
p_index = 1; % Parity bit = 0

% Extract data for parity = 0
Rb_fixedP = squeeze(Rb_table(:, p_index, :));   % [SF x MFSK]
BW_fixedP = squeeze(BW_table(:, p_index, :)); % [SF x MFSK]

% === Sensitivity Heatmap ===
figure;
imagesc(Rb_fixedP);
colormap(flipud(hot)); % Better contrast for negative dBm values
colorbar;
title('Bit Rate achievable(bps) at Parity = 0');
xlabel('Modulation');
ylabel('Spreading Factor (SF)');
xticks(1:length(M_labels));
xticklabels(M_labels);
yticks(1:length(SF_values));
yticklabels(arrayfun(@(x) sprintf('SF%d', x), SF_values, 'UniformOutput', false));
for i = 1:size(Rb_fixedP,1)
    for j = 1:size(Rb_fixedP,2)
        text(j, i, sprintf('%.1f', Rb_fixedP(i,j)), 'HorizontalAlignment', 'center', 'Color', 'k');
    end
end

% === Bandwidth Heatmap ===
figure;
imagesc(BW_fixedP);
colormap(parula);
colorbar;
title('Bandwidth needed (Hz) at Parity = 0');
xlabel('Modulation');
ylabel('Spreading Factor (SF)');
xticks(1:length(M_labels));
xticklabels(M_labels);
yticks(1:length(SF_values));
yticklabels(arrayfun(@(x) sprintf('SF%d', x), SF_values, 'UniformOutput', false));
for i = 1:size(BW_fixedP,1)
    for j = 1:size(BW_fixedP,2)
        text(j, i, sprintf('%.0f', BW_fixedP(i,j)), 'HorizontalAlignment', 'center', 'Color', 'k');
    end
end




%% OFDM-MFSK STUDY for sensitivity fixed like LoRa
division = zeros(length(SF_values), length(p_values), length(m_values));
subcarriers_needed = zeros(length(SF_values), length(p_values), length(m_values));
subcarriers_spacing = zeros(length(SF_values), length(p_values), length(m_values));
bandwidth_needed= zeros(length(SF_values), length(p_values), length(m_values));
can_achieve_lora_bit_rate = false(length(SF_values), length(p_values), length(m_values));
BW_used=zeros(length(SF_values), length(p_values), length(m_values));
remaining_bw = zeros(length(SF_values), length(p_values), length(m_values));
spec_eff_mfsk_improved= zeros(length(SF_values), length(p_values), length(m_values));
for x1 = 1:length(m_values)
    m = m_values(x1);
    for x2 = 1:length(p_values)
        p = p_values(x2);
        for x3 = 1:length(SF_values)
            SF = SF_values(x3);
            division(x3, x2, x1)=bit_rate(x3, x2)/Rb_table(x3, x2, x1);
            subcarriers_needed(x3, x2, x1) = ceil(division(x3, x2, x1)*M);
            subcarriers_spacing(x3, x2, x1) = BW_LoRa/subcarriers_needed(x3, x2, x1);
            bandwidth_needed(x3, x2, x1)=(Rb_table(x3, x2, x1) *2^m)/m;
            if subcarriers_spacing(x3, x2, x1)>=bandwidth_needed(x3, x2, x1)
                can_achieve_lora_bit_rate(x3, x2, x1)= true;
                if can_achieve_lora_bit_rate(x3, x2, x1) == true
                    BW_used(x3, x2, x1)=subcarriers_needed(x3, x2, x1)*(Rb_table(x3, x2, x1) *2^m)/m;
                    remaining_bw(x3, x2, x1)=BW_LoRa-BW_used(x3, x2, x1);
                    spec_eff_mfsk_improved(x3, x2, x1)=Rb_table(x3, x2, x1)/BW_used(x3, x2, x1);
                end

            end

        
        end

    end
end
disp(spec_eff_mfsk_improved)
%disp(division)
%disp(subcarriers_needed)
%disp(subcarriers_spacing)
%disp(can_achieve_lora_bit_rate(:, 1, :))
%disp(BW_used)
%disp(remaining_bw(:, 1, :))
% calculations 

%% display data

% Define your data
SF_values = 7:12;
MFSK_labels = {'BFSK','4-FSK','8-FSK','16-FSK','32-FSK'};

% Extract and scale the data for parity = 0
bw_data = squeeze(remaining_bw(:, 1, :)) / 1e3;  % size: 6 SF × 5 MFSK
cmap = [ ...
    0.2, 0.4, 0.8;  % Blue
    0.4, 0.3, 0.9;  % Indigo
    0.6, 0.2, 0.8;  % Purple
    0.7, 0.3, 0.6;  % Plum
    0.8, 0.4, 0.7]; % Violet

% Plot grouped bar chart
figure;
b = bar(bw_data, 'grouped');
for k = 1:length(b)
    b(k).FaceColor = cmap(k, :);
end
set(gca, 'XTickLabel', SF_values);
xlabel('Spreading Factor (SF)');
ylabel('Remaining Bandwidth [kHz]');
legend(MFSK_labels, 'Location', 'northwest');
title('Remaining Bandwidth vs MFSK (Parity = 0)');
grid on;
% Add numeric labels on top of bars
hold on;
[rows, cols] = size(bw_data);
for i = 1:cols
    for j = 1:rows
        x = j + (i - (cols+1)/2)*(0.8/cols); % x-offset for grouped bars
        y = bw_data(j, i);
        if y > 0
            text(x, y + 3, sprintf('%.2f', y), 'HorizontalAlignment', 'center', 'FontSize', 8);
        end
    end
end




%% OFDM-MFSK STUDY for Rb fixed like LoRa
%obtained sensitivity S and BW needed BW
%use sensitivity formula!







%% BW fixed to 125e3 Hz
BW=125e3;
m=0;
NF= 6;
Rb = zeros(1, length(m_values));
sensi = zeros(1, length(m_values));
for idx2 = 1:length(m_values)
    m = m_values(idx2);
    Rb(idx2)=(m*BW)/(2^m);
    sensi(idx2) = -174 + NF + 10*log10(Rb(idx2) * 10^(EbN0(idx2)/10));

end
display(Rb)
display(sensi)
%%
clc; clear; close all;

%%
clc; clear; close all;

%% Bit Rate AND BW fixed like LoRa SF12
Rb=366; %bps
BW =125e3;
m=0;
NF= 6;

target_Pb = 1e-5;
m_values = [1, 2, 3, 4, 5];
EbN0_dB_range = -5:0.1:20;
EbN0_lin_range = 10.^(EbN0_dB_range / 10);
SNR_lin_range = EbN0_lin_range.*Rb/BW;
SNR_dB_range = 10*log10(SNR_lin_range);

EbN0= zeros(1, length(m_values));
figure;
for idx = 1:length(m_values)
    m = m_values(idx);
    M = 2^m;
    Pb_values = compute_Pb(EbN0_lin_range, M);
    
    EbN0_estimated = interp1(Pb_values, EbN0_dB_range, target_Pb, 'linear', 'extrap');
    sensi = -174 + NF + 10*log10(BW)+ SNR_dB_range;
    semilogy(EbN0_dB_range, sensi, 'DisplayName', sprintf('%d-FSK', M)); hold on;

    fprintf('For %d-FSK, the estimated Eb/N0 is %.2f dB \n', M, EbN0_estimated);
end
EbN0
title('MFSK P_b as a function of Eb/N_0');
xlabel('E_b/N_0 (dB)');
ylabel('P_b (Average Bit Error Probability)');
legend show;
grid on;
ylim([1e-7, 1]);

%%
Rb = 366; % bps
BW = 125e3; % Hz
NF = 6; % dB

m_values = [1, 2, 3, 4, 5, 6]; % Corresponds to M = 2, 4, 8, 16, 32
EbN0_dB_range = -5:0.1:20;
EbN0_lin_range = 10.^(EbN0_dB_range / 10);



for idx = 1:length(m_values)
    m = m_values(idx);
    M = 2^m;

    % Compute BER (P_b) vs. Eb/N0
    Pb_values = compute_Pb(EbN0_lin_range, M);

    % Calculate SNR for each Eb/N0
    SNR_lin = EbN0_lin_range * Rb / BW;
    SNR_dB = 10*log10(SNR_lin);

    % Calculate sensitivity (in dBm) for each point
    sensitivity_dBm = -174 + 10*log10(BW) + NF + SNR_dB;

    % Plot Sensitivity vs. BER (Pb)
    plot(Pb_values, sensitivity_dBm, 'DisplayName', sprintf('%d-FSK', M), 'LineWidth', 1.5);
    hold on;
end

grid on;
xlabel('Bit Error Rate (P_b)');
ylabel('Sensitivity (dBm)');
title('Sensitivity vs. BER for different M-FSK modulations');
legend show;
xlim([0 1e-5]);
ylim([-140 -120]);
hold off;