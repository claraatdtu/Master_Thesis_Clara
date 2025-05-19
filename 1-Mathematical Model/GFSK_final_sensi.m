clc; clear; close all;
%% GFSK bit rate = same of lora for comparing
SF_values = 12; %or 7
p_values = 0:4;
BW_LoRa = 125e3;
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

%% assumptions: change of BW but rest fixed
%BW = 125e3;
Rb= bit_rate(1,1);
NF= 6;
h = 0.5; %usual modulation index

% %% PSK bit rate = same of lora for comparing
% SF_values = 12;
% p_values = 0:4;
% 
% % initialize 
% spec_eff = zeros(length(SF_values), length(p_values));
% 
% bit_rate = zeros(length(SF_values), length(p_values));
% % spectral efficiency values
% for i = 1:length(SF_values)
%     SF = SF_values(i);
%     for j = 1:length(p_values)
%         p = p_values(j);
%         spec_eff(i, j) = (SF / (2^SF)) * (4 / (4 + p));
%         bit_rate(i, j)=spec_eff(i, j).*BW;
%     end
% end


%% GFSK FUNCTION
function Pb = compute_Pb(EbN0, h)
    % Constants
    alpha = EbN0;  % Compute alpha
    
    % Compute a and b
    a = sqrt(alpha / 2) .* sqrt(1 - sqrt(1 - (sin(2 * pi * h) ./ (2 * pi * h)).^2));
    b = sqrt(alpha / 2) .* sqrt(1 + sqrt(1 - (sin(2 * pi * h) ./ (2 * pi * h)).^2));
    
    % Compute the summation term
    Pb = exp(-alpha / 2) .* (0.5 * besseli(0, a .* b)); % First term in brackets
    
    % Sum over k terms (truncated at k_max for approximation)
    k_max = 1000; % infinite series Truncated
    for k = 1:k_max
        Pb = Pb + ((a ./ b).^k) .* besseli(k, a .* b);
    end
end


%% use function iterations

EbN0_dB = -5:15;  %SNR values in dB
%bit_rate = [366];
%figure; 
%hold on;

Pb_values = compute_Pb(10.^(EbN0_dB/10), h); %convert dB to linear
semilogy(EbN0_dB, Pb_values,'DisplayName', sprintf('GFSK'), 'LineWidth', 1.5, 'Marker', 'o')

for idx = 1:length(p_values)
    %interpolation
    p = p_values(idx);
    target_Pb = 1e-3;
    EbN0_estimated = interp1(Pb_values, EbN0_dB, target_Pb, 'linear', 'extrap');
    Sensitivity_dBm = -174 + 10 * log10(bit_rate(1, p+1) * 10.^(EbN0_estimated/10)) + NF;
    fprintf('For GFSK at p= %d,', p)
    fprintf(' the estimated Eb/N0 is %d dB', EbN0_estimated)
    %fprintf(' Estimated GFSK SNR at P_b = 1e-5 is %.2f dB', EbN0_estimated);
    fprintf(' and the sensitivity is %d dBm. \n', Sensitivity_dBm)
    % for idx = 1:length(bit_rate)
    %     R_b = bit_rate(idx);
    % end
end
title('GFSK P_b as a function of Eb/N0')
xlabel('Eb/N0 (dB)');
ylabel('P_b (Average Bit Error Probability)');
legend show;
grid on;
hold off;

