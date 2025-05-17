clc; clear; close all;

%% assumptions
BW=125e3;
NF= 6;
marker_list = {'o', 's', 'd', '^', 'v', '>'}; % the markers for each SF
%% z function of Pb
% Define Pb range from 10^-321 to 1 (log scale)
Pb_values = linspace(10e-100, 10e-6 , 1000000); 

% Compute z from Pb
z_values = sqrt(2) * erfcinv(2 * Pb_values);

% Plot z vs. Pb
figure;
plot(Pb_values, z_values, 'b', 'LineWidth', 1.5);
xlabel('P_b');
ylabel('z');
title('z as a function of P_b');
grid on;

%% SNR function of z
z_values = linspace(0, 111, 1000); % Start from 0

p = 0; % Parity bits
SF_values = 7:12; % Spreading Factor

figure; 
hold on;
colors = lines(length(SF_values));

for idx = 1:length(SF_values)
    SF = SF_values(idx);
    
    % Compute SNR in linear scale
    SNR_linear = ((z_values * sqrt(2) * SF) / (log10(SF)/log10(12))).* (4 / (4 + p)) ./ (2.^SF);
    
    % Convert to dB
    SNR_dB = 10 * log10(SNR_linear);
    
    % Plot SNR vs. z
    plot(z_values, SNR_dB, 'Color', colors(idx, :), 'LineWidth', 1.5, 'DisplayName', sprintf('SF = %d', SF), 'Marker', marker_list{idx}, ...
        'MarkerIndices', 1:50:length(z_values));
end

xlabel('z');
ylabel('SNR (dB)');
title(['SNR as a function of z for different SF values for p=', num2str(p)]);
grid on;
legend show;
hold off;


%% SNR function of z
z_values = linspace(0, 111, 1000); % Start from 0

SF = 12; % Parity bits
p_values = 0:4; % Spreading Factor

figure; 
hold on;
colors = lines(length(p_values));

for idx = 1:length(p_values)
    p = p_values(idx);
    
    % Compute SNR in linear scale
    SNR_linear = ((z_values * sqrt(2) * SF) / (log10(SF)/log10(12))).* (4 / (4 + p)) ./ (2.^SF);
    
    % Convert to dB
    SNR_dB = 10 * log10(SNR_linear);
    
    % Plot SNR vs. z
    plot(z_values, SNR_dB, 'Color', colors(idx, :), 'LineWidth', 1.5, 'DisplayName', sprintf('p = %d', p), 'Marker', marker_list{idx},'MarkerIndices', 1:50:length(z_values));
end

xlabel('z');
ylabel('SNR (dB)');
title(['SNR as a function of z for different p values for SF=', num2str(SF)]);
grid on;
legend show;
hold off;

%% snr function of pb

Pb_values = linspace(0, 10e-3, 1000); 

% Compute z from Pb
z_values = sqrt(2) * erfcinv(2* Pb_values);

p = 0; % Parity bits
SF_values = 7:12; % Spreading Factor

figure; 
hold on;
colors = lines(length(SF_values));

for idx = 1:length(SF_values)
    SF = SF_values(idx);
    
    % Compute SNR in linear scale
    SNR_linear = ((z_values * sqrt(2) * SF) / (log10(SF)/log10(12))).* (4 / (4 + p)) ./ (2.^SF);
    
    
    % Convert to dB
    SNR_dB = 10 * log10(SNR_linear);
    
    % Plot SNR vs. z
    plot(Pb_values, SNR_dB, 'Color', colors(idx, :), 'LineWidth', 1.5, 'DisplayName', sprintf('SF = %d', SF),'Marker', marker_list{idx},'MarkerIndices', 1:50:length(z_values));
end

xlabel('P_{b}');
ylabel('SNR (dB)');
title(['SNR_{dB} as a function of P_{b} for different SF values for p=', num2str(p)]);
grid on;
legend show;
hold off;

%% snr function of pb for different p values, at sf=12

Pb_values = linspace(0, 10e-3, 1000); 

% Compute z from Pb
z_values = sqrt(2) * erfcinv(2* Pb_values);

SF = 12; % Parity bits
p_values = 0:4; % Spreading Factor


figure; 
hold on;
colors = lines(length(p_values));

for idx = 1:length(p_values)
    p = p_values(idx);
    
    % Compute SNR in linear scale
    SNR_linear = ((z_values * sqrt(2) * SF) / (log10(SF)/log10(12))).* (4 / (4 + p)) ./ (2.^SF);
    
    % Convert to dB
    SNR_dB = 10 * log10(SNR_linear);
    
    % Plot SNR vs. z
    plot(Pb_values, SNR_dB, 'Color', colors(idx, :), 'LineWidth', 1.5, 'DisplayName', sprintf('p = %d', p), 'Marker', marker_list{idx},'MarkerIndices', 1:50:length(z_values));
end

xlabel('P_{b}');
ylabel('SNR (dB)');
title(['SNR_{dB} as a function of P_{b} for different p values for SF=', num2str(SF)]);
grid on;
legend show;
hold off;
%% sensitivity sf varies
Pb_values = linspace(10e-10, 10e-4  , 1000);  %0.25

% Compute z from Pb
z_values = sqrt(2) * erfcinv(2 * Pb_values);

p = 0; % Parity bits
SF_values = 7:12; % Spreading Factor

figure; 
hold on;
colors = lines(length(SF_values));

for idx = 1:length(SF_values)
    SF = SF_values(idx);
    
    % Compute SNR in linear scale
    SNR_linear = ((z_values * sqrt(2) * SF) / (log10(SF)/log10(12))).* (4 / (4 + p)) ./ (2.^SF);
    
    % Convert to dB
    SNR_dB = 10 * log10(SNR_linear);
    Sensitivity_dBm = -174 + 10 * log10(BW) + NF + SNR_dB;
    % Plot SNR vs. z
    plot(Pb_values, Sensitivity_dBm, 'Color', colors(idx, :), 'LineWidth', 1.5, 'DisplayName', sprintf('SF = %d', SF), 'Marker', marker_list{idx},'MarkerIndices', 1:50:length(z_values));
end

xlabel('P_{b}');
ylabel('Sensitivity (dBm)');
title(['Sensitivity_{dB} as a function of P_{b} for different SF values for p=', num2str(p)]);
grid on;
legend show;
hold off;

%% sensitivity p varies
Pb_values = linspace(10e-10, 10e-4, 1000); % Probability of bit error

% Compute z from Pb
z_values = sqrt(2) * erfcinv(2 * Pb_values);

SF = 12; % Fixed Spreading Factor
p_values = 0:4; % Parity bits

figure;
hold on;
colors = lines(length(p_values));

for idx = 1:length(p_values)
    p = p_values(idx);
    
    % Compute SNR in linear scale
    SNR_linear = ((z_values * sqrt(2) * SF) / (log10(SF)/log10(12))).* (4 / (4 + p)) ./ (2.^SF);
    
    % Convert to dB
    SNR_dB = 10 * log10(SNR_linear);
    Sensitivity_dBm = -174 + 10 * log10(BW) + NF + SNR_dB;
    
    % Plot Sensitivity vs. Pb
    plot(Pb_values, Sensitivity_dBm, 'Color', colors(idx, :), 'LineWidth', 1.5, 'DisplayName', sprintf('p = %d', p), 'Marker', marker_list{idx},'MarkerIndices', 1:50:length(z_values));
end

xlabel('P_{b}');
ylabel('Sensitivity (dBm)');
title(['Sensitivity_{dB} as a function of P_{b} for SF = ', num2str(SF)]);
grid on;
legend show;
hold off;


%% TABLE SENSITIVITY FINAL

%% Final Sensitivity Table: Sensitivity (dBm) for each (SF, p) at fixed Pb
BW=125e3;
NF= 6;
Pb_target = 1e-3; % Target probability of bit error
z = sqrt(2) * erfcinv(2 * Pb_target); % Corresponding z value

SF_values = 7:12;
p_values = 0:4;

Sensitivity_table = zeros(length(SF_values), length(p_values)); % Initialize result matrix

for i = 1:length(SF_values)
    SF = SF_values(i);
    for j = 1:length(p_values)
        p = p_values(j);

        % Compute SNR in linear scale
        SNR_linear = ((z * sqrt(2) * SF) / (log10(SF)/log10(12))).* (4 / (4 + p)) ./ (2.^SF);
    
        % Convert to dB
        SNR_dB = 10 * log10(SNR_linear);

        % Compute Sensitivity
        Sensitivity_dBm = -174 + 10 * log10(BW) + NF + SNR_dB;

        % Store in table
        Sensitivity_table(i, j) = Sensitivity_dBm;
    end
end

% Display the table
fprintf('\nSensitivity (dBm) at P_b = %.1e\n', Pb_target);
fprintf('         p = 0     p = 1     p = 2     p = 3     p = 4\n');
for i = 1:length(SF_values)
    fprintf('SF = %2d ', SF_values(i));
    fprintf('%9.2f', Sensitivity_table(i, :));
    fprintf('\n');
end



%% change of formula
EbN0_dB_range = -5:0.5:15;
EbN0 = 10.^(EbN0_dB_range / 10);
cr = 0;
SF_values = 7:12;



% Initialize figure
figure;

% Plot Pb for each SF with different markers
for idx = 1:length(SF_values)
    sf = SF_values(idx);
    term1 = (log(sf) / log(12)) / (2);
    term2 = 4 / (4 + cr);
    Pb = 0.5 * (erfc(2*term1 * term2 * EbN0));
    
    semilogy(EbN0_dB_range, Pb, 'DisplayName', ['SF = ', num2str(sf)], 'LineWidth', 1.5,'Marker', marker_list{idx},'MarkerIndices', 1:3:length(EbN0_dB_range));  % Show markers at intervals
    hold on;
end

title('LoRa CSS P_b as a function of E_b/N_0');
xlabel('E_b/N_0 (dB)');
ylabel('P_b (Average Bit Error Probability)');
legend show;
grid on;
ylim([1e-7, 1]);
hold

%% Main script Pb function EbN0
% EbN0_dB_range = -5:1:15;
% alpha =  10.^(EbN0_dB_range / 10);
% 
% SF_values = 7:12;
% 
% % Initialize figure
% %figure;
% 
% % Plot Pb for each SF
% for sf = SF_values
%     Pb = 0.5 * erfc((1.28 * sqrt(sf * alpha) - 1.28 * sqrt(sf) + 0.4) / sqrt(2));
%     semilogy(EbN0_dB_range, Pb, 'DisplayName', ['SF = ', num2str(sf)]);
%     hold on;
% end
% title('LoRa CSS P_b as a function of E_b/N_0');
% xlabel('E_b/N_0 (dB)');
% ylabel('P_b (Average Bit Error Probability)');
% legend show;
% grid on;
% ylim([1e-7, 1]);
% hold off;
