%% LoRa sensitivity study
% AUTHOR: Clara SORRE
% This MATLAB code analyzes the bit error probability (P_b), the signal-to-noise ratio (SNR),
% and the receiver sensitivity for LoRa using the Chirp Spread Spectrum (CSS) modulation defined in the report.
% It shows the impact of spreading factor (SF) and coding rate (parity bits) on the performance.

clc; clear; close all;

%% 1- assumptions
BW=125e3;
NF= 6;
marker_list = {'o', 's', 'd', '^', 'v', '>'}; % the markers for each SF


%% 2- Q function plotted
x = -5:0.01:5; % define x for normal gaussian distribution
pdf = (1/sqrt(2*pi)) * exp(-x.^2 / 2);
z_vals = -5:0.01:5;% z values for Q-function
Q = @(z) 0.5 * erfc(z / sqrt(2));
q_vals = Q(z_vals);
%z-value to illustrate Q(z) area
z_example = 0;
figure;
% plot 1: Normal distribution with Q(z) tail
subplot(1,2,1);
plot(x, pdf, 'b', 'LineWidth', 2); hold on;
area_x = z_example:0.01:5;
area_y = (1/sqrt(2*pi)) * exp(-area_x.^2 / 2);
area(area_x, area_y, 'FaceColor', [1 0.6 0.6], 'EdgeColor', 'none');
line([z_example z_example], [0 (1/sqrt(2*pi)) * exp(-z_example^2 / 2)], 'Color', 'r', 'LineStyle', '--');
text(z_example + 0.2, 0.1, sprintf('Q(%.1f)', z_example), 'Color', 'r');
xlabel('x');
ylabel('Probability density function');
title('Standard Normal Distribution');
legend('Probability density function', 'Q(z) area', 'Location', 'northeast');
grid on;

% Plot 2: Q-function
subplot(1,2,2);
plot(z_vals, q_vals, 'r', 'LineWidth', 2);
hold on;
plot(z_example, Q(z_example), 'ko', 'MarkerFaceColor', 'k');
text(z_example + 0.2, Q(z_example), sprintf('Q(%.1f) = %.3f', z_example, Q(z_example)));
xlabel('z');
ylabel('Q(z)');
title('Q-function');
grid on;

%% 3- LoRa CSS P_b as a function of E_b/N_0
EbN0_dB_range = -5:0.5:15;
EbN0 = 10.^(EbN0_dB_range / 10);
cr = 0;
SF_values = 7:12;
figure;
for idx = 1:length(SF_values) %plot pb for each sf
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



%% 4- SNR as a function of z for different SF values
z_values = linspace(0, 111, 1000); % start from 0

p = 0; % parity bits
SF_values = 7:12; % spreading factor

figure; 
hold on;
colors = lines(length(SF_values));

for idx = 1:length(SF_values)
    SF = SF_values(idx);
    SNR_linear = ((z_values * sqrt(2) * SF) / (log10(SF)/log10(12))).* (4 / (4 + p)) ./ (2.^SF); % compute the SNR in linear scale
    SNR_dB = 10 * log10(SNR_linear); % convert to dB
    plot(z_values, SNR_dB, 'Color', colors(idx, :), 'LineWidth', 1.5, 'DisplayName', sprintf('SF = %d', SF), 'Marker', marker_list{idx}, ...
        'MarkerIndices', 1:50:length(z_values)); %plot SNR vs. z
end

xlabel('z');
ylabel('SNR (dB)');
title(['SNR as a function of z for different SF values for p=', num2str(p)]);
grid on;
legend show;
hold off;


%% 5- SNR as a function of z for different p values

z_values = linspace(0, 111, 1000); % start from 0
SF = 12; 
p_values = 0:4; % parity bits
figure; 
hold on;
colors = lines(length(p_values));

for idx = 1:length(p_values)
    p = p_values(idx);
    SNR_linear = ((z_values * sqrt(2) * SF) / (log10(SF)/log10(12))).* (4 / (4 + p)) ./ (2.^SF); % compute SNR in linear scale
    SNR_dB = 10 * log10(SNR_linear); % convert to dB
    plot(z_values, SNR_dB, 'Color', colors(idx, :), 'LineWidth', 1.5, 'DisplayName', sprintf('p = %d', p), 'Marker', marker_list{idx},'MarkerIndices', 1:50:length(z_values));
end

xlabel('z');
ylabel('SNR (dB)');
title(['SNR as a function of z for different p values for SF=', num2str(SF)]);
grid on;
legend show;
hold off;


%% 6- z function of Pb

Pb_values = linspace(10e-100, 10e-6 , 1000000); %Pb range
z_values = sqrt(2) * erfcinv(2 * Pb_values); % compute z from Pb
figure; %plot
plot(Pb_values, z_values, 'b', 'LineWidth', 1.5);
xlabel('P_b');
ylabel('z');
title('z as a function of P_b');
grid on;


%% 7- SNR_{dB} as a function of P_{b} for different SF values

Pb_values = linspace(0, 10e-3, 1000); 
z_values = sqrt(2) * erfcinv(2* Pb_values); % compute z from Pb
p = 0; % parity bits
SF_values = 7:12; % SF
figure; 
hold on;
colors = lines(length(SF_values));

for idx = 1:length(SF_values)
    SF = SF_values(idx);
    SNR_linear = ((z_values * sqrt(2) * SF) / (log10(SF)/log10(12))).* (4 / (4 + p)) ./ (2.^SF); % compute SNR in linear scale
    SNR_dB = 10 * log10(SNR_linear); %convert to dB
    plot(Pb_values, SNR_dB, 'Color', colors(idx, :), 'LineWidth', 1.5, 'DisplayName', sprintf('SF = %d', SF),'Marker', marker_list{idx},'MarkerIndices', 1:50:length(z_values));
end

xlabel('P_{b}');
ylabel('SNR (dB)');
title(['SNR_{dB} as a function of P_{b} for different SF values for p=', num2str(p)]);
grid on;
legend show;
hold off;

%% 8- SNR_{dB} as a function of P_{b} for different p values 
Pb_values = linspace(0, 10e-3, 1000); 
z_values = sqrt(2) * erfcinv(2* Pb_values);%compute z from Pb
SF = 12; % parity bits
p_values = 0:4; % SF
figure; 
hold on;
colors = lines(length(p_values));
for idx = 1:length(p_values)
    p = p_values(idx);
    SNR_linear = ((z_values * sqrt(2) * SF) / (log10(SF)/log10(12))).* (4 / (4 + p)) ./ (2.^SF);% compute SNR in linear scale
    SNR_dB = 10 * log10(SNR_linear);
    plot(Pb_values, SNR_dB, 'Color', colors(idx, :), 'LineWidth', 1.5, 'DisplayName', sprintf('p = %d', p), 'Marker', marker_list{idx},'MarkerIndices', 1:50:length(z_values));
end

xlabel('P_{b}');
ylabel('SNR (dB)');
title(['SNR_{dB} as a function of P_{b} for different p values for SF=', num2str(SF)]);
grid on;
legend show;
hold off;


%% 9- Sensitivity_{dB} as a function of P_{b} for different SF values
Pb_values = linspace(10e-10, 10e-4  , 1000);  
z_values = sqrt(2) * erfcinv(2 * Pb_values);% compute z from Pb
p = 0; % parity bits
SF_values = 7:12;
figure; 
hold on;
colors = lines(length(SF_values));
for idx = 1:length(SF_values)
    SF = SF_values(idx);
    SNR_linear = ((z_values * sqrt(2) * SF) / (log10(SF)/log10(12))).* (4 / (4 + p)) ./ (2.^SF); % compute SNR in linear scale
    % convert to dB
    SNR_dB = 10 * log10(SNR_linear);
    Sensitivity_dBm = -174 + 10 * log10(BW) + NF + SNR_dB;
    plot(Pb_values, Sensitivity_dBm, 'Color', colors(idx, :), 'LineWidth', 1.5, 'DisplayName', sprintf('SF = %d', SF), 'Marker', marker_list{idx},'MarkerIndices', 1:50:length(z_values));
end

xlabel('P_{b}');
ylabel('Sensitivity (dBm)');
title(['Sensitivity_{dB} as a function of P_{b} for different SF values for p=', num2str(p)]);
grid on;
legend show;
hold off;

%% 10- Sensitivity_{dB} as a function of P_{b} for p that varies
Pb_values = linspace(10e-10, 10e-4, 1000); % probability of bit error
z_values = sqrt(2) * erfcinv(2 * Pb_values);
SF = 12; % fixed Spreading Factor
p_values = 0:4; % parity bits
figure;
hold on;
colors = lines(length(p_values));
for idx = 1:length(p_values)
    p = p_values(idx);
    SNR_linear = ((z_values * sqrt(2) * SF) / (log10(SF)/log10(12))).* (4 / (4 + p)) ./ (2.^SF); %compute SNR in linear scale
    SNR_dB = 10 * log10(SNR_linear);
    Sensitivity_dBm = -174 + 10 * log10(BW) + NF + SNR_dB;
    plot(Pb_values, Sensitivity_dBm, 'Color', colors(idx, :), 'LineWidth', 1.5, 'DisplayName', sprintf('p = %d', p), 'Marker', marker_list{idx},'MarkerIndices', 1:50:length(z_values));
end

xlabel('P_{b}');
ylabel('Sensitivity (dBm)');
title(['Sensitivity_{dB} as a function of P_{b} for SF = ', num2str(SF)]);
grid on;
legend show;
hold off;


%% 11- Final Sensitivity Table: sensitivity (dBm) for each (SF, p) at fixed Pb
BW=125e3;
NF= 6;
Pb_target = 1e-3; % target probability of bit error
z = sqrt(2) * erfcinv(2 * Pb_target); % z value

SF_values = 7:12;
p_values = 0:4;

Sensitivity_table = zeros(length(SF_values), length(p_values)); % initialize result of the matrix

for i = 1:length(SF_values)
    SF = SF_values(i);
    for j = 1:length(p_values)
        p = p_values(j);
        SNR_linear = ((z * sqrt(2) * SF) / (log10(SF)/log10(12))).* (4 / (4 + p)) ./ (2.^SF);% compute SNR in linear scale
        SNR_dB = 10 * log10(SNR_linear);
        Sensitivity_dBm = -174 + 10 * log10(BW) + NF + SNR_dB;
        Sensitivity_table(i, j) = Sensitivity_dBm;  % store in table
    end
end

fprintf('\nSensitivity (dBm) at P_b = %.1e\n', Pb_target);
fprintf('         p = 0     p = 1     p = 2     p = 3     p = 4\n');
for i = 1:length(SF_values)
    fprintf('SF = %2d ', SF_values(i));
    fprintf('%9.2f', Sensitivity_table(i, :));
    fprintf('\n');
end




