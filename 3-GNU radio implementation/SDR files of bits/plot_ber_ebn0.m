%% GNU Radio BER vs Eb/N0
% AUTHOR: Clara SORRE
% This MATLAB code ...........

eb_n0_dB = [-5, 0, 5, 10, 15];  % Eb/N0 values in dB 
target_ber = 1e-3;
interpolate_target_ber = @(ebn0, ber, target) interp1(log10(ber), ebn0, log10(target), 'linear', 'extrap'); %interpolation function
% BER results measured gnu radio
ber_lorasf7 = [0.396657, 0.229680, 0.005685, 0.00002, 0.00000001];   %when ber=0: artificially 0.00000001
ber_bfsk = [0.42288, 0.34419, 0.15655, 0.01876, 0.00004];  
ber_4fsk = [0.39636, 0.29288, 0.1293, 0.01482, 0.00004];   
ber_8fsk = [0.315653, 0.201418, 0.066603, 0.00315, 0.0000001]; %when ber=0: artificially 0.0000001
ber_16fsk = [0.24627, 0.1374, 0.03234, 0.00069, 0.0000001]; %when ber=0: artificially 0.0000001
ber_qpsk = [0.372220, 0.159881, 0.059364, 0.000780, 0.000010]; 
ber_lorasf12 = [0.213394, 0.167567,0.000310,  0.00000002, 0.00000001];
schemes = {
    'LoRa SF=7', ber_lorasf7;
    'LoRa SF=12', ber_lorasf12;
    'BFSK', ber_bfsk;
    '4-FSK', ber_4fsk;
    '8-FSK', ber_8fsk;
    '16-FSK', ber_16fsk;
    'QPSK', ber_qpsk
};
fprintf('Estimated E_b/N_0 for BER = %.1e:\n', target_ber);
for i = 1:size(schemes,1)
    name = schemes{i,1};
    ber = schemes{i,2};
    est = interpolate_target_ber(eb_n0_dB, ber, target_ber);
    fprintf('%s: %.2f dB\n', name, est);
end
% plot
figure;
semilogy(eb_n0_dB, ber_lorasf7, 'o-', 'LineWidth', 1.5,'DisplayName', sprintf('LoRa SF=7'));
hold on;
semilogy(eb_n0_dB, ber_lorasf12, '<-', 'LineWidth', 1.5,'DisplayName', sprintf('LoRa SF=12'));
semilogy(eb_n0_dB, ber_bfsk, 's-', 'LineWidth', 1.5,'DisplayName', sprintf('BFSK'));
semilogy(eb_n0_dB, ber_4fsk, 'd-', 'LineWidth', 1.5,'DisplayName', sprintf('4-FSK'));
semilogy(eb_n0_dB, ber_8fsk, '^-', 'LineWidth', 1.5,'DisplayName', sprintf('8-FSK'));
semilogy(eb_n0_dB, ber_16fsk, 'v-', 'LineWidth', 1.5,'DisplayName', sprintf('16-FSK'));
semilogy(eb_n0_dB, ber_qpsk, '>-', 'LineWidth', 1.5,'DisplayName', sprintf('QPSK'));
grid on;
legend show;
xlabel('E_b/N_0 (dB)');
ylabel('Bit Error Rate (BER)');
title('BER vs. E_b/N_0 for modulation schemes');
xlim([min(eb_n0_dB)-1, max(eb_n0_dB)+1]);
ylim([1e-4 1]);

%display the data on the graph
for i = 1:length(eb_n0_dB)
    text(eb_n0_dB(i), ber_lorasf7(i)*1.2, sprintf('%.1e', ber_lorasf7(i)),'HorizontalAlignment', 'center', 'FontSize', 8);
    text(eb_n0_dB(i), ber_lorasf12(i)*1.2, sprintf('%.1e', ber_lorasf12(i)),'HorizontalAlignment', 'center', 'FontSize', 8);
    text(eb_n0_dB(i), ber_bfsk(i)*1.2, sprintf('%.1e', ber_bfsk(i)),'HorizontalAlignment', 'center', 'FontSize', 8);
    text(eb_n0_dB(i), ber_4fsk(i)*1.2, sprintf('%.1e', ber_4fsk(i)),'HorizontalAlignment', 'center', 'FontSize', 8);
    text(eb_n0_dB(i), ber_8fsk(i)*1.2, sprintf('%.1e', ber_8fsk(i)),'HorizontalAlignment', 'center', 'FontSize', 8);
    text(eb_n0_dB(i), ber_16fsk(i)*1.2, sprintf('%.1e', ber_16fsk(i)),'HorizontalAlignment', 'center', 'FontSize', 8);
    text(eb_n0_dB(i), ber_qpsk(i)*1.2, sprintf('%.1e', ber_qpsk(i)),'HorizontalAlignment', 'center', 'FontSize', 8);
end


