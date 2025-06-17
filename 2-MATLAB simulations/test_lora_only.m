%% TEST LoRa modulation and demodulation of a message of bits
% AUTHOR: Clara SORRE
% DESCRIPTION OF THE CODE: This MATLAB code simulates the BER of LoRa CSS 
% modulation under AWGN for SF from 7 to 12.
% The simulation focuses on the modulation and demodulation only, the
% preprocessing such as Whitening, Gray mapping, Hamming encoding and also
% interleaving are not considered here.
% IN INPUT: msg_bits_length = the total number of random input bits to simulate
% IN OUTPUT: the BER vs Eb/N0 plot (both theoretical and simulated graphs)

function test_lora_only(msg_bits_length) % LoRa BER Simulation
    marker_list = {'o', 's', 'd', '^', 'v', '>'}; % the markers for each SF
    tic; % start the timing
    n = msg_bits_length; % nb of bits to simulate
    b = randi([0, 1], 1, n); % random bit sequence
    bw = 125e3;  % LoRa bandwidth in Hz
    EbN0_dB = -5:0.5:15; 
    SF_range = 7:12;
    colors = lines(length(SF_range));
    ber_matrix = zeros(length(SF_range), length(EbN0_dB)); % rows = SF, columns = Eb/N0
    for sf_idx = 1:length(SF_range)
        SF = SF_range(sf_idx);
        M = 2^SF; %number of possible symbols
        Ns = M; % samples per symbol
        Rs = bw / M; % symbol rate
        %Rb = Rs * SF; % Bit rate
        num_symbols = floor(n / SF);  % full LoRa symbols
        b_cut = b(1:num_symbols * SF); % trim the bits
        b_matrix = reshape(b_cut, SF, []).'; % SF bits per row
        data_tx = bi2de(b_matrix, 'left-msb'); % convert to integers
        t = (0:Ns-1) / bw; % Time vector
        Ts = 1 / Rs;
        downchirp = exp(-1j * 2 * pi * ((-bw/2)*t + (bw/(2*Ts))*t.^2)); %downchirp
        base_chirp = exp(1j * 2 * pi * (bw/(2*Ts)) * t.^2); %base up chirp 
        f0_all = -bw/2 + (bw/M)*(0:M-1); % frequency offsets
        for eb_idx = 1:length(EbN0_dB)
            EbN0 = EbN0_dB(eb_idx);
            chirps = base_chirp .* exp(1j * 2 * pi * f0_all.' * t); % all the reference chirps (M x Ns)
            tx_lora_matrix = chirps(data_tx + 1, :); %chirps for each symbol
            tx_lora = reshape(tx_lora_matrix.', 1, []); % flatten to 1D waveform
            ebn0_lin = 10^(EbN0 / 10);  % convert Eb/N0 to linear scale 
            rb = bw / (2^SF) * SF;  % estimated bit rate
            snr_lin = ebn0_lin *  (rb / bw);  % convert Eb/N0 to SNR per symbol
            snr_dB = 10 * log10(snr_lin);
            tx_lora_noisy = awgn(tx_lora, snr_dB, 'measured'); %function awgn matlab
            rx_symbols_matrix = reshape(tx_lora_noisy, Ns, []); % receiver processing: reshape
            dechirped = rx_symbols_matrix .* downchirp.';% dechirp
            fft_out = abs(fft(dechirped));%fft
            [~, k_hat] = max(fft_out, [], 1);  % symbol detection by FFT peak
            data_rx = mod(k_hat - 1, M); %due to indexing
            rx_bits_matrix = de2bi(data_rx, SF, 'left-msb');%back to bits
            rx_bits = reshape(rx_bits_matrix.', 1, []);
            bit_errors = sum(rx_bits ~= b_cut); % BER
            ber_matrix(sf_idx, eb_idx) = bit_errors / length(b_cut);
        end
    end
    grid on; % plot Simulated BER
    for sf_idx = 1:length(SF_range)
        semilogy(EbN0_dB, ber_matrix(sf_idx, :), 'LineStyle','--', 'LineWidth', 1.5,'DisplayName', ['Sim SF = ' num2str(SF_range(sf_idx))], 'Marker', marker_list{sf_idx}, 'Color', colors(sf_idx, :), 'MarkerIndices', 1:5:length(EbN0_dB));
        hold on;
    end
    cr=0;
    EbN0 = 10.^(EbN0_dB / 10);% plot Theoretical P_b
    sf_idx2=1;
    for sf = SF_range
        term1= (log10(sf)/log10(12)) / (2);
        term2 = 4 / (4 + cr);
        Pb = 0.5 * (erfc(2*term1 * term2 * EbN0));
        semilogy(EbN0_dB, Pb,'LineWidth', 1.5, 'DisplayName', ['Theory SF = ' num2str(sf)], 'Marker', marker_list{sf_idx2},'Color', colors(sf_idx2, :),  'MarkerIndices', 1:5:length(EbN0_dB));
        sf_idx2= sf_idx2 +1;
        hold on;
    end
    xlabel('E_b/N_0 [dB]');
    ylabel('Bit Error Rate (BER)');
    title('LoRa Simulated and Theoretical BER vs E_b/N_0 for different SF');
    legend('show');
    grid on;
    ylim([1e-5, 1]);
    hold off;
    toc; % end timing
end
