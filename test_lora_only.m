function test_lora_only(msg_bits_length) % LoRa BER Simulation
    %colors = {'r', 'g', 'b', 'c', 'm', 'y'};  % predefined color name% plot colors
    % colors = [1.0, 0.0, 0.0;   % red
    % 0.0, 1.0, 0.0;   % green
    % 0.0, 0.0, 1.0;   % blue
    % 0.0, 1.0, 1.0;   % cyan
    % 1.0, 0.0, 1.0;   % magenta
    % 1.0, 1.0, 0.0;   % yellow
    % 0.0, 0.0, 0.0];  % black

    tic; % start the timing
    n = msg_bits_length; % nb of bits to simulate
    b = randi([0, 1], 1, n); % random bit sequence
    bw = 125e3;  % LoRa bandwidth in Hz
    fs = bw;     % sampling frequency
    EbN0_dB = -5:1:15; 
    SF_range = 7:12;
    ber_matrix = zeros(length(SF_range), length(EbN0_dB)); % Rows = SF, Cols = Eb/N0
    for sf_idx = 1:length(SF_range)
        SF = SF_range(sf_idx);
        M = 2^SF;
        Ns = M; % samples per symbol
        Rs = bw / M; % symbol rate
        %Rb = Rs * SF; % Bit rate
        num_symbols = floor(n / SF);  % full LoRa symbols
        b_cut = b(1:num_symbols * SF); % trim bits
        b_matrix = reshape(b_cut, SF, []).'; % SF bits per row
        data_tx = bi2de(b_matrix, 'left-msb'); % convert to integers

        t = (0:Ns-1) / fs; % Time vector
        Ts = 1 / Rs;
        downchirp = exp(-1j * 2 * pi * ((-bw/2)*t + (bw/(2*Ts))*t.^2));
        base_chirp = exp(1j * 2 * pi * (bw/(2*Ts)) * t.^2); 
        f0_all = -bw/2 + (bw/M)*(0:M-1); % frequency offsets

        for eb_idx = 1:length(EbN0_dB)
            EbN0 = EbN0_dB(eb_idx);
            chirps = base_chirp .* exp(1j * 2 * pi * f0_all.' * t); % all reference chirps (M x Ns)
            tx_lora_matrix = chirps(data_tx + 1, :); 
            tx_lora = reshape(tx_lora_matrix.', 1, []); % flatten to 1D

           
            % %Eb = sum(abs(tx_lora).^2) / (num_symbols * SF); 
            % Es = sum(abs(tx_lora).^2) / length(data_tx);   % Energy per symbol
            % %Rb = Rs * SF;                                  % Bit rate
            % Eb = Es / SF;                                  % Energy per bit (since SF bits per symbol)
            % 
            % 
            % N0 = Eb / (10^(EbN0 / 10));
            % noise = sqrt(N0 / 2) * (randn(1, length(tx_lora)) + 1j * randn(1, length(tx_lora)));
            % tx_lora_noisy = tx_lora + noise;
            ebn0_lin = 10^(EbN0 / 10);  % Convert Eb/N0 to linear scale
            
            bits_per_sym = SF;  % For LoRa, each symbol carries SF bits
            rb = bw / (2^SF) * SF;  % Estimated bit rate
            %rs = bw_values(idx) / (2^SF);       % Symbol rate
            snr_lin = ebn0_lin *  (rb / bw);  % Convert Eb/N0 to SNR per symbol bits_per_sym *
            snr_dB = 10 * log10(snr_lin);
            
            tx_lora_noisy = awgn(tx_lora, snr_dB, 'measured');


            rx_symbols_matrix = reshape(tx_lora_noisy, Ns, []);
            dechirped = rx_symbols_matrix .* downchirp.';
            fft_out = abs(fft(dechirped));
            [~, k_hat] = max(fft_out, [], 1);
            data_rx = mod(k_hat - 1, M);
            rx_bits_matrix = de2bi(data_rx, SF, 'left-msb');
            rx_bits = reshape(rx_bits_matrix.', 1, []);
            bit_errors = sum(rx_bits ~= b_cut); % BER
            ber_matrix(sf_idx, eb_idx) = bit_errors / length(b_cut);
        end
    end
    color_idx1=1;
    grid on; % plot Simulated BER
    for sf_idx = 1:length(SF_range)
        semilogy(EbN0_dB, ber_matrix(sf_idx, :),  '--','DisplayName', ['Sim SF = ' num2str(SF_range(sf_idx))]); %,'Color', colors(color_idx1)
        color_idx1=color_idx1+1;
        hold on;
    end
    cr=0;
    color_idx=1;
    EbN0 = 10.^(EbN0_dB / 10);% plot Theoretical P_b
    for sf = SF_range
        term1= (log10(sf)/log10(12)) / (2);
        term2 = 4 / (4 + cr);
        %argument = term1 * term2 * EbN0;
        Pb = 0.5 * (erfc(2*term1 * term2 * EbN0));
        %Pb = 0.5 * erfc((1.28 * sqrt(sf * alpha) - 1.28 * sqrt(sf) + 0.4) / sqrt(2));
        semilogy(EbN0_dB, Pb, 'DisplayName', ['Theory SF = ' num2str(sf)]); %,'Color', colors(color_idx)
        color_idx=color_idx+1;
        hold on;
    end
    xlabel('E_b/N_0 [dB]');
    ylabel('Bit Error Rate (BER)');
    title('LoRa Simulated and Theoretical BER vs E_b/N_0 for different SF');
    legend('show');
    grid on;
    ylim([1e-7, 1]);

    hold off;
    toc; % end timing
end
