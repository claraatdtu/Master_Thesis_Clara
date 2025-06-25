%% TEST MFSK modulation and demodulation of a message of bits
% AUTHOR: Clara SORRE
% DESCRIPTION OF THE CODE: This MATLAB code simulates the BER of MFSK
% modulation schemes for M = 2,4,8 and 16, under AWGN.
% The simulation focuses on the modulation, the transmission and demodulation. 
% The center frequency used here is 1 (simulated directly in the baseband, without simulating a high-frequency RF carrier) 
% with the assumption that they would later be reconverted to the carrier frequency used in a the real system.
% IN INPUT: msg_bits_length = the total number of random input bits to simulate
% IN OUTPUT: the BER vs Eb/N0 plot (both theoretical and simulated graphs)

function test_mfsk_only(msg_bits_length)
    tic; %elapsed time 
    n = msg_bits_length;% number of bits
    b = randi([0, 1], 1, n); % random binary message
    samples_per_bit_map = containers.Map( {2, 4, 8, 16},{28, 28, 37, 64} ); %add calibrated values of samples per bit for trade off time/orthogonality
    %a map is used: associate modulations M to samples per bit
    bitrate = 366.2; % bitrate in bits per second HERE FOR SF=7
    ebn0_range = -5:1:15;% Eb/N0 range in dB
    BER_sim = struct(); % simulated BER
    BER_th = struct(); % theoretical BER
    grid on;
    colors = ['b', 'r', 'g', 'm'];% plot colors
    marker_list = {'o', 's', 'd', '^'};
    M_values = [2, 4, 8, 16];
    bw_values = zeros(1, length(M_values));  % preallocate bw
    %loop through each M order
    for idx = 1:length(M_values) 
        M = M_values(idx);
        bits_per_sym = log2(M); %number of bits per FSK symbol
        bw_values(idx) = (M * bitrate) / bits_per_sym; %bw associated 
        nb = floor(n / bits_per_sym) * bits_per_sym;
        bM = b(1:nb);  % bits to a multiple of bits_per_sym: bit stream evenly divided into symbols
        current_samples_per_bit = samples_per_bit_map(M);
        t_current = linspace(0, 1, current_samples_per_bit + 1); %time vector for 1 symbol
        t_current(end) = [];  % remove the last point to have exactly current_samples_per_bit samples
        symM = zeros(current_samples_per_bit, M); %M orthogonal waveforms
        for k = 0:M-1
            fk = (1+k);  % frequencies: 1 to M: simulate in base band and then upconverted for transmission
            wave = sin(2*pi*fk*t_current);%sinus wave
            symM(:, k+1) = wave / norm(wave);  % normalized waveform
        end
        bM_matrix = reshape(bM, bits_per_sym, []).'; %reshape the bits into symbols
        idxM = bi2de(bM_matrix, 'left-msb') + 1; %map the bits to symbol indices
        tx_signal = symM(:, idxM);  % each column is a symbol
        tx_signal = tx_signal(:).';  % serialize to 1D
        BER = zeros(size(ebn0_range));
        for i = 1:length(ebn0_range)
            ebn0_dB = ebn0_range(i);
            ebn0_lin = 10^(ebn0_dB / 10); %convert to linear
            snr_sym = ebn0_lin * bits_per_sym * (bitrate / bw_values(idx)); %obtain snr dB instead of EbN0
            snr_sym_dB = 10*log10(snr_sym); 
            rx_signal = awgn(tx_signal, snr_sym_dB, 'measured'); %add the awgn
            rx_matrix = reshape(rx_signal, current_samples_per_bit, []); %reshape into symbols
            energies = zeros(M, size(rx_matrix, 2)); %energy detection
            for m = 1:M
                template = symM(:, m);
                energies(m, :) = sum((rx_matrix .* template).^2, 1); %correlation
            end
            [~, id_detected] = max(energies, [], 1);%pick the symbol with the highest energy
            rx_bits = de2bi(id_detected - 1, bits_per_sym, 'left-msb').'; %demap symbols back to bits
            rx_bits = rx_bits(:).';%flatten to 1D
            BER(i) = mean(rx_bits ~= bM); % BER
        end
        BER_sim.(['F', num2str(M)]) = BER;
        semilogy(ebn0_range, BER, '--', 'LineWidth', 1.5,'DisplayName', [num2str(M) '-FSK Sim'], 'Marker', marker_list{idx},'Color', colors(idx)); % Plot the simulated BER
        hold on;
    end
    for idx = 1:length(M_values) % Theoretical BER (non-coherent orthogonal MFSK)
        M = M_values(idx);
        bits_per_sym = log2(M);
        Pb = zeros(size(ebn0_range));
        for k = 1:length(ebn0_range)
            EbN0 = 10^(ebn0_range(k) / 10);  % linear scale
            Pb_k = 0;
            for n = 1:M-1
                try
                    C = nchoosek(M - 1, n); %binomial coeff
                catch
                    C = 0; % overflow protection
                end
                if C == 0 || ~isfinite(C)
                    continue;
                end
                exponent = (-n * bits_per_sym * EbN0) / (n + 1);
                if exponent < -700  % prevent underflow in exp
                    exp_term = 0;
                else
                    exp_term = exp(exponent);
                end
                term = ((M / 2) / (M - 1)) * ((-1)^(n + 1) / (n + 1)) * C * exp_term;
                Pb_k = Pb_k + term;
            end
            Pb(k) = max(min(Pb_k, 1), 1e-10);  % avoid NaN/Inf
        end
        BER_th.(['F', num2str(M)]) = Pb;
        semilogy(ebn0_range, Pb, 'LineWidth', 1.5,'DisplayName', [num2str(M) '-FSK Theory'], 'Marker', marker_list{idx},'Color', colors(idx));%plot ber
    end
    xlabel('E_b/N_0 [dB]');     % plot
    ylabel('Bit Error Rate (BER)');
    legend('Location', 'southwest');
    title('Simulated vs. Theoretical BER for Non-coherent MFSK');
    ylim([1e-5, 1]);
    grid on;
    hold off;
    toc;
end









        % try
        %     tx_signal = symM(:, idxM);  % each column is a symbol
        % catch
        %     error('Indexing into symM failed. Check dimensions: M=%d, max(idxM)=%d', M, max(idxM));
        % end

 % 
        % if any(idxM > M)
        %     error('Index exceeds number of MFSK symbols.');
        % end

        % Create symbol stream

        % symM = zeros(samples_per_bit, M);
        % if M==16
        %     for k = 0:M-1
        %         fk = 1 + k;  % frequencies: 1, 2, ..., M
        % 
        %         wave = sin(2*pi*fk*t2);
        %         symM(:, k+1) = wave / norm(wave);  % normalized waveforms
        %     end
        % else
        %     for k = 0:M-1
        %         fk = 1 + k;  % frequencies: 1, 2, ..., M
        % 
        %         wave = sin(2*pi*fk*t);
        %         symM(:, k+1) = wave / norm(wave);  % normalized waveforms
        %     end
        % end
        % Mapping bits to symbols
        % bM_matrix = reshape(bM, bits_per_sym, []).';
        % idxM = bi2de(bM_matrix, 'left-msb') + 1;
        % tx_signal = symM(:, idxM);  % symbols as columns
        % tx_signal = tx_signal(:).';  % serialize


 % samples_per_bit = 64;              % oversampling factor
    % samples_per_bit2 = 38;
    % samples_per_bit3 = 28;
    % samples_per_bit4 = 28;
    % t = 0:1/samples_per_bit:1 - 1/samples_per_bit;  % time vector for 1 symbol
    % t2 =  0:1/samples_per_bit2:1 - 1/samples_per_bit2;
    % t3 =  0:1/samples_per_bit3:1 - 1/samples_per_bit3;
    % t4 =  0:1/samples_per_bit4:1 - 1/samples_per_bit4;

    %cycles_per_symbol = 4; % Adjust this (2–5) to trade off symbol duration vs. orthogonality

    % Precompute samples_per_bit for each M
    % samples_per_bit_map = containers.Map('KeyType', 'double', 'ValueType', 'double');
    % for M = [2, 4, 8, 16]
    %     samples_per_bit_map(M) = M * cycles_per_symbol;
    % end
% Calibrated values for good orthogonality (empirically chosen) % M values % Corresponding samples_per_bit




        % if isempty(bM)
        %     error('Empty bit stream after trimming. Choose larger msg_bits_length.');
        % end

        % Generate orthogonal MFSK waveforms
            % Choose proper sample rate for current M
        % if M == 16
        %     current_samples_per_bit = samples_per_bit;
        %     t_current = t;
        % end
        % if M== 8
        %     current_samples_per_bit = samples_per_bit2;
        %     t_current = t2;
        % end
        % if M== 4
        %     current_samples_per_bit = samples_per_bit3;
        %     t_current = t3;
        % end
        % if M== 2
        %     current_samples_per_bit = samples_per_bit4;
        %     t_current = t4;
        % end



% function test_mfsk_only(msg_bits_length)
%     tic;
%     % Simulation parameters
%     n = msg_bits_length;                % number of bits
%     b = randi([0, 1], 1, n);            % random binary message
%     samples_per_bit = 40;              % oversampling factor
%     bitrate = 366.2;                   % bitrate in bits per second (e.g., LoRa-like)
%     %fs = bitrate * samples_per_bit;   % sampling frequency
%     t = 0:1/samples_per_bit:1 - 1/samples_per_bit;  % time vector for 1 symbol
% 
%     ebn0_range = -5:1:15;              % Eb/N0 range in dB
%     BER_sim = struct();                % simulated BER
%     BER_th = struct();                 % theoretical BER
%     %figure; 
%     grid on;  % <<== This ensures all plots are shown on the same figure
%     bw = [0 0 0 0];
%     for M = [2, 4, 8, 16]
%         bw(log2(M))=(M*bitrate)/log2(M)
%         bits_per_sym = log2(M);
%         nb = floor(n / bits_per_sym) * bits_per_sym;
%         bM = b(1:nb);  % trim bits to multiple of bits_per_sym
% 
%         symM = zeros(samples_per_bit, M);
%         for k = 0:M-1
%             fk = 1 + k;  % frequencies: 1, 2, ..., M
%             symM(:, k+1) = sin(2*pi*fk*t) / norm(sin(2*pi*fk*t));
%         end
% 
%         bM_matrix = reshape(bM, bits_per_sym, []).';
%         idxM = bi2de(bM_matrix, 'left-msb') + 1;
%         tx_signal = symM(:, idxM);
%         tx_signal = tx_signal(:).';  % flatten
% 
%         % Simulated BER
%         BER = zeros(size(ebn0_range));
%         for i = 1:length(ebn0_range)
%             ebn0 = ebn0_range(i);
%             snr = (10^(ebn0/10))*(bitrate/bw(log2(M)));
%             snr_adj= 10*log10(snr)+ 10*log10(bits_per_sym);% 10*log10(samples_per_bit);
%             %ebn0_adjusted = ebn0 + 10*log10(bits_per_sym) ;%- 10*log10(samples_per_bit);
% 
%             rx_signal = awgn(tx_signal, snr_adj, 'measured');
% 
%             rx_matrix = reshape(rx_signal, samples_per_bit, []);
%             corrs = symM.' * rx_matrix;
%             [~, id_detected] = max(corrs, [], 1);
%             rx_bits = de2bi(id_detected - 1, bits_per_sym, 'left-msb').';
%             rx_bits = rx_bits(:).';
%             BER(i) = mean(rx_bits ~= bM);
%         end
%         BER_sim.(['F', num2str(M)]) = BER;
% 
%         % Theoretical BER (approximation for orthogonal MFSK, coherent detection)
%         EbN0_lin = 10.^(ebn0_range/10);
%         Pb = (M-1)/(2*log2(M)) * erfc(sqrt(log2(M)*EbN0_lin/2));
%         BER_th.(['F', num2str(M)]) = Pb;
% 
%         % Plot
%         semilogy(ebn0_range, BER, 'DisplayName', [num2str(M) '-FSK Sim']);
% 
%         hold on; 
% 
%         semilogy(ebn0_range, Pb, '--', 'DisplayName', [num2str(M) '-FSK Theory']);
%         %hold on;
%     end
% 
%     xlabel('E_b/N_0 [dB]');
%     ylabel('Bit Error Rate (BER)');
%     legend('Location', 'southwest');
%     title('Simulated vs. Theoretical BER for MFSK');
%     grid on;
% 
%     ylim([1e-7, 1]);
%     hold off;
%     toc;
% end


% 
% function test_mfsk_only(msg_bits_length)
%     tic;
%     % Simulation parameters
%     n = msg_bits_length;                % number of bits
%     b = randi([0, 1], 1, n);            % random binary message
%     samples_per_bit = 40;              % oversampling factor
%     bitrate = 366.2;                   % bitrate in bits per second (e.g., LoRa-like)
%     t = 0:1/samples_per_bit:1 - 1/samples_per_bit;  % time vector for 1 symbol
% 
%     ebn0_range = -5:1:15;              % Eb/N0 range in dB
%     BER_sim = struct();                % simulated BER
%     BER_th = struct();                 % theoretical BER
% 
%     grid on;
%     bw = [0 0 0 0];
%     for M = [2, 4, 8, 16]
%         bw(log2(M)) = (M * bitrate) / log2(M);
%         bits_per_sym = log2(M);
%         nb = floor(n / bits_per_sym) * bits_per_sym;
%         bM = b(1:nb);  % trim bits to multiple of bits_per_sym
% 
%         % Generate orthogonal MFSK waveforms
%         symM = zeros(samples_per_bit, M);
%         for k = 0:M-1
%             fk = 1 + k;  % frequencies: 1, 2, ..., M
%             wave = sin(2*pi*fk*t);
%             symM(:, k+1) = wave / norm(wave);  % normalize energy
%         end
% 
%         % Mapping bits to symbols
%         bM_matrix = reshape(bM, bits_per_sym, []).';
%         idxM = bi2de(bM_matrix, 'left-msb') + 1;
%         tx_signal = symM(:, idxM);  % each column is a symbol
%         tx_signal = tx_signal(:).';  % serialize
%         disp(bw(log2(M)))
%         % Simulated BER
%         BER = zeros(size(ebn0_range));
%         for i = 1:length(ebn0_range)
%             ebn0 = ebn0_range(i);
%             snr = (10^(ebn0/10))*(bitrate/bw(log2(M)));
%             snr_adj= 10*log10(snr)+ 10*log10(bits_per_sym);
% %             snr_adj= 10*log10(snr)+ 10*log10(bits_per_sym);% 10*log10(samples_per_bit);
%             %snr = ebn0 + 10*log10(bits_per_sym);  % SNR per symbol (non-coherent)
%             rx_signal = awgn(tx_signal, snr_adj, 'measured');
% 
%             % Reshape into symbol-long chunks
%             rx_matrix = reshape(rx_signal, samples_per_bit, []);
% 
%             % Non-coherent detection: energy-based
%             energies = zeros(M, size(rx_matrix, 2));
%             for m = 1:M
%                 template = symM(:, m);  % reference waveform
%                 energies(m, :) = sum((rx_matrix .* template).^2, 1);  % energy detection
%             end
%             [~, id_detected] = max(energies, [], 1);
% 
%             % Demap symbols back to bits
%             rx_bits = de2bi(id_detected - 1, bits_per_sym, 'left-msb').';
%             rx_bits = rx_bits(:).';
% 
%             % Compute BER
%             BER(i) = mean(rx_bits ~= bM);
%         end
%         BER_sim.(['F', num2str(M)]) = BER;
% 
%         % Theoretical BER (your original formula, assuming non coherent detection)
%         EbN0_lin = 10.^(ebn0_range/10);
%         Pb = (M - 1)/(2 * log2(M)) * erfc(sqrt(log2(M) * EbN0_lin / 2));
%         BER_th.(['F', num2str(M)]) = Pb;
% 
%         % Plot
%         semilogy(ebn0_range, BER, 'DisplayName', [num2str(M) '-FSK Sim']);
%         hold on; 
%         semilogy(ebn0_range, Pb, '--', 'DisplayName', [num2str(M) '-FSK Theory']);
%     end
% 
%     xlabel('E_b/N_0 [dB]');
%     ylabel('Bit Error Rate (BER)');
%     legend('Location', 'southwest');
%     title('Simulated vs. Theoretical BER for MFSK (Non-coherent detection)');
%     grid on;
%     ylim([1e-7, 1]);
%     hold off;
%     toc;
% end
% 
% 
% function test_mfsk_only(msg_bits_length)
%     tic;
% 
%     % Simulation parameters
%     n = msg_bits_length;                % number of bits
%     b = randi([0, 1], 1, n);            % random binary message
%     samples_per_bit = 40;              % oversampling factor
%     bitrate = 366.2;                   % bitrate in bits per second
%     t = 0:1/samples_per_bit:1 - 1/samples_per_bit;  % time vector for 1 symbol
%     ebn0_range = -5:1:15;              % Eb/N0 range in dB
% 
%     BER_sim = struct();                % simulated BER results
%     BER_th = struct();                 % theoretical BER results
% 
%     %figure;
%     grid on;
% 
% 
%     colors = ['b', 'r', 'g', 'm'];     % plot colors
% 
%     M_values = [2, 4, 8, 16];
%     bw_values = zeros(1, length(M_values));  % preallocate BW array
% 
%     for mi = 1:length(M_values)
%         M = M_values(mi);
%         bits_per_sym = log2(M);
%         bw_values(mi) = (M * bitrate) / bits_per_sym;
% 
%         nb = floor(n / bits_per_sym) * bits_per_sym;
%         bM = b(1:nb);  % trim bits to multiple of bits_per_sym
% 
%         if isempty(bM)
%             error('Empty bit stream after trimming. Choose larger msg_bits_length.');
%         end
% 
%         % Generate orthogonal MFSK waveforms
%         symM = zeros(samples_per_bit, M);
%         for k = 0:M-1
%             fk = 1 + k;  % frequencies: 1, 2, ..., M
%             wave = sin(2*pi*fk*t);
%             symM(:, k+1) = wave / norm(wave);  % normalized waveforms
%         end
% 
%         % Mapping bits to symbols
%         bM_matrix = reshape(bM, bits_per_sym, []).';
%         idxM = bi2de(bM_matrix, 'left-msb') + 1;
%         tx_signal = symM(:, idxM);  % symbols as columns
%         tx_signal = tx_signal(:).';  % serialize
% 
%         BER = zeros(size(ebn0_range));
%         for i = 1:length(ebn0_range)
%             ebn0_dB = ebn0_range(i);
%             ebn0_lin = 10^(ebn0_dB / 10);
% 
%             % Adjusted SNR for symbol-level AWGN injection
%             snr_sym = ebn0_lin * bits_per_sym * (bitrate / bw_values(mi));
%             snr_sym_dB = 10*log10(snr_sym);
% 
%             rx_signal = awgn(tx_signal, snr_sym_dB, 'measured');
% 
%             % Reshape received signal
%             rx_matrix = reshape(rx_signal, samples_per_bit, []);
% 
%             % Energy detection
%             energies = zeros(M, size(rx_matrix, 2));
%             for m = 1:M
%                 template = symM(:, m);
%                 energies(m, :) = sum((rx_matrix .* template).^2, 1);
%             end
%             [~, id_detected] = max(energies, [], 1);
% 
%             % Demapping
%             rx_bits = de2bi(id_detected - 1, bits_per_sym, 'left-msb').';
%             rx_bits = rx_bits(:).';
% 
%             % BER
%             BER(i) = mean(rx_bits ~= bM);
%         end
% 
%         BER_sim.(['F', num2str(M)]) = BER;
% 
%         % Plot simulated BER
%         semilogy(ebn0_range, BER, '-', 'DisplayName', [num2str(M) '-FSK Sim'], 'Color', colors(mi));
%         hold on;
%     end
% 
%     % Theoretical BER (non-coherent orthogonal MFSK)
%     for mi = 1:length(M_values)
%         M = M_values(mi);
%         bits_per_sym = log2(M);
%         Pb = zeros(size(ebn0_range));
% 
%         for k = 1:length(ebn0_range)
%             EbN0 = 10^(ebn0_range(k) / 10);  % linear scale
%             Pb_k = 0;
%             for n = 1:M-1
%                 coeff = ((M / 2) / (M - 1)) * ((-1)^(n+1) / (n + 1)) * nchoosek(M - 1, n);
%                 Pb_k = Pb_k + coeff * exp((-n * bits_per_sym * EbN0) / (n + 1));
%             end
%             Pb(k) = Pb_k;
%         end
% 
%         BER_th.(['F', num2str(M)]) = Pb;
%         semilogy(ebn0_range, Pb, '--', 'DisplayName', [num2str(M) '-FSK Theory'], 'Color', colors(mi));
%     end
% 
%     % Plot formatting
%     xlabel('E_b/N_0 [dB]');
%     ylabel('Bit Error Rate (BER)');
%     legend('Location', 'southwest');
%     title('Simulated vs. Theoretical BER for Non-coherent MFSK');
%     ylim([1e-7, 1]);
%     grid on;
%     hold off;
%     toc;
% end

