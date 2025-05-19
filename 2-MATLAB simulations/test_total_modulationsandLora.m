function test_total_modulationsandLora(msg_bits_length)

    tic; %timer to see how long is the simulation
    n = msg_bits_length; %number of bits to simulate
    b = randi([0, 1], 1, n); %generate the binary message b
    samples_per_bit = 30; %how many samples in one bit => oversampling factor
    bitrate = 366.2; % Example bitrate (from LoRa study) => take the precise one
    BW = 125e3; %bandwidth available (in Hz)
    fs = bitrate * samples_per_bit; %calculates the sampling frequency based on the bit rate and the oversampling
    t = 0:1/samples_per_bit:1 - 1/samples_per_bit; %time vector for one symbol duration

    %% BPSK waveform setup
    f = 1; % carrier frequency (normalized to 1 unit)
    sp0 = -sin(2*pi*f*t); %BPSK reference signal for '0' 
    sp0 = sp0 / norm(sp0); %normalized to have unit energy
    sp1 = sin(2*pi*f*t);  %for '1'
    sp1 = sp1 / norm(sp1); %normalized

    %% QPSK waveforms: EXPLAIN MORE
    qpsk_syms = [pi/4, 3*pi/4, 5*pi/4, 7*pi/4]; %4 possible phase angles for QPSK
    num_symbols = n/2;
    qpsk = zeros(1, num_symbols * samples_per_bit);  % preallocate

    for sym_idx = 1:num_symbols
        i = (sym_idx-1)*2 + 1;  % current bit index
        idx = b(i)*2 + b(i+1);
        phase = qpsk_syms(idx + 1);
        symbol = cos(2*pi*f*t + phase);
        qpsk((sym_idx-1)*samples_per_bit + 1 : sym_idx*samples_per_bit) = symbol / norm(symbol);
    end
    %% BPSK modulation
    bpsk = zeros(1, n * samples_per_bit);  % preallocate

    for bit_idx = 1:n
        if b(bit_idx) == 1
            bpsk((bit_idx-1)*samples_per_bit + 1 : bit_idx*samples_per_bit) = sp1;
        else
            bpsk((bit_idx-1)*samples_per_bit + 1 : bit_idx*samples_per_bit) = sp0;
        end
    end


    %% the BER simulation setup
    ebn0_range = -5:1:15; %range of E_b/N_0 values
    ebn0_linear = 10.^(ebn0_range/10); % linear scale
    % preallocates all the BER vectors
    BER_B = zeros(size(ebn0_range));
    BER_Q = zeros(size(ebn0_range));
    BER_F2 = zeros(size(ebn0_range));
    BER_F4 = zeros(size(ebn0_range));
    BER_F8 = zeros(size(ebn0_range));
    BER_F16 = zeros(size(ebn0_range));

    for idx = 1:length(ebn0_range) %the BER loop for each E_b/N_0 point
        ebno = ebn0_range(idx);

        %% BPSK detection
        ebno_adjusted_B = ebno - 10*log10(samples_per_bit); %adjusts the E_b/N_0 because of oversampling (the total energy per bit spread on the samples)
        bpskn = awgn(bpsk, ebno_adjusted_B, 'measured'); %awgn needs per sample: white gaussian noise
        bpskn_matrix = reshape(bpskn, samples_per_bit, []); % each column = one bit segment

        sp1 = sp1(:); % make sure sp1 is a column vector
        sp0 = sp0(:); % make sure sp0 is a column vector
        
        corr_sp1 = sum(bpskn_matrix .* sp1, 1); % correlate with sp1
        corr_sp0 = sum(bpskn_matrix .* sp0, 1); % correlate with sp0
        
        D = double(corr_sp1 > corr_sp0); % decide bit 1 if corr_sp1 > corr_sp0
        BER_B(idx) = mean(D ~= b);

        %% QPSK detection
        ebno_adjusted_Q = ebno + 10*log10(2) - 10*log10(samples_per_bit);
        qpskn = awgn(qpsk, ebno_adjusted_Q, 'measured');
        ref_matrix = zeros(samples_per_bit, 4); % each column is a reference symbol
        for k = 1:4
            ref = cos(2*pi*f*t + qpsk_syms(k));
            ref_matrix(:,k) = ref(:) / norm(ref);
        end
        
        qpskn_matrix = reshape(qpskn, samples_per_bit, []); % reshape the received QPSK signal into columns: samples_per_bit x num_symbols
        corrs = ref_matrix' * qpskn_matrix; % Compute correlations (matrix multiplication) 4 x num_symbols
        [~, idx_syms] = max(corrs, [], 1); % decide based on the max correlation: indices of detected symbols
        Dq = reshape(dec2bin(idx_syms - 1, 2).' - '0', 1, []); % Map detected symbols back to bits: (2 bits per symbol, column major order)
        BER_Q(idx) = mean(Dq ~= b(1:length(Dq)));% Calculate BER

        %% FSK simulations (2-, 4-, 8-, 16-FSK)
        for M = [2, 4, 8, 16] %loop for every MFSK configuration
            bits_per_sym = log2(M); %adjust the number of bits to a multiple of bits_per_sym
            nb = floor(n / bits_per_sym) * bits_per_sym;
            bM = b(1:nb);
            symM = zeros(samples_per_bit, M); % preallocate
            for k = 0:M-1 
                fk = 1 + k;
                symM(:,k+1) = sin(2*pi*fk*t) / norm(sin(2*pi*fk*t)); %normalize
            end
        
            bM_matrix = reshape(bM, bits_per_sym, []).'; % map the bits into fsk symbols: each row = bits for 1 symbol
            idxM = bi2de(bM_matrix, 'left-msb') + 1; % index for each symbol 
            fskM = symM(:,idxM); % select the columns according to idxM
            fskM = fskM(:).'; % flatten to a row vector
            ebno_adjustedM = ebno + 10*log10(bits_per_sym) - 10*log10(samples_per_bit); %add the noise adjusted
            fsknM = awgn(fskM, ebno_adjustedM, 'measured');
        
            fsknM_matrix = reshape(fsknM, samples_per_bit, []); % detection: each column = one symbol
            corrs = symM.' * fsknM_matrix; % correlation matrix: M × num_symbols
            [~, id_detected] = max(corrs, [], 1); % detected symbol indices
            DcfM = de2bi(id_detected-1, bits_per_sym, 'left-msb').'; % convert back to the bits: 2D bits matrix
            DcfM = DcfM(:).'; % flatten into a 1D array

            BER_val = mean(DcfM ~= bM); % Calculate BER
        
            switch M
                case 2
                    BER_F2(idx) = BER_val;
                case 4
                    BER_F4(idx) = BER_val;
                case 8
                    BER_F8(idx) = BER_val;
                case 16
                    BER_F16(idx) = BER_val;
            end
        end

    end

    %% --- LoRa Simulation ---
    fs = BW; %the sampling frequency equals LoRa BW
   
    EbN0_dB = -5:1:15; 
    SF_range = 7:12;
    ber_matrix = zeros(length(SF_range), length(EbN0_dB)); %rows= SF, columns=EbN0 points
    
    for sf_idx = 1:length(SF_range)
        SF = SF_range(sf_idx);
        M = 2^SF;%number of chirps
        Ns = M;%samples per chirp
        Rs = BW / M;%symbol rate
        Rb = Rs * SF; %bit rate

        num_symbols = floor(n / SF);  % number of LoRa symbols from bits
        b_cut = b(1:num_symbols*SF); % keep necessary bits
        b_matrix = reshape(b_cut, SF, []).';% reshape the bits: each row=symbol
        %data_tx = bit2int(b_matrix);    
        data_tx = bi2de(b_matrix, 'left-msb'); % each row= integer 0 to M-1 %REPLACE BY bit2int
    
        t = (0:Ns-1)/fs; %time vector for 1 symbol
        Ts = 1/Rs;
        downchirp = exp(-1j*2*pi*((-BW/2)*t + (BW/(2*Ts))*t.^2));
        base_chirp = exp(1j*2*pi*(BW/(2*Ts))*t.^2); % base chirp (phase common for all chirps)
        f0_all = -BW/2 + (BW/M)*(0:M-1);% frequency steps
    
        for eb_idx = 1:length(EbN0_dB)
            EbN0 = EbN0_dB(eb_idx);
            chirps = base_chirp .* exp(1j*2*pi*f0_all.' * t); % all chirps: size: M x Ns
            tx_lora_matrix = chirps(data_tx+1, :); % select the transmitted chirps based on data each row = one symbol
            %tx_lora = reshape(tx_lora_matrix.', 1, []);
            tx_lora = tx_lora_matrix.'; 
            tx_lora = tx_lora(:).'; % flatten into row vector
    
            Eb = sum(abs(tx_lora).^2) / (num_symbols * SF); %calculate energy per bit
            N0 = Eb / (10^(EbN0/10));
            noise = sqrt(N0/2) * (randn(1, length(tx_lora)) + 1j*randn(1, length(tx_lora)));% Add AWGN
            tx_lora_noisy = tx_lora + noise;
    
            % receiver
            rx_symbols_matrix = reshape(tx_lora_noisy, Ns, []); % each column = one symbol
            dechirped = rx_symbols_matrix .* downchirp.'; % apply downchirp
            fft_out = abs(fft(dechirped)); % size: Ns x num_symbols
            [~, k_hat] = max(fft_out, [], 1);
            data_rx = mod(k_hat-1, M);
            rx_bits_matrix = de2bi(data_rx, SF, 'left-msb'); % map the detected symbols back to bits (vectorized)
            rx_bits = rx_bits_matrix.';
            rx_bits = rx_bits(:).'; % flatten
    
            % BER
            bit_errors = sum(rx_bits ~= b_cut);
            ber_matrix(sf_idx, eb_idx) = bit_errors / length(b_cut);
        end
    end
        
        

    %% Plot results
    figure;
    semilogy(ebn0_range, BER_F2, 'r-', ...
             ebn0_range, BER_F4, 'b-', ...
             ebn0_range, BER_F8, 'g-', ...
             ebn0_range, BER_F16, 'c-', ...
             ebn0_range, BER_B, 'k--', ...
             ebn0_range, BER_Q, 'm--', ...
             ebn0_range, ber_matrix(1,:), 'o--', ...
             ebn0_range, ber_matrix(6,:), 'x--', 'LineWidth', 1)

             %ebn0_range, ber_matrix(2,:), 'x--', ...
             %ebn0_range, ber_matrix(3,:), 'o--', ...
             %ebn0_range, ber_matrix(4,:), 'x--', ...
             %ebn0_range, ber_matrix(5,:), 'o--', ...
    grid on;
    xlabel('E_b/N_0 (dB)');
    ylabel('Bit Error Rate (BER)');
    legend('2-FSK', '4-FSK', '8-FSK', '16-FSK', 'BPSK', 'QPSK', 'LoRa SF7', 'LoRa SF12');
    title('BER vs. E_b/N_0 for Modulation Schemes');
    toc;
end


%% comments
 %num_symbols = msg_length;
    % qpsk = [];
    % for i = 1:2:n %group of bit in pairs
    %     idx = b(i)*2 + b(i+1); 
    %     phase = qpsk_syms(idx + 1); 
    %     symbol = cos(2*pi*f*t + phase);
    %     qpsk = [qpsk symbol / norm(symbol)];
    % end

        % bpsk = [];
    % for i = 1:n
    %     bpsk = [bpsk (b(i) == 1)*sp1 + (b(i) == 0)*sp0];
    % end
    
            % D = [];
        % for i = 1:samples_per_bit:length(bpskn)
        %     segment = bpskn(i:i+samples_per_bit-1); %segments the received signal
        %     bit = double(sum(segment .* sp1) > sum(segment .* sp0)); %correlates with the references sp0 and sp1
        %     D = [D bit]; %decision of the bit
        % end
        % BER_B(idx) = sum(D ~= b) / n;




        % snr_Q = ebno + 10*log10(2) - 10*log10(samples_per_bit);
        % qpskn = awgn(qpsk, snr_Q, 'measured');
        % Dq = [];
        % for i = 1:samples_per_bit:length(qpskn)
        %     segment = qpskn(i:i+samples_per_bit-1);
        %     corrs = zeros(1,4);
        %     for k = 1:4
        %         ref = cos(2*pi*f*t + qpsk_syms(k));
        %         corrs(k) = sum(segment .* (ref / norm(ref)));
        %     end
        %     [~, idx_sym] = max(corrs);
        %     bits = dec2bin(idx_sym - 1, 2) - '0';
        %     Dq = [Dq bits];
        % end
        % BER_Q(idx) = sum(Dq ~= b(1:length(Dq))) / length(Dq);


            %         symM = cell(M,1);
    %         fskM = [];
    %         for k = 0:M-1 
    %             fk = 1 + k;
    %             symM{k+1} = sin(2*pi*fk*t) / norm(sin(2*pi*fk*t));
    %         end
    %         for i = 1:bits_per_sym:nb %map the bits into FSK symbols
    %             idxM = bin2dec(num2str(bM(i:i+bits_per_sym-1)));
    %             fskM = [fskM symM{idxM+1}];
    %         end
    %         ebno_adjustedM = ebno + 10*log10(bits_per_sym) - 10*log10(samples_per_bit); %noise
    %         fsknM = awgn(fskM, ebno_adjustedM, 'measured'); %to make more efficient, variable changes size each op
    %         DcfM = [];
    %         for i = 1:samples_per_bit:length(fsknM) %detect with correlation
    %             seg = fsknM(i:i+samples_per_bit-1);
    %             corr = cellfun(@(x) sum(seg .* x), symM);
    %             [~, id] = max(corr);
    %             bits = dec2bin(id-1, bits_per_sym) - '0';
    %             DcfM = [DcfM bits];
    %         end
    %         BER_val = sum(DcfM ~= bM) / length(bM); %calculate ber
    %         switch M
    %             case 2
    %                 BER_F2(idx) = BER_val;
    %             case 4
    %                 BER_F4(idx) = BER_val;
    %             case 8
    %                 BER_F8(idx) = BER_val;
    %             case 16
    %                 BER_F16(idx) = BER_val;
    %         end
    %     end



    %     data_tx = zeros(1, num_symbols);
    %     for i = 1:num_symbols
    %         bits_chunk = b((i-1)*SF + 1 : i*SF);
    %         data_tx(i) = bin2dec(num2str(bits_chunk));
    %     end
    % 
    %     t = (0:Ns-1)/fs;
    %     Ts = 1/Rs;
    %     downchirp = exp(-1j*2*pi*((-BW/2)*t + (BW/(2*Ts))*t.^2));
    % 
    %     for eb_idx = 1:length(EbN0_dB)
    %         EbN0 = EbN0_dB(eb_idx);
    %         tx_lora = [];
    %         for i = 1:num_symbols
    %             k = data_tx(i);
    %             f0 = -BW/2 + k * (BW / M);
    %             chirp = exp(1j*2*pi*(f0*t + (BW/(2*Ts))*t.^2));
    %             tx_lora = [tx_lora chirp]; %not efficient, variable changes size each operation
    %         end
    %         Eb = sum(abs(tx_lora).^2) / (num_symbols * SF);
    %         N0 = Eb / (10^(EbN0/10));
    %         noise = sqrt(N0/2) * (randn(1, length(tx_lora)) + 1j * randn(1, length(tx_lora)));
    %         tx_lora_noisy = tx_lora + noise;
    % 
    %         data_rx = zeros(1, num_symbols);
    %         for i = 1:num_symbols
    %             rx_symbol = tx_lora_noisy((i-1)*Ns+1 : i*Ns);
    %             dechirped = rx_symbol .* downchirp;
    %             fft_out = abs(fft(dechirped));
    %             [~, k_hat] = max(fft_out);
    %             data_rx(i) = mod(k_hat-1, M);
    %         end
    %         errors = sum(data_rx ~= data_tx);
    %         % ber_matrix(sf_idx, eb_idx) = errors / (num_symbols * SF);
    %         %ber_matrix(sf_idx, eb_idx) = (errors * SF) / (num_symbols * SF);  % same as errors / num_symbols
    %         % Reconstruct bits from symbols
    %         rx_bits = [];
    %         for i = 1:num_symbols
    %             rx_bits = [rx_bits dec2bin(data_rx(i), SF) - '0'];
    %         end
    %         orig_bits = b(1:num_symbols * SF);
    %         bit_errors = sum(rx_bits ~= orig_bits);
    %         ber_matrix(sf_idx, eb_idx) = bit_errors / length(orig_bits);
    % 
    % 
    %     end
    % end