function test_psk_only(msg_bits_length)
    tic; %for the elapsed time
    n = msg_bits_length; %length of message
    b = randi([0, 1], 1, n); %random bit stream
    samples_per_bit = 16; %sampling resolution per bit SHOULD BE THAT!!
    %bits per symbol should be 1 or 2 according to BPSK or QPSK...
    bitrate = 366.2;
    %fs = bitrate * samples_per_bit;
    t = 0:1/samples_per_bit:1 - 1/samples_per_bit; %time vector for 1 bit/symbol duration
    f = 1;%868e6;%carrier frequency
    bw_b=20*bitrate; %assumption bandwidth bpsk
    bw_q=10*bitrate;%bandwidth qpsk

    %% BPSK signal
    sp0 = -sin(2*pi*f*t); %bit 0
    sp0 = sp0 / sqrt(sum(sp0.^2)); %normalized
    sp1 =  sin(2*pi*f*t); %bit 1
    sp1 = sp1 / sqrt(sum(sp1.^2));%normalized
    bpsk = zeros(1, n * samples_per_bit); %preallocation
    for bit_idx = 1:n %assign sp0 or sp1 waveform for each bit
        if b(bit_idx) == 1
            bpsk((bit_idx-1)*samples_per_bit + 1 : bit_idx*samples_per_bit) = sp1;
        else
            bpsk((bit_idx-1)*samples_per_bit + 1 : bit_idx*samples_per_bit) = sp0;
        end
    end
    %% QPSK signal
    qpsk_syms = [pi/4, 3*pi/4, 5*pi/4, 7*pi/4]; %gray coded phase mapping for qpsk
    num_symbols = floor(n/2); 
    b_qpsk = b(1:num_symbols*2); %needs even number of bits: 2 bits/symbol

    qpsk = zeros(1, num_symbols * samples_per_bit); %preallocation qpsk vector
    for sym_idx = 1:num_symbols
        i = (sym_idx-1)*2 + 1;
        idx = b_qpsk(i)*2 + b_qpsk(i+1); %bit pair
        phase = qpsk_syms(idx + 1); %map each bit pair to a qpsk phase
        symbol = cos(2*pi*f*t + phase);  % construct signal for one symbol
        qpsk((sym_idx-1)*samples_per_bit + 1 : sym_idx*samples_per_bit) = symbol;
    end

     qpsk= qpsk/norm(qpsk) ;%* 1000*sqrt(2); %normalize qpsk to unit power




    %% BER simulation
    ebn0_range = -5:1:15;
    BER_B = zeros(size(ebn0_range));
    BER_Q = zeros(size(ebn0_range));

    for idx = 1:length(ebn0_range)
        ebno = ebn0_range(idx);

        % BPSK
        %ebno_adj_B = ebno - 10*log10(samples_per_bit);
        snr_B = (10^((ebno)/10))*(bitrate/bw_b); %converts EbN0 to SNR using specific bandwidth
        snr_adj_B= 10*log10(snr_B);%- 10*log10(samples_per_bit);%-10*log10(samples_per_bit)
       
        bpskn = awgn(bpsk, snr_adj_B, 'measured'); %add awgn to the bpsk signal
        bpskn_matrix = reshape(bpskn, samples_per_bit, []); %reshape
        corr_sp1 = sum(bpskn_matrix .* sp1(:), 1); %correlation with sp0 and sp1
        corr_sp0 = sum(bpskn_matrix .* sp0(:), 1);
        D_B = double(corr_sp1 > corr_sp0); %choose the highest correlation
        BER_B(idx) = mean(D_B ~= b);%compare with the b transmitted bits

        % QPSK 
        
        snr_Q = (10^(ebno/10))*(bitrate/bw_q);
        snr_adj_Q= 10*log10(snr_Q)+10*log10(2); %same snr conversion
        %ebno_adj_Q = ebno + 10*log10(2) - 10*log10(samples_per_bit);
        qpskn = awgn(qpsk, snr_adj_Q, 'measured'); %noise addition
        ref_matrix = zeros(samples_per_bit, 4); %reference waveforms for each qpsk phase
        for k = 1:4
            ref = cos(2*pi*f*t + qpsk_syms(k));
            ref_matrix(:,k) = ref(:) / norm(ref);
        end

        qpskn_matrix = reshape(qpskn, samples_per_bit, []);
        corrs = ref_matrix' * qpskn_matrix; %correlate received signal with each reference
        [~, idx_syms] = max(corrs, [], 1); %choose max index
        Dq = reshape(dec2bin(idx_syms - 1, 2).' - '0', 1, []); %convert indices back to bit pair
        BER_Q(idx) = mean(Dq ~= b_qpsk(1:length(Dq))); %compute ber
    end

    %% theoretical BER 
    EbN0_dB = -5:1:15;
    EbN0_lin = 10.^(EbN0_dB / 10);
    %snr_lin=(10^(ebno/10))*bitrate/bw_b;
    Pb_BPSK = 0.5 * erfc(sqrt(EbN0_lin)); % same for BPSK and QPSK under AWGN
    

    %% plot
    figure;
    semilogy(EbN0_dB, Pb_BPSK, 'DisplayName', 'Theoretical BPSK/QPSK');
    hold on;
    semilogy(ebn0_range, BER_B,'DisplayName', 'Sim BPSK');
    semilogy(ebn0_range, BER_Q,  'DisplayName', 'Sim QPSK');
    grid on;
    xlabel('E_b/N_0 (dB)');
    ylabel('BER (P_b)');
    title('BPSK and QPSK BER Theory and Simulation');
    legend show;
    ylim([1e-7 1]);
    xlim([-5 20]);
    hold off;
    toc;
end
