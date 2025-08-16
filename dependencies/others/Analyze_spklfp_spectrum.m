function Analyze_spklfp_spectrum(n, y, params)
    %INPUT:
    %   n: 0/1 spike data where rows are trials and columns are timestamps
    %   y: LFP data (mV) where rows are trials and columns are timestamps
    %   t: timestamp vector
    %   params: .Fs - sampling freq
    %           .time - [xmin xmax]
    
    %Get Useful Measures
    K = size(y, 2);%Trial
    N = size(y,1);
    y=y';
    %sample_interval = t(2) - t(1);
    %sample_freq = 1/sample_interval;
    sample_freq = params.Fs;
    nyquist_freq = sample_freq/2;
    t=(params.time(1)):(1/sample_freq):(params.time(2));
    

    %Get Spectrum
    TW = 3;
    ntapers = 2*TW - 1;
    params.Fs = sample_freq;
    params.tapers = [TW ntapers];
    params.pad = -1;
    params.trialave = 1;
    [C, ~, ~, Syy, Snn, f] = coherencycpb(transpose(y), transpose(n), params);
    
    %Visualize Spectrum
    figure()
    subplot(1, 3, 1)
    plot(f, Snn, 'k', 'LineWidth', 2)
    xlim([0 120])
    xlabel('Frequency')
    ylabel('Power (Hz)')
    title('Snn')
    set(gca, 'FontSize', 14)
    subplot(1, 3, 2)
    plot(f, 10*log10(Syy), 'k', 'LineWidth', 2)
    xlim([0 120])
    xlabel('Frequency')
    ylabel('Power (Hz)')
    title('Syy')
    set(gca, 'FontSize', 14)
    subplot(1,3,3);
    plot(f, C, 'k', 'LineWidth', 2)
    xlim([0 120])
    xlabel('Frequency')
    ylabel('Coherence')
    title('Coherence')
    set(gca, 'FontSize', 14)
end