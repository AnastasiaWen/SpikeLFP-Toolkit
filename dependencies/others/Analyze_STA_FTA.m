function Analyze_STA_FTA(n, y, params)
    %INPUT:
    %   n: 0/1 spike data where rows are trials and columns are timestamps
    %   y: LFP data (mV) where rows are trials and columns are timestamps
    %   wn: fpass-FTA
    %   params: .Fs - sampling freq
    %           .time - [xmin xmax]
    
    %Get Useful Measures
    K = size(y, 2);%Trial
    N = size(y,1);
    y=y';
    sample_freq = params.Fs;
    nyquist_freq = sample_freq/2;
    t=(params.time(1)):(1/sample_freq):(params.time(2));
    

    %Compute Spike-Triggered Average
    win = 100;
    STA = zeros(K, 2*win+1);
    for k = 1:K
        spks = find(n(k, :) == 1);
        for i = 1:length(spks)
            if (spks(i) > win & spks(i) + win < N)
                STA(k, :) = STA(k, :) + y(k, spks(i) - win:spks(i) + win)/length(spks);
            end
        end
    end
    
    
    %Visualize Spike-Triggered Average
    figure();plot(-100:100, STA, 'k', 'LineWidth', 2)
    xlabel('Time (milliseconds)');ylabel('Voltage (mV)');title('Spike-Triggered per trial');
    set(gca, 'FontSize', 14)
    figure();plot(-100:100, mean(STA, 1), 'k', 'LineWidth', 2);
    xlabel('Time (milliseconds)');ylabel('Voltage (mV)');title('Spike-Triggered Average');
    set(gca, 'FontSize', 14)
    
    %Compute Field Triggered Average
    
    Wn = [44 46]/nyquist_freq;
    ord = 100;
    b = fir1(ord, Wn);
    FTA = zeros(K, N);
    for k = 1:K
        Vlo = filtfilt(b, 1, y(k, :));
        phi = angle(hilbert(Vlo));
        [~, indices] = sort(phi);
        FTA(k, :) = n(k, indices);
    end
    
    %Visualize Field Triggered Average
    figure();plot(linspace(-pi, pi, N), mean(FTA, 1), 'k', 'LineWidth', 2);
    xlabel('Phase');ylabel('FTA');title('FTA of Spikes Using 9-11 Hz Filter');
    xlim([-pi pi]);ylim([0 0.3]);set(gca, 'FontSize', 14)
    
end