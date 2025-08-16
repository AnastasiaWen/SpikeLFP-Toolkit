function [FTA,phi]=fta(data_spike, data_lfp, fs,wb,t,T)
    %INPUT:
    %   n: 0/1 spike data where rows are trials and columns are timestamps
    %   y: LFP data (mV) where rows are trials and columns are timestamps
    %   wn: fpass-FTA
    %   params: .Fs - sampling freq
    %           .time - [xmin xmax]
    
    
    
    if isstruct(data_spike)
        data_spike=spk_ts2bin(data_spike,t);
    end
    data_lfp=data_lfp';
    
    if exist('T','var')
    smp=1/fs;
    indx = find(t>=T(1));
    indx=indx(1);
    t1 = [T(1):smp:T(2)+eps];
    indx = indx:(indx+length(t1)-1);
    data_lfp = data_lfp(:,indx);
    data_spike=data_spike(:,indx);
    end
    
    %Get Useful Measure
    K = size(data_lfp, 1);%Trial
    N = size(data_lfp,2); 
    sample_freq = fs;
    nyquist_freq = sample_freq/2;
    %t=(params.time(1)):(1/sample_freq):(params.time(2));
    
    %Compute Field Triggered Average
    Wn = wb/nyquist_freq;
    if Wn(1)==0
        [B,A]=butter(3,Wn(2));
    else
        [B,A]=butter(3,Wn);
    end
    FTA = zeros(K, N);
    for k = 1:K
        Vlo = filter(B, A, data_lfp(k, :));
        phi = angle(hilbert(Vlo));
        [~, indices] = sort(phi);
        FTA(k, :) = data_spike(k, indices);
    end
    phi=linspace(-pi, pi, N);
    
    %Visualize Field Triggered Average
%     figure();plot(linspace(-pi, pi, N), mean(FTA, 1), 'k', 'LineWidth', 2);
%     xlabel('Phase');ylabel('FTA');title('FTA of Spikes Using 9-11 Hz Filter');
%     xlim([-pi pi]);ylim([0 0.3]);set(gca, 'FontSize', 14)
    
end