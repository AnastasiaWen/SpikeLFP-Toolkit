function LFP=lfpextract(data,plotname,opn)
Fs=1000;
%% downsample
data=downsample(data,40); %40000-1000
wb=data;

%% optional notch
if exist('opn')
    for i=1:length(opn)
        ffilt = opn(i); % 欲消除的频率为50Hz
        Q = 100; % 带宽参数
        w = ffilt/(Fs/2);
        bw = w/Q;
        [b, a] = iirnotch(w, bw);
        data=filtfilt(b,a,data);
    end
end

%% 50Hz notch
ffilt = 50; % 欲消除的频率为50Hz
Q = 100; % 带宽参数
w = ffilt/(Fs/2);
bw = w/Q;
[b, a] = iirnotch(w, bw);
data=filtfilt(b,a,data);


ffilt = 100; % 欲消除的频率为50Hz
Q = 100; % 带宽参数
w = ffilt/(Fs/2);
bw = w/Q;
[b, a] = iirnotch(w, bw);
data=filtfilt(b,a,data);

ffilt = 150; % 欲消除的频率为50Hz
Q = 100; % 带宽参数
w = ffilt/(Fs/2);
bw = w/Q;
[b, a] = iirnotch(w, bw);
data=filtfilt(b,a,data);

ffilt = 200; % 欲消除的频率为50Hz
Q = 100; % 带宽参数
w = ffilt/(Fs/2);
bw = w/Q;
[b, a] = iirnotch(w, bw);
data=filtfilt(b,a,data);

ffilt = 250; % 欲消除的频率为50Hz
Q = 100; % 带宽参数
w = ffilt/(Fs/2);
bw = w/Q;
[b, a] = iirnotch(w, bw);
data=filtfilt(b,a,data);

%%
% [b,a]=designNotchPeakIIR(FilterOrder=2,...
%     Response='notch',CenterFrequency=w,QualityFactor=100);
%% 带通滤波 0.5Hz-300Hz
d = designfilt('bandpassfir', 'FilterOrder', 300, ...
       'CutoffFrequency1', 0.5, 'CutoffFrequency2', 300,...
       'SampleRate', Fs,'Window','hamming');
% d=designfilt('bandpassiir', ...
%     'FilterOrder',4,'HalfPowerFrequency1',0.5, ...
%     'HalfPowerFrequency2',300,'SampleRate',Fs);
data=filtfilt(d,data);

%% detrend
LFP=detrend(data);

LFP_before_emd=LFP;
%% EMD分离
Fs = 1000; % 采样频率

L = length(LFP); % 数据长度
t = 0:L-1; % 时间向量
dt=1/Fs;
t=t/Fs;
T=t(end);

imf=emd(LFP);
emd_visu(LFP(1:10000),t(1:10000),imf(:,1:10000),1);  % EMD专用画图函数
saveas(gcf, ['channel',plotname(3:end),'_emd.png']);
%% EMD合成
nimf=size(imf,1)-1;
figure;
tiledlayout(ceil(nimf/4),4);
set(gcf,'unit','normalized','position',[0.2,0.2,0.64,1]);
for i=1:nimf
     [faxis,Sxx]=fftspecc(imf(i,:)');
     findx=faxis<=300;
     faxis=faxis(findx);
     Sxx=Sxx(findx);
     [peakS,peakindx]=max(Sxx);
     if faxis(peakindx)>1
        emd_inv_indx=i;
     end
     nexttile(i);
     plot(faxis,10*log10(Sxx));
     xlabel('Freq(Hz)')
     ylabel('dB')
     xlim([0 300]);
     ylim([-150 0]);
     title(['imf',num2str(i)])
end
sgtitle(['Use imf 1:',num2str(emd_inv_indx)]);
saveas(gcf, ['channel',plotname(3:end),'_lfp_emd_inv_spectrum.png']);
LFP=sum(imf(1:emd_inv_indx,:),1);

%% plot
[faxis,Sxx]=fftspecc(wb(1:10000));
figure;
%befre
subplot(3,2,1);
plot(faxis,10*log10(Sxx)) 
xlim([0 200]);ylim([-150 0])
title('Spectrum WB')
xlabel('Freq (Hz)')
ylabel('dB')

subplot(3,2,2);
plot(t(1:10000),wb(1:10000));
title('WB')
xlabel('time(s)')
ylabel('Volt')

[faxis,Sxx]=fftspecc(LFP_before_emd(1:10000));

%befre
subplot(3,2,3);
plot(faxis,10*log10(Sxx)) 
xlim([0 200]);ylim([-150 0])
title('Spectrum LFP before')
xlabel('Freq (Hz)')
ylabel('dB')

subplot(3,2,4);
plot(t(1:10000),LFP_before_emd(1:10000));
title('LFP befre')
xlabel('time(s)')
ylabel('Volt')

[faxis,Sxx]=fftspecc(LFP(1:10000));
subplot(3,2,5);
plot(faxis,10*log10(Sxx)) 
xlim([0 200]);ylim([-150 0])
title('Spectrum LFP after')
xlabel('Freq (Hz)')
ylabel('dB')

subplot(3,2,6);
plot(t(1:10000),LFP(1:10000));
title('LFP after')
xlabel('time(s)')
ylabel('Volt')

saveas(gcf, ['channel',plotname(3:end),'_spectrum_compare.png']);
end
