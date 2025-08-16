% Description:
%   this script is used to view all wb spectrum in a .pl2
%   to find whether there is any noise needed to notch

close all
WBnames={'WB02','WB03','WB04','WB05','WB06','WB07','WB08','WB09','WB10','WB11','WB12','WB13','WB14','WB15','WB16'};
for i=1:length(WBnames)
    eval(['data=',WBnames{i},';']);
    Fs=1000;
%% downsample
data=downsample(data,40); %40000-1000
wb=data;

%% optional notch
%% optional notch
% ffilt = 50.3; % 欲消除的频率为50Hz
% Q = 100; % 带宽参数
% w = ffilt/(Fs/2);
% bw = w/Q;
% [b, a] = iirnotch(w, bw);
% data=filtfilt(b,a,data);
% 
% ffilt = 69.1; % 欲消除的频率为50Hz
% Q = 100; % 带宽参数
% w = ffilt/(Fs/2);
% bw = w/Q;
% [b, a] = iirnotch(w, bw);
% data=filtfilt(b,a,data);
% 
% ffilt = 110.8; % 欲消除的频率为50Hz
% Q = 100; % 带宽参数
% w = ffilt/(Fs/2);
% bw = w/Q;
% [b, a] = iirnotch(w, bw);
% data=filtfilt(b,a,data);
% 
% ffilt = 167.9; % 欲消除的频率为50Hz
% Q = 100; % 带宽参数
% w = ffilt/(Fs/2);
% bw = w/Q;
% [b, a] = iirnotch(w, bw);
% data=filtfilt(b,a,data);
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

[faxis,Sxx]=fftspecc(wb(1:10000));
figure;
%befre
subplot(2,2,1);
plot(faxis,10*log10(Sxx)) 
xlim([0 200]);ylim([-150 0])
title('Spectrum WB')
xlabel('Freq (Hz)')
ylabel('dB')

subplot(2,2,2);
plot(t(1:10000),wb(1:10000));
title('WB')
xlabel('time(s)')
ylabel('Volt')

[faxis,Sxx]=fftspecc(LFP_before_emd(1:10000));

%befre
subplot(2,2,3);
plot(faxis,10*log10(Sxx)) 
xlim([0 200]);ylim([-150 0])
title('Spectrum LFP before')
xlabel('Freq (Hz)')
ylabel('dB')

subplot(2,2,4);
plot(t(1:10000),LFP_before_emd(1:10000));
title('LFP befre')
xlabel('time(s)')
ylabel('Volt')

end

% %%
% clear;clc;
% opn=[76.5,79.3];
% trialpath='G:\2023_Mn_14325_RightBrain\Week07\D2_site04\D2_site04_TRIAL01';
% lfp_extract_only1(trialpath,opn);
% 
