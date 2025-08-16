function FTA_site(sitepath)
% %coorperate with lfp_wavelet_spectro.m
% close all;
% combinepath=fullfile(sitepath,'combined');
% cd(combinepath);
% load p_data.mat
% 
%  savePath = fullfile(sitepath,'results\spectro\wavelet');
% rmdir(savePath,'s'); % 选择
% savePath = fullfile(sitepath,'results\spectro\wavelet\nocut');
% mkdir(savePath);

% coorperate with lfp_wavelet_spectro.m
close all;
combinepath=fullfile(sitepath,'combined');
cd(combinepath);
if isfile('Combined_cut_new.mat')
    load Combined_cut_new.mat
else
    load Combined_cut.mat
end


savePath = fullfile(sitepath,'results\fta\fta30to50');
if isfolder(savePath)
    rmdir(fullfile(sitepath,'results\fta'),'s');
end
mkdir(savePath);



%Set the parameters of the MTM.
Fs = 1000;                  %Sampling frequency
TW = 3;						%Choose time-bandwidth product of 3.
ntapers = 2*TW-1;			%Choose the # of tapers.
fpass = [0 200];
pad = 2;

params.Fs = Fs;				%... sampling frequency
params.tapers=[TW,ntapers];	%... time-bandwidth product & tapers.
params.fpass=fpass;       %... fpass
params.pad = pad;             %... padding.
params.trialave = 1;		%... trial average.

%speparam=params;
%speparam.fpass=[0 120];

%时间窗参数length steplength
%movingwin=[0.3 0.03];
%%
%% 计算channel base要，但是比较channel不运行这段
% aquire 0.2s-1s 实际是1.2s-2s 
recordtime=-1:1/Fs:1.5-1/Fs; % 实际变为0-2.5s
timeindx=recordtime>=0.2 & recordtime <=1 | abs(recordtime-0.2)<eps | abs(recordtime-1)<eps;
needtime=recordtime(timeindx);
vars_lfp = who('-regexp', '^LFP');
vars_spk = who('-regexp', '^SPKC');
for k=1:numel(vars_spk)
    eval(['lfp=',vars_lfp{k},';']);
    eval(['spk=',vars_spk{k},';']);
    for kk=1:numel(spk)
        tmpts=spk(kk).times;
        tmpindx=tmpts>=1.2 & tmpts<=2 | abs(tmpts-1.2)<eps | abs(tmpts-2)<eps;
        tmpts=tmpts(tmpindx);
        % extratime=1.2;%onset变为0
        % spk(kk).times=tmpts-extratime;   
        spk(kk).times=tmpts;   
    end
    lfp=lfp(timeindx,:);
    eval([vars_lfp{k},'=lfp;']);
    eval([vars_spk{k},'=spk;']);
end
disp('now ts is 0.2-1s!(actually 1.2-2s) ')
%disp('now ts is 0.2-1s! and onset is 0 point!!')

%%
Fs=1000;

vars = who;
vars_lfp = who('-regexp', '^LFP');
vars_spk = who('-regexp', '^SPKC');

totallfp_fta=struct();
for i=1:numel(vars_lfp)
    lfp=eval(vars_lfp{i});
    lfp=lfp_norm(lfp); %norm lfp
    spk=eval(vars_spk{i});

    nyquist_freq = Fs/2;
    Wn = [20 40]/nyquist_freq;
    ord = 100;
    b = fir1(ord, Wn);
    K=size(lfp,2);
    N=size(lfp,1);
    FTA = zeros(K, N);
    n=binspikes(spk,Fs,[1.2 2]);
    for k = 1:K
        Vlo = filtfilt(b, 1, lfp(:,k));
        phi = angle(hilbert(Vlo));
        [~, indices] = sort(phi);
        FTA(k, :) = n(indices,k);
    end
    
    %Visualize Field Triggered Average
    figure();plot(linspace(-pi, pi, N), mean(FTA, 1), 'k', 'LineWidth', 2);
    xlabel('Phase');ylabel('FTA');title('FTA of Spikes Using 0-30 Hz Filter');
    xlim([-pi pi]);
    saveas(gcf, fullfile(savePath, [vars_spk{i},'_',vars_lfp{i},'.png']));
    %ylim([0 0.3]);
    % set(gca, 'FontSize', 14)
     totallfp_fta.(vars_spk{i})=FTA;
    %[delta,theta,alpha,beta,gamma]=detaild_filter(s);
end
save(fullfile(savePath,'totallfp_fta20to40.mat'),'-struct','totallfp_fta');
% save(fullfile(savePath,'lowlfp_spec.mat'),'-struct','lowlfp_spec');
