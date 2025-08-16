function lfp_carspec(trialpath,dataname)

close all
load(fullfile(trialpath,dataname));
savepath=fullfile(trialpath,['\',dataname(1:end-4)]);
if isdir(savepath)
    rmdir(savepath,"s");
end

mkdir(savepath);
%%
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
movingwin=[0.3 0.03];
%%
Fs=1000;

vars = who;
vars_lfp = who('-regexp', '^LFP');

for i=1:numel(vars_lfp)
    eval(['lfp=',vars_lfp{i},';']);

        % 找出不全是 NaN 的列
validColumns = ~all(isnan(lfp));

% 只保留不全是 NaN 的列
lfp = lfp(:, validColumns);


    lfp=lfp_norm(lfp); %norm lfp



    spectro=[];

    % multi taper method:
    [spectro,t,frequency]=mtspecgramc(lfp,movingwin,params);
    spectro=spectro';
    meanidx=t>=0.5 &t<=0.75;
    meanspec=spectro(:,meanidx);
    meanspec=mean(meanspec,2);
    spectro=spectro./meanspec;% normalize


    figure;
    set(gcf,'unit','normalized','position',[0,0,1,0.5]);
    subplot(1,2,1);
    imagesc(t-1,frequency,10*log10(spectro));
    %freq = 2.^(round(log2(min(frequency))):round(log2(max(frequency))));
    AX=gca;
    %set(gca,"YScale","log")
    AX.YTickLabelMode = "auto";
        AX.YTick = [0,2,4,8,16,32,64,128,200];;
        hold on
    line([0.2,0.2],ylim,'Color','red','LineWidth',2)
    line([1,1],ylim,'Color','blue','LineWidth',2)
    hold off
    h = colorbar;
    h.Label.String = 'dB';
    xlabel('Time/s');
    xticks(-0.5:0.5:1.5);
    ylabel('Frequency/Hz');
    set(AX, 'YDir', 'normal');
    title([vars_lfp{i},' Total Spectrogram']);

    subplot(1,2,2);
    imagesc(t-1,frequency,10*log10(spectro));
    %freq = 2.^(round(log2(min(frequency))):round(log2(max(frequency))));
    AX=gca;
    %set(gca,"YScale","log")
    AX.YTickLabelMode = "auto";

    hold on
    line([0.2,0.2],ylim,'Color','red','LineWidth',2)
    line([1,1],ylim,'Color','blue','LineWidth',2)
    hold off
    h = colorbar;
    h.Label.String = 'dB';
    xlabel('Time/s');
    xticks(-0.5:0.5:1.5);
    ylabel('Frequency/Hz');
    ylim([0 30]);
    yticks([0,2,4,8,12,16,21,30]);
    set(AX, 'YDir', 'normal');
    title([vars_lfp{i},' 0-30Hz Spectrogram']);
%%
    saveas(gcf,fullfile(savepath,[vars_lfp{i},'.png']));
end