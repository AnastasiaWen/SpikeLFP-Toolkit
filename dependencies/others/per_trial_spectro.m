function per_trial_spectro(trialpath)
close all
load(fullfile(trialpath,'data_sd.mat'));
savepath=fullfile(trialpath,'neurospec');
if isdir(savepath)
    rmdir(savepath,"s");
end

mkdir(savepath);
%%
vars = who;
vars_lfp = who('-regexp', '^LFP');
vars_spk = who('-regexp', '^SPKC');

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



%生成timestamp结构体 ( 1 x trial)
for i = 1:numel(vars_spk)
    deltrial=[];
    eval([ 'tmpspk=',vars_spk{i},';']);
    tmpname=vars_spk{i};
    for j=1:numel(tmpspk)
        values = tmpspk(j).times;
        if isempty(values)
            deltrial=[deltrial,j];
        end
    end
    if ~exist(['LFP',tmpname(5:6)])
        continue;
    end
    tmplfp=eval(['LFP',tmpname(5:6)]);
    tmplfp(:,deltrial)=[];
    if isempty(tmplfp)
        continue;
    end
    lfp=tmplfp;

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
    title([vars_spk{i},' Total Spectrogram']);

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
    title([vars_spk{i},' 0-30Hz Spectrogram']);

    saveas(gcf,fullfile(savepath,[vars_spk{i},'.png']));
end