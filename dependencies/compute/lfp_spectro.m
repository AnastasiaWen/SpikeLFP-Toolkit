function lfp_spectro(sitepath)
%计算lfp的spectrogram

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
% if isfile('Combined_cut_new.mat')
%     load Combined_cut_new.mat
% else
    load p_data_cut_inc.mat
% end
% load p_data_detail.mat


savePath = fullfile(sitepath,'results\spectro\multitaper\new');
% if isfolder(fullfile(sitepath,'results\spectro\multitaper\nocut'))
% rmdir(fullfile(sitepath,'results\spectro'),'s');
% end
if isfolder(savePath)
rmdir(fullfile(sitepath,'results\spectro'),'s');
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
movingwin=[0.3 0.03];
%%
Fs=1000;

vars = who;
vars_lfp = who('-regexp', '^LFP');
vars_spk = who('-regexp', '^SPKC');

totallfp_spec=struct();
lowlfp_spec=struct();
% vars=increaseName;
for i=1:numel(vars_lfp)
    tmpname=vars_lfp{i};
    % tmpname=['LFP',tmpname(end-2:end)];
    eval(['lfp=',tmpname,';']);
    colsToRemove = all(isnan(lfp), 1);
    lfp(:, colsToRemove) = [];
    % eval(['lfp=',vars_lfp{i},';']);
    lfp=lfp_norm(lfp); %norm lfp
    
    spectro=[];
    % wavelet method:
    % for j=1:size(lfp,2)
    %     tmplfp=lfp(:,j);
    %     [scalogram,frequency,t]=lfp_wavelet_spectro(tmplfp,[0 200]);
    %     spectro(:,:,j)=scalogram;
    % end
    % spectro=mean(spectro,3);
    % meanspec=spectro(:,501:751);%取-0.5~-0.25s
    % meanspec=mean(meanspec,2);
    % spectro=spectro./meanspec;% normalize

    % multi taper method:
    [spectro,t,frequency]=mtspecgramc(lfp,movingwin,params);
    spectro=spectro';
    meanidx=t>=0.5 &t<=0.75;
    meanspec=spectro(:,meanidx);
    meanspec=mean(meanspec,2);
    spectro=spectro./meanspec;% normalize
    

    eval(['totallfp_spec.',vars_lfp{i},'=spectro;']);
    figure;
    set(gcf,'unit','normalized','position',[0,0,1,0.5]);
    subplot(1,2,1);
    imagesc(t-1,frequency,10*log10(spectro));
    % caxis([-1 1]);
 % colormap('Turbo')
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
    saveas(gcf, fullfile(savePath, ['total_',vars_lfp{i},'.png']));

    % spectro=[];
    % for j=1:size(lfp,2)
    %     tmplfp=lfp(:,j);
    %     [scalogram,frequency,t]=lfp_wavelet_spectro(tmplfp,[0 30]);
    %     spectro(:,:,j)=scalogram;
    % end
    % spectro=mean(spectro,3);
    % figure;
    % imagesc(t-1,frequency,10*log10(spectro));
    % hold on
    % line([0.2,0.2],ylim,'Color','red','LineWidth',2)
    % line([1,1],ylim,'Color','blue','LineWidth',2)
    % hold off
    % h = colorbar;
    % h.Label.String = 'dB';
    % xlabel('Time/s');
    % xticks(-1:0.5:1.5);
    % ylabel('Frequency/Hz');
    % %ylim([0 30]);
    % yticks([4,8,12,16,21,30]);
    % title([vars_lfp{i},' Low-f Spectrogram']);
    % saveas(gcf, fullfile(savePath, ['lowf_',vars_lfp{i},'.png']));
    
    close all;
end
save(fullfile(savePath,'totallfp_spec.mat'),'-struct','totallfp_spec');
% save(fullfile(savePath,'lowlfp_spec.mat'),'-struct','lowlfp_spec');
