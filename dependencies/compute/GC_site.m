function GC_site(sitepath)
%granger causality计算同channel的spk和lfp结果


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
    load p_data_cut_inc.mat




savePath = fullfile(sitepath,'results\GC');
%if isfolder(fullfile(sitepath,'results\spectro\multitaper\nocut'))
% rmdir(fullfile(sitepath,'results\spectro'),'s');
%end
if isfolder(savePath)
    rmdir(savePath,'s');
end
mkdir(savePath);


%Set the parameters of the MTM.
Fs = 1000;                  %Sampling frequency
TW = 3;						%Choose time-bandwidth product of 3.
ntapers = 2*TW-1;			%Choose the # of tapers.
fpass = [0 50];
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

%% aquire -0.5s-1.5s 实际是0.5s-2.5s 
recordtime=-1:1/Fs:1.5-1/Fs; % 实际变为0-2.5s
timeindx=recordtime>=-0.5 & recordtime <=1.5 | abs(recordtime-0.5)<eps | abs(recordtime-1.5)<eps;
needtime=recordtime(timeindx);
vars_lfp = who('-regexp', '^LFP');
vars_spk = who('-regexp', '^SPKC');
% for k=1:numel(vars_spk)
%     eval(['lfp=',vars_lfp{k},';']);
%     eval(['spk=',vars_spk{k},';']);
%     for kk=1:numel(spk)
%         tmpts=spk(kk).times;
%         tmpindx=tmpts>=0.5 & tmpts<=2.5 | abs(tmpts-0.5)<eps | abs(tmpts-2.5)<eps;
%         tmpts=tmpts(tmpindx);
%         extratime=0.5;%onset变为0
%         spk(kk).times=tmpts-extratime;   
%         %spk(kk).times=tmpts;   
%     end
%     lfp=lfp(timeindx,:);
%     eval([vars_lfp{k},'=lfp;']);
%     eval([vars_spk{k},'=spk;']);
% end
for k=1:numel(vars_spk)
    eval(['spk=',vars_spk{k},';']);
    for kk=1:numel(spk)
        tmpts=spk(kk).times;
        tmpindx=tmpts>=0.5 & tmpts<=2.5 | abs(tmpts-0.5)<eps | abs(tmpts-2.5)<eps;
        tmpts=tmpts(tmpindx);
        extratime=0.5;%onset变为0
        spk(kk).times=tmpts-extratime;   
        %spk(kk).times=tmpts;   
    end;
    eval([vars_spk{k},'=spk;']);
end
for k=1:numel(vars_lfp)
    eval(['lfp=',vars_lfp{k},';']);
    lfp=lfp(timeindx,:);
    eval([vars_lfp{k},'=lfp;']);
end
disp('now ts is -0.5-1.5s!(actually 0.5-2.5s) ')
%%

totalGCs2p=struct();
totalGCp2s=struct();
for i=1:numel(vars_lfp)
    eval(['lfp=',vars_lfp{i},';']);
    lfp=lfp_norm(lfp); %norm lfp
    spk=eval(vars_spk{i});
    [~,~,~,~,tt,f,GCs2p,GCp2s,ctl]=spk_lfp_GC(lfp,spk,params,movingwin);
    if ctl==1
        continue;
    end
    eval(['totalGCs2p.',vars_spk{i},'=GCs2p;']);
    eval(['totalGCp2s.',vars_spk{i},'=GCp2s;']);
    figure;
    subplot(1,2,1);
    plot(tt-0.5, sum(GCs2p,2));
    %title([vars_spk{i},'+',vars_lfp{i},'Time-GC']);
    xlabel('time (s)');ylabel('GC');hold on;
    plot(tt-0.5, sum(GCp2s,2));
    %title([vars_spk{i},'+',vars_lfp{i}, 'Time-GC']);
     xlabel('time (s)');ylabel('GC');hold off;
    legend('GCsp2','GCp2s');
    xlim([-0.1 1.1]);
    xticks(-0.5:0.5:1.5);
    subplot(1,2,2);
    plot(f, sum(GCs2p,1));
    %title([vars_spk{i},'+',vars_lfp{i},'Freq-GC']);
    xlabel('Frequency (Hz)');ylabel('GC');hold on;
    plot(f, sum(GCp2s,1));
    %title([vars_spk{i},'+',vars_lfp{i}, 'Freq-GC']);
    xlabel('Frequency (Hz)');ylabel('GC');hold off;
    %xticks(-0.5:0.5:1.5);
    xlim([0 100]);
    legend('GCsp2','GCp2s');

    
    saveas(gcf, fullfile(savePath, [vars_spk{i},'_',vars_lfp{i},'.png']));
    
    close all;
end
savename={'tt','f','totalGCp2s','totalGCs2p'}
save(fullfile(savePath,'totallGC.mat'),savename{:});
% save(fullfile(savePath,'lowlfp_spec.mat'),'-struct','lowlfp_spec');
