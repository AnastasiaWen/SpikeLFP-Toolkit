function channel_cohegram(sitepath)
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




savePath = fullfile(sitepath,'results\channel_cohe_ch');
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
if numel(vars_lfp)==1
    return
end
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
totalC=struct();
totalphi=struct();
num_permutations=1000;
for i=1:numel(vars_lfp)
    eval(['lfp=',vars_lfp{i},';']);
    lfp=lfp_norm(lfp); %norm lfp
    lfporder_name=vars_lfp{i};
    lfporder_name=['Order',lfporder_name(end-1:end)];
    lfporder=eval(lfporder_name);

    tmpname=vars_lfp{i};tmpname=tmpname(end-1:end);tmpname=['^SPKC',tmpname];
    vars_sp = vars_spk(~cellfun('isempty', regexp(vars_spk, tmpname)));

    for ss=1:numel(vars_sp)
        SPKCdata=eval(vars_sp{ss});
        tmpname=vars_sp{ss};tmpname=tmpname(end-2:end);
        orderspk=eval(['Order',tmpname]);
        [tmplfp,tmpspk] = match_trials(lfp, SPKCdata, lfporder, orderspk);
        if isempty(tmplfp) | isempty(tmpspk)
            continue;
        end
        % if i==1& ss==1
            % [tt,f,p_value_fisher, GCs2p,GCp2s] = permutation_test_granger(tmplfp, tmpspk, num_permutations,1);
        % else
        % [tt,f,p_value_fisher, GCs2p,GCp2s] = permutation_test_granger(tmplfp, tmpspk, num_permutations,0);
        % end
        [C,phi,~,~,~,tt,f]=cohgramcpt(tmplfp, tmpspk,movingwin,params)
        eval(['totalC.',vars_lfp{i},'_',vars_sp{ss},'=C;']);
        eval(['totalphi.',vars_lfp{i},'_',vars_sp{ss},'=phi;']);
        % eval(['pGCp2s.',vars_lfp{i},'_',vars_sp{ss},'=p_value_fisher(1);']);
        % eval(['pGCs2p.',vars_lfp{i},'_',vars_sp{ss},'=p_value_fisher(2);']);
        figure;
         subplot(1,2,1);
        plot(tt-0.5, mean(C,2));
        %title([vars_spk{i},'+',vars_lfp{i},'Time-GC']);
        xlabel('time (s)');ylabel('C');hold on;
        subplot(1,2,2);
        plot(tt-0.5, mean(phi,2));
        %title([vars_spk{i},'+',vars_lfp{i}, 'Time-GC']);
         xlabel('time (s)');ylabel('Phi');hold off;

        xlim([-0.2 1.2]);
        xticks([0:0.2:1]);
        % title(['pGCs2p=',num2str(p_value_fisher(1)),',','pGCp2s=',num2str(p_value_fisher(2))])
        % subplot(1,2,2);
        % plot(f, sum(GCs2p,1));
        % %title([vars_spk{i},'+',vars_lfp{i},'Freq-GC']);
        % xlabel('Frequency (Hz)');ylabel('GC');hold on;
        % plot(f, sum(GCp2s,1));
        % %title([vars_spk{i},'+',vars_lfp{i}, 'Freq-GC']);
        % xlabel('Frequency (Hz)');ylabel('GC');hold off;
        % %xticks(-0.5:0.5:1.5);
        % xlim([0 100]);
        % legend('GCsp2','GCp2s');
    
        
        saveas(gcf, fullfile(savePath, [vars_lfp{i},'_',vars_sp{ss},'.png']));
       
    end
    close all;
end
if ~exist('tt')
    return
end
savename={'tt','f','totalC','totalphi'};
save(fullfile(savePath,'totalCoh.mat'),savename{:});
% save(fullfile(savePath,'lowlfp_spec.mat'),'-struct','lowlfp_spec');
