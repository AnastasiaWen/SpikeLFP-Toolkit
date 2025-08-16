function LFP_car(trialpath)
close all;
cd(trialpath);
savepath=fullfile(trialpath,'carsig');
mkdir(savepath);
load data_delch.mat
%%
Fs=1000;
detail_change=0;
vars = who;
vars_lfp = who('-regexp', '^LFP');
vars_spk = who('-regexp', '^SPKC');
vars_ang = who('-regexp', '^StimAngle');
vars_spd = who('-regexp', '^StimSpeed');
vars_wf = who('-regexp', '^Waveform');
vars_wft = who('-regexp', '^Wavetime');

outliers=[];
    d = designfilt('bandpassfir', 'FilterOrder', 300, ...
       'CutoffFrequency1', 1, 'CutoffFrequency2', 300,...
       'SampleRate', Fs,'Window','hamming');
for i=1:size(LFP02,2)
    close all;
    com=[];
    for j=1:numel(vars_lfp)
    LFP_data=eval(vars_lfp{j});
    LFP_data=LFP_data(:,i);
    if ~isnan(LFP_data)
    % LFP_data=filtfilt(d,LFP_data);
    %  eval([vars_lfp{j},'(:,i)=LFP_data;'])
            com=[com,LFP_data];
    end

    end


    % hw=figure;
    % set(gcf,'unit','normalized','position',[0,0,0.4,1]);
    % ha=tight_subplot(8,2,0.04);
    % 
    % % figure;
    % hold on
    % cm=[];
    % for q=1:numel(vars_lfp)
    %     eval(['lfp=',vars_lfp{q},';']);
    %     % plot(lfp(:,1));
    %     cm=[cm,lfp(:,i)];
    % 
    %     axes(ha(q));
    %     plot(-1:1/1000:1.4999,lfp(:,i));
    %         title(vars_lfp{q},'Color','black');
    %         ylim([-0.18 0.18]);
    % end
    % 
    % saveas(gcf,fullfile(savepath,['trial',num2str(i),'_orig.png']));
    % 
    % 
    % close all



    mean_values = mean(com, 2);             % 计算每个 trial 的均值 (1x120)
    std_values = std(com, 0, 2);            % 计算每个 trial 的标准差 (1x120)
    
    % 定义阈值范围
    lower_bound = mean_values - 10 * std_values;  % 下限
    upper_bound = mean_values + 10 * std_values;  % 上限
    
    % 找出超出范围的 trial
    outliers = false(1, size(com, 2));      % 创建逻辑数组标记异常trial

    for j=1:size(com,2)
        tmpdata=com(:,j);
         % if any(tmpdata(:) < lower_bound(j) | tmpdata(:) > upper_bound(j))
         %    outliers(j) = true;
         % end
         if any(tmpdata(:) < -0.4 | tmpdata(:) > 0.4 )
            outliers(j) = true;
         end
    end
         
    removed_trials_index = find(outliers);

    com(:,removed_trials_index)=[];

    
    commean=mean(com,2);
    if isempty(com)
        commean=zeros(2500,1);
    end

    for j=1:numel(vars_lfp)
    LFP_data=eval(vars_lfp{j});
    LFP_data=LFP_data(:,i);
    if ~isnan(LFP_data)
    eval([vars_lfp{j},'(:,i)=',vars_lfp{j},'(:,i)-commean;'])
    end
    end


% hw=figure;
% set(gcf,'unit','normalized','position',[0,0,0.4,1]);
% ha=tight_subplot(8,2,0.04);
% 
% % figure;
% hold on
% cm=[];
% for q=1:numel(vars_lfp)
%     eval(['lfp=',vars_lfp{q},';']);
%     % plot(lfp(:,1));
%     cm=[cm,lfp(:,i)];
% 
%     axes(ha(q));
%     plot(-1:1/1000:1.4999,lfp(:,i));
%     if ~isempty(find(removed_trials_index==q))
%     title(vars_lfp{q},'Color','red');
%     else
%         title(vars_lfp{q},'Color','black');
%     end
%     ylim([-0.18 0.18]);
% end
% 
% saveas(gcf,fullfile(savepath,['trial',num2str(i),'.png']));

end

savevar=[vars_lfp;vars_spk;vars_ang;vars_spd;vars_wf;vars_wft];

save(fullfile(trialpath,'data_car'), savevar{:});

end
