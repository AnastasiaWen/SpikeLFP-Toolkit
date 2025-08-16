%画图的 total lfp average power 
%%
clear;
clc;

nowpath=pwd;

workpath='G:\2023_Mn_14325_RightBrain\validData';
cd(workpath);   

weeklist = dir(workpath);
weeklist =weeklist ([weeklist .isdir]);



num=1;

for j=3:numel(weeklist)
    weekname=weeklist(j).name;
    weekpath=fullfile(workpath,weekname);
    disp(['Pre-processing...',weekname]);
    cd(weekpath);   

    sitelist = dir(weekpath);
    sitelist =sitelist ([sitelist .isdir]);

    for i=3:numel(sitelist)
        sitename=sitelist(i).name;
        sitepath=fullfile(weekpath,sitename);
        disp(['Pre-processing...',sitename]);
        matData = load(fullfile(sitepath,'results/spectro/multitaper/new/totallfp_spec.mat'));
        if exist(fullfile(sitepath,'combined/p_data_detail2.mat'))
        load(fullfile(sitepath,'combined/p_data_detail2.mat'));
        else
            load(fullfile(sitepath,'combined/p_data_detail.mat'));
        end
        varNames = fieldnames(matData);
        for k = 1:length(varNames)
            varName = varNames{k};
            specData = matData.(varName);
            totalspec(:,:,num)= normalize(zscore(specData),'range',[-1,1]);
            num=num+1;
        end
    end
end
totalspec=mean(totalspec,3);
disp(['combine lfp spec Finished!']);
cd(workpath);   

spec23in=totalspec;
%%
clearvars -except totalspec workpath

load G:\Global_Results\lfp\multitaper\spec_f.mat
load G:\Global_Results\lfp\multitaper\spec_t.mat

%%
figure;
indx=frequency>=1&frequency<=3;
tmpspec=spec23in(indx,:);
meandelta=mean(tmpspec,1);
deltain=meandelta;
sem=std(tmpspec,1)/sqrt(size(tmpspec,1));
shadedErrorBar(t-1,meandelta,sem,'lineProps',{'Color',[1 0.3 0.5]});
hold on
tmpspec=spec23de(indx,:);
meandelta=mean(tmpspec,1);
sem=std(tmpspec,1)/sqrt(size(tmpspec,1));
shadedErrorBar(t-1,meandelta,sem,'lineProps',{'Color',[0.5 0.8 0.82]});
deltade=meandelta;
tmpspec=spec23total(indx,:);
meandelta=mean(tmpspec,1);
deltat=meandelta;
sem=std(tmpspec,1)/sqrt(size(tmpspec,1));
shadedErrorBar(t-1,meandelta,sem,'lineProps','black');
%legend('Inc','Dec','Total');
xlim([-0.5,1.2]);
xticks([0,1]);
ylim([-1 .5])
set(gca,'FontSize',15)
%%
alpha=mean(alpha,1);
delta=mean(delta,1);
theta=mean(theta,1);
beta=mean(beta,1);
gamma=mean(gamma,1);
%%
clc
indx1=t>(-0.5+1)&t<(-0.1+1);
indx=t>(0.2+1)&t<(0.4+1);
% indx=t>(0.6+1)&t<(1+1);
[h,p]=ttest2(delta(indx1),delta(indx))
[h,p]=ttest2(theta(indx1),theta(indx))
[h,p]=ttest2(alpha(indx1),alpha(indx))
[h,p]=ttest2(beta(indx1),beta(indx))
[h,p]=ttest2(gamma(indx1),gamma(indx))
%%
indx=13:26;
time1_data=[deltain(indx)',deltade(indx)',deltat(indx)'];
% indx=t>(0.2+1)&t<(0.6+1);
indx=t>(0.6+1)&t<(1+1);
time2_data=[deltain(indx)',deltade(indx)',deltat(indx)'];
[h,p]=ttest2(deltain(13:26),deltain(indx))
[h,p]=ttest2(deltade(13:26),deltade(indx))
[h,p]=ttest2(deltat(13:26),deltat(indx))
%%
% 假设 time1_data 和 time2_data 分别为三个条件下两个时间段的数据矩阵 (74x3)

% 将数据合并为一个表格
data = [time1_data(:,1), time2_data(:,1), time1_data(:,2), time2_data(:,2), time1_data(:,3), time2_data(:,3)];
dataTable = array2table(data, 'VariableNames', {'Time1_Cond1', 'Time2_Cond1', 'Time1_Cond2', 'Time2_Cond2', 'Time1_Cond3', 'Time2_Cond3'});

% 定义重复测量模型
rm = fitrm(dataTable, 'Time1_Cond1-Time2_Cond3 ~ 1', 'WithinDesign', table([1 2 1 2 1 2]', [1 1 2 2 3 3]', 'VariableNames', {'Time', 'Condition'}));

% 执行重复测量 ANOVA
ranovaResults = ranova(rm, 'WithinModel', 'Time*Condition');

% 显示结果
disp(ranovaResults);

%%


figure;
set(gcf,'unit','normalized','position',[0,0,0.25,0.3]);
imagesc(t-1,frequency,totalspec);
%imagesc(t-1,frequency,normalize(zscore(totalspec),'range',[-1 1]));
%freq = 2.^(round(log2(min(frequency))):round(log2(max(frequency))));
AX=gca;
% set(gca,"YScale","log")
%set(gca,'LineWidth',1);
%AX.YTickLabelMode = "auto";
ylim([0 100]);
AX.YTick = [0,50,100];
hold on
line([0.2,0.2],ylim,'Color','red','LineWidth',2)
line([1,1],ylim,'Color','blue','LineWidth',2)
hold off
xlim([-0.2 1.2]);
h = colorbar;
%h.Label.String = 'Z-Score';
AX.FontSize=8.5;
AX.FontWeight='bold';
%xlabel('Time/s');
xticks(-0.5:0.5:1.5);
%ylabel('Frequency/Hz');
set(AX, 'YDir', 'normal');
%set(AX,'TickDir','none');
set(AX,'FontSize',20)
set(AX,'YMinorTick',false)
caxis([-1 1])
colormap('Turbo')
%title(' Total Spectrogram for Mn14365');
title(' Mn2');
xlim([-0.2 1.2]);
box off     % 取消边框
ax1 = axes('Position',get(gca,'Position'),'XAxisLocation','top',...
    'YAxisLocation','right','Color','none','XColor','k','YColor','k');  % 设置坐标区
set(ax1,'XTick', [],'YTick', []);   % 去掉xy轴刻度
colormap('Turbo')
% saveas(gcf,fullfile('G:\Global_Results\lfp\multitaper\mn2','incspec_cut.png'));
% save(fullfile('G:\Global_Results\lfp\multitaper\mn2','inc_lfp_spec_cut.mat'),"totalspec");

%%
% clear;
% clc;

nowpath=pwd;

workpath='G:\2023_Mn_14325_RightBrain\validData';
cd(workpath);   

weeklist = dir(workpath);
weeklist =weeklist ([weeklist .isdir]);



num=1;

for j=3:numel(weeklist)
    weekname=weeklist(j).name;
    weekpath=fullfile(workpath,weekname);
    disp(['Pre-processing...',weekname]);
    cd(weekpath);   

    sitelist = dir(weekpath);
    sitelist =sitelist ([sitelist .isdir]);

    for i=3:numel(sitelist)
        sitename=sitelist(i).name;
        sitepath=fullfile(weekpath,sitename);
        disp(['Pre-processing...',sitename]);
        matData = load(fullfile(sitepath,'results/spectro/multitaper/cut2/totallfp_spec.mat'));
        if exist(fullfile(sitepath,'combined/p_data_detail2.mat'))
        load(fullfile(sitepath,'combined/p_data_detail2.mat'));
        else
            load(fullfile(sitepath,'combined/p_data_detail.mat'));
        end
        varNames = decreaseName;
        for k = 1:length(varNames)
            varName = varNames{k};
            specData = matData.(['LFP',varName(end-2:end)]);
            totalspec(:,:,num)= normalize(zscore(specData),'range',[-1,1]);
            num=num+1;
        end
    end
end
totalspec=mean(totalspec,3);
disp(['combine lfp spec Finished!']);
cd(workpath);   
spec23de=totalspec;
%%
clearvars -except totalspec workpath

load G:\Global_Results\lfp\multitaper\spec_f.mat
load G:\Global_Results\lfp\multitaper\spec_t.mat

figure;
set(gcf,'unit','normalized','position',[0,0,0.25,0.3]);
imagesc(t-1,frequency,totalspec);
%imagesc(t-1,frequency,normalize(zscore(totalspec),'range',[-1 1]));
%freq = 2.^(round(log2(min(frequency))):round(log2(max(frequency))));
AX=gca;
% set(gca,"YScale","log")
%set(gca,'LineWidth',1);
%AX.YTickLabelMode = "auto";
ylim([0 100]);
AX.YTick = [0,50,100];
hold on
line([0.2,0.2],ylim,'Color','red','LineWidth',2)
line([1,1],ylim,'Color','blue','LineWidth',2)
hold off
xlim([-0.2 1.2]);
h = colorbar;
%h.Label.String = 'Z-Score';
AX.FontSize=8.5;
AX.FontWeight='bold';
%xlabel('Time/s');
xticks(-0.5:0.5:1.5);
%ylabel('Frequency/Hz');
set(AX, 'YDir', 'normal');
%set(AX,'TickDir','none');
set(AX,'FontSize',20)
set(AX,'YMinorTick',false)
caxis([-1 1])
colormap('Turbo')
%title(' Total Spectrogram for Mn14365');
title(' Mn2');
box off     % 取消边框
ax1 = axes('Position',get(gca,'Position'),'XAxisLocation','top',...
    'YAxisLocation','right','Color','none','XColor','k','YColor','k');  % 设置坐标区
set(ax1,'XTick', [],'YTick', []);   % 去掉xy轴刻度
colormap('Turbo')
% saveas(gcf,fullfile('G:\Global_Results\lfp\multitaper\mn2','decspec_cut.png'));
% save(fullfile('G:\Global_Results\lfp\multitaper\mn2','dec_lfp_spec_cut.mat'),"totalspec");
%%
clear;
clc;

nowpath=pwd;

workpath='G:\2022_Mn_14365\validData';
cd(workpath);   

weeklist = dir(workpath);
weeklist =weeklist ([weeklist .isdir]);



num=1;

for j=3:numel(weeklist)
    weekname=weeklist(j).name;
    weekpath=fullfile(workpath,weekname);
    disp(['Pre-processing...',weekname]);
    cd(weekpath);   

    sitelist = dir(weekpath);
    sitelist =sitelist ([sitelist .isdir]);

    for i=3:numel(sitelist)
        sitename=sitelist(i).name;
        sitepath=fullfile(weekpath,sitename);
        disp(['Pre-processing...',sitename]);
        matData = load(fullfile(sitepath,'results/spectro/multitaper/cut2/totallfp_spec.mat'));
        if exist(fullfile(sitepath,'combined/p_data_detail2.mat'))
        load(fullfile(sitepath,'combined/p_data_detail2.mat'));
        else
            load(fullfile(sitepath,'combined/p_data_detail.mat'));
        end
        varNames = decreaseName;
        for k = 1:length(varNames)
            varName = varNames{k};
            specData = matData.(['LFP',varName(end-2:end)]);
            totalspec(:,:,num)= normalize(zscore(specData),'range',[-1,1]);
            num=num+1;  
        end
    end
end
totalspec=mean(totalspec,3);
disp(['combine lfp spec Finished!']);
spec22de=totalspec;
%%
cd(workpath);   

clearvars -except totalspec workpath
load G:\Global_Results\lfp\multitaper\spec_f.mat
load G:\Global_Results\lfp\multitaper\spec_t.mat

figure;
set(gcf,'unit','normalized','position',[0,0,0.25,0.3]);
imagesc(t-1,frequency,totalspec);
%imagesc(t-1,frequency,normalize(zscore(totalspec),'range',[-1 1]));
%freq = 2.^(round(log2(min(frequency))):round(log2(max(frequency))));
AX=gca;
% set(gca,"YScale","log")
%set(gca,'LineWidth',1);
%AX.YTickLabelMode = "auto";
ylim([0 100]);
AX.YTick = [0,50,100];
hold on
line([0.2,0.2],ylim,'Color','red','LineWidth',2)
line([1,1],ylim,'Color','blue','LineWidth',2)
hold off
xlim([-0.2 1.2]);
h = colorbar;
%h.Label.String = 'Z-Score';
AX.FontSize=8.5;
AX.FontWeight='bold';
%xlabel('Time/s');
xticks(-0.5:0.5:1.5);
%ylabel('Frequency/Hz');
set(AX, 'YDir', 'normal');
%set(AX,'TickDir','none');
set(AX,'FontSize',20)
set(AX,'YMinorTick',false)
caxis([-1 1])
colormap('Turbo')
%title(' Total Spectrogram for Mn14365');
title(' Mn1');
colormap('Turbo')
box off     % 取消边框
ax1 = axes('Position',get(gca,'Position'),'XAxisLocation','top',...
    'YAxisLocation','right','Color','none','XColor','k','YColor','k');  % 设置坐标区
set(ax1,'XTick', [],'YTick', []);   % 去掉xy轴刻度

% saveas(gcf,fullfile('G:\Global_Results\lfp\multitaper\mn1','decspec_cut.png'));
% save(fullfile('G:\Global_Results\lfp\multitaper\mn1','dec_lfp_spec_cut.mat'),"totalspec");
%%
% clear;
% clc;

nowpath=pwd;

workpath='G:\2022_Mn_14365\validData';
cd(workpath);   

weeklist = dir(workpath);
weeklist =weeklist ([weeklist .isdir]);



num=1;

for j=3:numel(weeklist)
    weekname=weeklist(j).name;
    weekpath=fullfile(workpath,weekname);
    disp(['Pre-processing...',weekname]);
    cd(weekpath);   

    sitelist = dir(weekpath);
    sitelist =sitelist ([sitelist .isdir]);

    for i=3:numel(sitelist)
        sitename=sitelist(i).name;
        sitepath=fullfile(weekpath,sitename);
        disp(['Pre-processing...',sitename]);
        matData = load(fullfile(sitepath,'results/spectro/multitaper/cut2/totallfp_spec.mat'));
        if exist(fullfile(sitepath,'combined/p_data_detail2.mat'))
        load(fullfile(sitepath,'combined/p_data_detail2.mat'));
        else
            load(fullfile(sitepath,'combined/p_data_detail.mat'));
        end
        varNames = increaseName;
        for k = 1:length(varNames)
            varName = varNames{k};
            specData = matData.(['LFP',varName(end-2:end)]);
            % totalspec(:,:,num)= normalize(zscore(specData),'range',[-1,1]);
            totalspec(:,:,num)= 10*log10(specData);
            num=num+1;  
        end
    end
end
totalspec=mean(totalspec,3);
disp(['combine lfp spec Finished!']);
spec22in=totalspec;
%%
cd(workpath);   

clearvars -except totalspec workpath
load G:\Global_Results\lfp\multitaper\spec_f.mat
load G:\Global_Results\lfp\multitaper\spec_t.mat


figure;
set(gcf,'unit','normalized','position',[0,0,0.25,0.3]);

imagesc(t-1,frequency,totalspec);
%imagesc(t-1,frequency,normalize(zscore(totalspec),'range',[-1 1]));
%freq = 2.^(round(log2(min(frequency))):round(log2(max(frequency))));
AX=gca;
% set(gca,"YScale","log")
%set(gca,'LineWidth',1);
%AX.YTickLabelMode = "auto";
ylim([0 100]);
AX.YTick = [0,50,100];
hold on
line([0.2,0.2],ylim,'Color','red','LineWidth',2)
line([1,1],ylim,'Color','blue','LineWidth',2)
hold off
xlim([-0.2 1.2]);
h = colorbar;
%h.Label.String = 'Z-Score';
AX.FontSize=8.5;
AX.FontWeight='bold';
%xlabel('Time/s');
xticks(-0.5:0.5:1.5);
%ylabel('Frequency/Hz');
set(AX, 'YDir', 'normal');
%set(AX,'TickDir','none');
set(AX,'FontSize',20)
set(AX,'YMinorTick',false)
caxis([-1 1])
colormap('Turbo')
%title(' Total Spectrogram for Mn14365');
title(' Mn1');
box off     % 取消边框
ax1 = axes('Position',get(gca,'Position'),'XAxisLocation','top',...
    'YAxisLocation','right','Color','none','XColor','k','YColor','k');  % 设置坐标区
set(ax1,'XTick', [],'YTick', []);   % 去掉xy轴刻度

% saveas(gcf,fullfile('G:\Global_Results\lfp\multitaper\mn1','incspec_cut.png'));
% save(fullfile('G:\Global_Results\lfp\multitaper\mn1','inc_lfp_spec_cut.mat'),"totalspec");
% %%



%% 画total
clear;
clc;

nowpath=pwd;
% workpath='G:\2023_Mn_14325_RightBrain\validData'
workpath='G:\2022_Mn_14365\validData';
cd(workpath);   

weeklist = dir(workpath);
weeklist =weeklist ([weeklist .isdir]);



num=1;

for j=3:numel(weeklist)
    weekname=weeklist(j).name;
    weekpath=fullfile(workpath,weekname);
    disp(['Pre-processing...',weekname]);
    cd(weekpath);   

    sitelist = dir(weekpath);
    sitelist =sitelist ([sitelist .isdir]);

    for i=3:numel(sitelist)
        sitename=sitelist(i).name;
        sitepath=fullfile(weekpath,sitename);
        disp(['Pre-processing...',sitename]);
        matData = load(fullfile(sitepath,'results/spectro/multitaper/new/totallfp_spec.mat'));
        matData2 = load(fullfile(sitepath,'results/spk_feature/fr.mat'));
        %load(fullfile(sitepath,'combined/p_data_detail.mat'));
        varNames = fieldnames(matData);
        for k = 1:length(varNames)
            varName = varNames{k};
            specData = matData.(varName);
            % totalspec(:,:,num)= normalize(zscore(specData),'range',[-1,1]);
            totalspec(:,:,num)= 10*log10(specData);
            % tmpname=['SPKC',varName(end-1:end)];
            % frData=matData2.(tmpname);
            % totalfr(:,num)=frData;
            % num=num+1;

        end
    end
end

totalspec=mean(totalspec,3);
% totalfr=mean(totalfr,2);
disp(['combine lfp spec Finished!']);
cd(workpath);   
spec23total=totalspec;
%%
clearvars -except totalspec workpath totalfr spec23total

load G:\Global_Results\lfp\multitaper\spec_f.mat
load G:\Global_Results\lfp\multitaper\spec_t.mat
windowSize = 8; % 平滑窗口大小，可以根据需要调整
totalspec = transpose(smoothSpectrum(spec23total', windowSize));
windowSize = 5; % 平滑窗口大小，可以根据需要调整
totalspec = smoothSpectrum(totalspec, windowSize);
windowSize = 4; % 平滑窗口大小，可以根据需要调整
totalspec = transpose(smoothSpectrum(totalspec', windowSize));
totalspec=smoothSpectrumWithGaussian(totalspec, 2);
% totalspec=spec23total;
figure;
 set(gcf,'unit','normalized','position',[0,0,0.24,0.26]);
imagesc(t-1,frequency,totalspec);

%imagesc(t-1,frequency,normalize(zscore(totalspec),'range',[-1 1]));
%freq = 2.^(round(log2(min(frequency))):round(log2(max(frequency))));
AX=gca;
% set(gca,"YScale","log")
set(gca,'LineWidth',1.5);
xtickangle(45);
%AX.YTickLabelMode = "auto";
ylim([0 100]);
AX.YTick = [0,50,100];
hold on

contour(t-1,frequency,totalspec, [1.8, 1.8], 'LineColor', 'red', 'LineWidth', 2);
% line([0,0],ylim,'Color','red','LineWidth',2)
% line([1,1],ylim,'Color','blue','LineWidth',2)
hold off
xlim([-0.2 1.2]);
h = colorbar;
%h.Label.String = 'Z-Score';
AX.FontSize=8.5;
AX.FontWeight='bold';
%xlabel('Time/s');
xticks([-0.2:0.2:1.2]);
%ylabel('Frequency/Hz');
set(AX, 'YDir', 'normal');
%set(AX,'TickDir','none');
set(AX,'FontSize',15)
set(AX,'YMinorTick',false)
 caxis([-1 2])
 % colormap('Turbo')
%title(' Total Spectrogram for Mn14365');
% title(' Mn2');
% box off     % 取消边框
% ax1 = axes('Position',get(gca,'Position'),'XAxisLocation','top',...
%     'YAxisLocation','right','Color','none','XColor','k','YColor','k');  % 设置坐标区
% set(ax1,'XTick', [],'YTick', []);   % 去掉xy轴刻度

% yyaxis left; % 使用右侧的 y 轴
% plot(-1:1/50:1.5,totalfr,'Color', 'black', 'LineWidth', 2); % 绘制发放率变化曲线
% ylabel('Firing Rate (Hz)');
% ylim([1, 2.5]); 
% colormap('Turbo')
% saveas(gcf,fullfile('G:\Global_Results\lfp\multitaper\mn2','total_cut.png'));
% save(fullfile('G:\Global_Results\lfp\multitaper\mn2','total_spec_cut.mat'),"totalspec");
%%
% clear;
% clc;

nowpath=pwd;

workpath='G:\2022_Mn_14365\validData';
cd(workpath);   

weeklist = dir(workpath);
weeklist =weeklist ([weeklist .isdir]);



num=1;

for j=3:numel(weeklist)
    weekname=weeklist(j).name;
    weekpath=fullfile(workpath,weekname);
    disp(['Pre-processing...',weekname]);
    cd(weekpath);   

    sitelist = dir(weekpath);
    sitelist =sitelist ([sitelist .isdir]);

    for i=3:numel(sitelist)
        sitename=sitelist(i).name;
        sitepath=fullfile(weekpath,sitename);
        disp(['Pre-processing...',sitename]);
        matData = load(fullfile(sitepath,'results/spectro/multitaper/cut/totallfp_spec.mat'));
        matData2 = load(fullfile(sitepath,'results/spk_feature/fr.mat'));
        %load(fullfile(sitepath,'combined/p_data_detail.mat'));
        varNames = fieldnames(matData);
        for k = 1:length(varNames)
            varName = varNames{k};
            specData = matData.(varName);
            % totalspec(:,:,num)= normalize(zscore(specData),'range',[-1,1]);
            totalspec(:,:,num)= 10*log10(specData);
            tmpname=['SPKC',varName(end-2:end)];
            frData=matData2.(tmpname);
            totalfr(:,num)=frData;
            num=num+1;

        end
    end
end

totalspec=mean(totalspec,3);
totalfr=mean(totalfr,2);
disp(['combine lfp spec Finished!']);
cd(workpath);   
spec22total=totalspec;
%%
cd(workpath);   

clearvars -except totalspec workpath totalfr
load G:\Global_Results\lfp\multitaper\spec_f.mat
load G:\Global_Results\lfp\multitaper\spec_t.mat

figure;
% set(gcf,'unit','normalized','position',[0,0,0.5,0.5]);
imagesc(t-1,frequency,totalspec);

%imagesc(t-1,frequency,normalize(zscore(totalspec),'range',[-1 1]));
%freq = 2.^(round(log2(min(frequency))):round(log2(max(frequency))));
AX=gca;
% set(gca,"YScale","log")
%set(gca,'LineWidth',1);
%AX.YTickLabelMode = "auto";
ylim([0 100]);
AX.YTick = [0,50,100];
hold on
 contour(t-1,frequency,totalspec, [1.6, 1.6], 'LineColor', 'red', 'LineWidth', 2);
% line([0,0],ylim,'Color','red','LineWidth',2)
% line([1,1],ylim,'Color','blue','LineWidth',2)
hold off
xlim([-0.2 1.2]);
h = colorbar;
%h.Label.String = 'Z-Score';
AX.FontSize=8.5;
AX.FontWeight='bold';
%xlabel('Time/s');
xticks([-0.2:0.2:1.2]);
%ylabel('Frequency/Hz');
set(AX, 'YDir', 'normal');
%set(AX,'TickDir','none');
set(AX,'FontSize',20)
set(AX,'YMinorTick',false)
  caxis([-1 2])
% colormap('Turbo')
%title(' Total Spectrogram for Mn14365');
% title(' Mn1');
% colormap('Turbo')
% box off     % 取消边框
% ax1 = axes('Position',get(gca,'Position'),'XAxisLocation','top',...
%     'YAxisLocation','right','Color','none','XColor','k','YColor','k');  % 设置坐标区
% set(ax1,'XTick', [],'YTick', []);   % 去掉xy轴刻度
% yyaxis left; % 使用右侧的 y 轴
% plot(-1:1/50:1.5,totalfr,'Color', 'black', 'LineWidth', 2); % 绘制发放率变化曲线
% ylabel('Firing Rate (Hz)');
% ylim([2.5, 6]); 
% saveas(gcf,fullfile('G:\Global_Results\lfp\multitaper\mn1','total_cut.png'));
% save(fullfile('G:\Global_Results\lfp\multitaper\mn1','total_spec_cut.mat'),"totalspec");


