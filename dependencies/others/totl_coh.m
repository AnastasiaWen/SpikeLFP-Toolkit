%%
clear;
clc;
load('G:\2022_Mn_14365\script\GC_f_tt.mat')
totals2p=[];
totalp2s=[];
num=0;
neuron={'neuron0fr','neuron45fr','neuron90fr','neuron135fr','neuron180fr'};
neurondi={'neuron0di','neuron45di','neuron90di','neuron135di','neuron180di'};
%%

nowpath=pwd;
% workpath='G:\2022_Mn_14365\validData';
workpath='G:\2023_Mn_14325_RightBrain\validData';
% workpath='G:\2022_Mn_14365\validData';
cd(workpath);   

weeklist = dir(workpath);
weeklist =weeklist ([weeklist .isdir]);



% num=1;

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
        if isfile(fullfile(sitepath,'results/channel_cohe_ch/totalCoh.mat'))
            if isfile(fullfile(sitepath,'results/channel_cohe_ch/totalCoh2.mat'))
                matData = load(fullfile(sitepath,'results/channel_cohe_ch/totalCoh2.mat'));
            else
                matData = load(fullfile(sitepath,'results/channel_cohe_ch/totalCoh.mat')); 
            end
        else
            continue
        end
        matDataC=matData.totalC;
        matDataPhi=matData.totalphi;
        % load(fullfile(sitepath,'combined/p_data_detail.mat'));
        varNames = fieldnames(matDataPhi);
        for k = 1:length(varNames)
            varName = varNames{k};
            
            if ~isempty(matDataC.(varName))
            totalC(:,:,num) = matDataC.(varName);
            totalPhi(:,:,num) = matDataPhi.(varName);
            num=num+1;
            end
        end
    end
end

dtotalC=totalC;
dtotalPhi=totalPhi;
cd(nowpath)

%%

nowpath=pwd;

% workpath='G:\2022_Mn_14365\validData';
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
        if isfile(fullfile(sitepath,'results/newGC/totallGC.mat'))
            if isfile(fullfile(sitepath,'results/newGC/totalGC.mat'))
                matData = load(fullfile(sitepath,'results/newGC/totalGC.mat'));
            else
                matData = load(fullfile(sitepath,'results/newGC/totallGC.mat')); 
            end
        else
            continue
        end
        matDatap2s=matData.totalGCp2s;
        matDataPhi=matData.totalGCs2p;
        % load(fullfile(sitepath,'combined/p_data_detail.mat'));
        varNames = fieldnames(matDataPhi);
        for k = 1:length(varNames)
            varName = varNames{k};
            
            if ~isempty(matDatap2s.(varName))
                tmp = matDatap2s.(varName);
                tmp=tmp./mean(tmp(6:12,:));
            totalp2s(:,:,num) = tmp;
            tmp = matDataPhi.(varName);
               tmp=tmp./mean(tmp(6:12,:));
            totals2p(:,:,num) = tmp;
            num=num+1;
            end
        end
    end
end

dtotalp2s=totalp2s;
dtotals2p=totals2p;
cd(nowpath)
%%
close all
tt=matData.tt;
t=tt-0.5;
f=matData.f;
%%
delta=squeeze(mean(dtotalC(:,f>1&f<4,:),2));
theta=squeeze(mean(dtotalC(:,f>4&f<8,:),2));
alpha=squeeze(mean(dtotalC(:,f>8&f<12,:),2));
beta=squeeze(mean(dtotalC(:,f>12&f<30,:),2));
gamma=squeeze(mean(dtotalC(:,f>30&f<50,:),2));
delta=delta';
theta=theta';
alpha=alpha';
beta=beta';
gamma=gamma';
deltaPhi=squeeze(mean(dtotalPhi(:,f>1&f<4,:),2));
thetaPhi=squeeze(mean(dtotalPhi(:,f>4&f<8,:),2));
alphaPhi=squeeze(mean(dtotalPhi(:,f>8&f<12,:),2));
betaPhi=squeeze(mean(dtotalPhi(:,f>12&f<30,:),2));
gammaPhi=squeeze(mean(dtotalPhi(:,f>30&f<50,:),2));
deltaPhi=deltaPhi';
thetaPhi=thetaPhi';
alphaPhi=alphaPhi';
betaPhi=betaPhi';
gammaPhi=gammaPhi';

%%
close all
% 假设 delta, theta, alpha, beta, gamma, deltaPhi, thetaPhi, alphaPhi, betaPhi, gammaPhi 和 t 已经在工作区中定义
bands = {'delta', 'theta', 'alpha', 'beta', 'gamma'};  % 各个频段名称
PhiBands = {'deltaPhi', 'thetaPhi', 'alphaPhi', 'betaPhi', 'gammaPhi'};  % 各个频段的相位名称

% 设置直方图分辨率为 0.05
binWidth = 0.05;

% 遍历每个频段
for i = 1:length(bands)
    band = eval(bands{i});  % 获取当前频段的coherence矩阵
    PhiBand = eval(PhiBands{i});  % 获取当前频段的相位矩阵
    
    maxCoherence = [];  % 存放当前频段的最大coherence值
    timePoints = [];    % 存放当前频段的对应时间点
    phaseValues = [];   % 存放当前频段的对应相位值
    
    % 遍历每一行
    for row = 1:size(band, 1)
        [maxVal, idx] = max(band(row, :));  % 找到当前行的最大值及其索引
        maxCoherence = [maxCoherence; maxVal];  % 保存最大值
        timePoints = [timePoints; t(idx)];  % 保存对应的时间点
        phaseValues = [phaseValues; PhiBand(row, idx)];  % 保存对应的相位值
    end
    
    % 绘制当前频段的直方图
    figure;
    set(gcf,'unit','normalized','position',[0,0,0.2,0.25]);

    
    histogram(maxCoherence, 'BinWidth', 0.004);
    mean(maxCoherence)
    % title([bands{i} ' Max Coherence Values=',num2str(mean(maxCoherence))]);
    xlim([0 0.2])
    % xlabel('Peak of Coherence');
    % ylabel('Frequency');
    box off;
    set(gca,'LineWidth',2)
    set(gca,'FontSize',15)
% set(gca,'FontWeight','bold')
ylim([0 20])
    figure;
    set(gcf,'unit','normalized','position',[0,0,0.2,0.3]);

    set(gca,'LineWidth',2)
    set(gca,'FontSize',15)
    histogram(timePoints, 'BinWidth', binWidth,'FaceColor',[0.1 0 0.2]);
    % title([bands{i} ' Corresponding Time Points']);
    xlabel('Time');
    % ylabel('Frequency');
    set(gca,'LineWidth',2)
    set(gca,'FontSize',15)
% set(gca,'FontWeight','bold')
xlim([-0.2 1.2])
xticks([-0.2:0.2:1.2])
xtickangle(45)
   box off;
ylim([0 20])
    % 绘制相位值的雷达图
    figure;
    phaseBins = -pi:pi/18:pi;  % 将相位分成12个部分，等间隔0到2pi
    phaseHist = histcounts(phaseValues, phaseBins);  % 计算相位值的频数
    csvwrite('my_array.csv', phaseHist/size(band, 1)); 
    phaseBinsCenter = phaseBins(1:end-1) + diff(phaseBins)/2;  % 获取每个bin的中心点
    
    % 雷达图数据准备
    raderp = [phaseBinsCenter phaseBinsCenter(1)];  % 补回起点以封闭雷达图
    rho = [phaseHist phaseHist(1)];  % 补回起点的值
    
    % 绘制雷达图
    polarplot(raderp, rho, '-o');
    title([bands{i} ' Phase Values (Radar Plot)']);
end

%%
figure;
set(gcf,'unit','normalized','position',[0,0,0.20,0.25]);
hold on
meandelta=mean(delta,1);
sem=std(delta,1)/sqrt(size(delta,1));
shadedErrorBar(tt-0.5,meandelta,sem,'lineProps',{'Color',[0.5 0.8 0.82],'LineWidth',2,'LineStyle','-','Marker', 'none'});
meantheta=mean(theta,1);
sem=std(theta,1)/sqrt(size(theta,1));
shadedErrorBar(tt-0.5,meantheta,sem,'lineProps',{'Color',[1 0.3 0.5],'LineWidth',2,'LineStyle','-','Marker', 'none'});
meanalpha=mean(alpha,1);
sem=std(alpha,1)/sqrt(size(alpha,1));
shadedErrorBar(tt-0.5,meanalpha,sem,'lineProps',{'Color',[0.8 0.8 0.2],'LineWidth',2,'LineStyle','-','Marker', 'none'});
meanbeta=mean(beta,1);
sem=std(beta,1)/sqrt(size(beta,1));
shadedErrorBar(tt-0.5,meanbeta,sem,'lineProps',{'Color',[0.5 0.3 0.2],'LineWidth',2,'LineStyle','-','Marker', 'none'});
meangamma=mean(gamma,1);
sem=std(gamma,1)/sqrt(size(gamma,1));
shadedErrorBar(tt-0.5,meangamma,sem,'lineProps',{'Color',[0.1 0.8 0.2],'LineWidth',2,'LineStyle','-','Marker', 'none'});
hold off
xlim([-0.2,1.2])
xticks([-0.2:0.2:1.2]);
yticks([0])
% yticks([0.02:0.02:0.06])
ylim([0.02 0.06])
ax = gca;        % 获取当前轴
ax.YColor = 'k'; 
ax.LineWidth=2;
ax.FontWeight="bold"
ax.FontSize=15;
ax.Box='on'

%%

%%
idx=f>4&f<8;

totalC=dtotalC(:,idx,:);
% totalp2s=dtotalp2s(:,idx,:);

totalC=squeeze(mean(totalC,2));
% totalp2s=squeeze(mean(totalp2s,2));
totalC=totalC';
% totalp2s=totalp2s';

%%
% 绘制曲线
figure;
set(gcf,'unit','normalized','position',[0,0,0.2,0.25]);

% subplot(1,2,1);
hold on
meanC=mean(totalC,1);
sem=std(totalC,1)/sqrt(size(totalC,1));
shadedErrorBar(tt-0.5,meanC,sem,'lineProps',{'Color','red','LineWidth',2,'LineStyle','-','Marker', 'none'});
% meanp2s=mean(totalp2s,1);
% sem=std(totalp2s,1)/sqrt(size(totalp2s,1));
% shadedErrorBar(tt-0.5,meanp2s,sem,'lineProps',{'Color','blue','LineWidth',2,'LineStyle','-','Marker', 'none'});
% plot(tt-0.5, means2p,'LineWidth',2,'Color','red');
% % title('Time-GC');xlabel('time (s)');ylabel('GC');
% hold on;
% plot(tt-0.5, meanp2s,'LineWidth',2,'Color','blue');
% title('Time-GC');xlabel('time (s)');ylabel('GC');
hold off;
set(gca,'LineWidth',3);
% plot(tt-0.5, sum(totals2p(:,idx),2));title('Time-GC');xlabel('time (s)');ylabel('GC');hold on;
% plot(tt-0.5, sum(totalp2s(:,idx),2));title('Time-GC');xlabel('time (s)');ylabel('GC');hold off;
% l1=legend('GCs2p','GCp2s');
% l1.Box='off';
box off
% l1.LineWidth=10;
%xlim([-0.5 1.5]);
% ytickformat('%.1e'); 
xticks(-0.2:0.2:1.2);
xtickangle(45)
xlim([-0.2,1.2])
 % ylim([0  0.002])
 % yticks([0 50 150]);
 % ylim([0  3])
 % ylim([0.0004  0.0015])
 % yticks([0 3 6]);
set(gca,'FontSize',20)
set(gca,'FontWeight','bold')
% subplot(1,2,2);
% plot(f, sum(totals2p,1));title('Freq-GC');xlabel('Frequency (Hz)');ylabel('GC');hold on;
% plot(f, sum(totalp2s,1));title( 'Freq-GC');xlabel('Frequency (Hz)');ylabel('GC');hold off;
% xlim([-0.5 1.5]);
% xticks([0,50,100]);
% yticks(ylim);
% l2=legend('GCs2p','GCp2s');
% l2.LineWidth=10;
% l2.Box='off';
% box off
% set(gca,'FontSize',15)
%%
idx=f>12&f<24;

totals2p=dtotals2p(:,idx,:);
totalp2s=dtotalp2s(:,idx,:);

totals2p=squeeze(mean(totals2p,2));
totalp2s=squeeze(mean(totalp2s,2));
totals2p=totals2p';
totalp2s=totalp2s';

%%
% 绘制曲线
figure;
set(gcf,'unit','normalized','position',[0,0,0.2,0.25]);

% subplot(1,2,1);
hold on
means2p=mean(totals2p,1);
sem=std(totals2p,1)/sqrt(size(totals2p,1));
shadedErrorBar(tt-0.5,means2p,sem,'lineProps',{'Color','red','LineWidth',2,'LineStyle','-','Marker', 'none'});
meanp2s=mean(totalp2s,1);
sem=std(totalp2s,1)/sqrt(size(totalp2s,1));
shadedErrorBar(tt-0.5,meanp2s,sem,'lineProps',{'Color','blue','LineWidth',2,'LineStyle','-','Marker', 'none'});
% plot(tt-0.5, means2p,'LineWidth',2,'Color','red');
% % title('Time-GC');xlabel('time (s)');ylabel('GC');
% hold on;
% plot(tt-0.5, meanp2s,'LineWidth',2,'Color','blue');
% title('Time-GC');xlabel('time (s)');ylabel('GC');
hold off;
set(gca,'LineWidth',3);
% plot(tt-0.5, sum(totals2p(:,idx),2));title('Time-GC');xlabel('time (s)');ylabel('GC');hold on;
% plot(tt-0.5, sum(totalp2s(:,idx),2));title('Time-GC');xlabel('time (s)');ylabel('GC');hold off;
% l1=legend('GCs2p','GCp2s');
% l1.Box='off';
box off
% l1.LineWidth=10;
%xlim([-0.5 1.5]);
% ytickformat('%.1e'); 
xticks(-0.2:0.2:1.2);
xtickangle(45)
xlim([-0.2,1.2])
 % ylim([0  0.002])
 % yticks([0 50 150]);
 % ylim([0  3])
 % yticks([0 3 6]);
set(gca,'FontSize',20)
set(gca,'FontWeight','bold')
% subplot(1,2,2);
% plot(f, sum(totals2p,1));title('Freq-GC');xlabel('Frequency (Hz)');ylabel('GC');hold on;
% plot(f, sum(totalp2s,1));title( 'Freq-GC');xlabel('Frequency (Hz)');ylabel('GC');hold off;
% xlim([-0.5 1.5]);
 ylim([0.0004  0.0015])
% xticks([0,50,100]);
% yticks(ylim);
% l2=legend('GCs2p','GCp2s');
% l2.LineWidth=10;
% l2.Box='off';
% box off
% set(gca,'FontSize',15)
%%
idx=f>8&f<12;

totals2p=dtotals2p(:,idx,:);
totalp2s=dtotalp2s(:,idx,:);

totals2p=squeeze(mean(totals2p,2));
totalp2s=squeeze(mean(totalp2s,2));
totals2p=totals2p';
totalp2s=totalp2s';

%%
% 绘制曲线
figure;
set(gcf,'unit','normalized','position',[0,0,0.2,0.25]);

% subplot(1,2,1);
hold on
means2p=mean(totals2p,1);
sem=std(totals2p,1)/sqrt(size(totals2p,1));
shadedErrorBar(tt-0.5,means2p,sem,'lineProps',{'Color','red','LineWidth',2,'LineStyle','-','Marker', 'none'});
meanp2s=mean(totalp2s,1);
sem=std(totalp2s,1)/sqrt(size(totalp2s,1));
shadedErrorBar(tt-0.5,meanp2s,sem,'lineProps',{'Color','blue','LineWidth',2,'LineStyle','-','Marker', 'none'});
% plot(tt-0.5, means2p,'LineWidth',2,'Color','red');
% % title('Time-GC');xlabel('time (s)');ylabel('GC');
% hold on;
% plot(tt-0.5, meanp2s,'LineWidth',2,'Color','blue');
% title('Time-GC');xlabel('time (s)');ylabel('GC');
hold off;
set(gca,'LineWidth',3);
% plot(tt-0.5, sum(totals2p(:,idx),2));title('Time-GC');xlabel('time (s)');ylabel('GC');hold on;
% plot(tt-0.5, sum(totalp2s(:,idx),2));title('Time-GC');xlabel('time (s)');ylabel('GC');hold off;
% l1=legend('GCs2p','GCp2s');
% l1.Box='off';
box off
% l1.LineWidth=10;
%xlim([-0.5 1.5]);
% ytickformat('%.1e'); 
xticks(-0.2:0.2:1.2);
xtickangle(45)
xlim([-0.2,1.2])
 % ylim([0  0.002])
 % yticks([0 50 150]);
 % ylim([0  3])
 % yticks([0 3 6]);
 ylim([0.0004  0.0015])
set(gca,'FontSize',20)
set(gca,'FontWeight','bold')
% subplot(1,2,2);
% plot(f, sum(totals2p,1));title('Freq-GC');xlabel('Frequency (Hz)');ylabel('GC');hold on;
% plot(f, sum(totalp2s,1));title( 'Freq-GC');xlabel('Frequency (Hz)');ylabel('GC');hold off;
% xlim([-0.5 1.5]);
% xticks([0,50,100]);
% yticks(ylim);
% l2=legend('GCs2p','GCp2s');
% l2.LineWidth=10;
% l2.Box='off';
% box off
% set(gca,'FontSize',15)
%%
idx=f>4&f<8;

totals2p=dtotals2p(:,idx,:);
totalp2s=dtotalp2s(:,idx,:);

totals2p=squeeze(mean(totals2p,2));
totalp2s=squeeze(mean(totalp2s,2));
totals2p=totals2p';
totalp2s=totalp2s';

%%
% 绘制曲线
figure;
set(gcf,'unit','normalized','position',[0,0,0.2,0.25]);

% subplot(1,2,1);
hold on
means2p=mean(totals2p,1);
sem=std(totals2p,1)/sqrt(size(totals2p,1));
shadedErrorBar(tt-0.5,means2p,sem,'lineProps',{'Color','red','LineWidth',2,'LineStyle','-','Marker', 'none'});
meanp2s=mean(totalp2s,1);
sem=std(totalp2s,1)/sqrt(size(totalp2s,1));
shadedErrorBar(tt-0.5,meanp2s,sem,'lineProps',{'Color','blue','LineWidth',2,'LineStyle','-','Marker', 'none'});
% plot(tt-0.5, means2p,'LineWidth',2,'Color','red');
% % title('Time-GC');xlabel('time (s)');ylabel('GC');
% hold on;
% plot(tt-0.5, meanp2s,'LineWidth',2,'Color','blue');
% title('Time-GC');xlabel('time (s)');ylabel('GC');
hold off;
set(gca,'LineWidth',3);
% plot(tt-0.5, sum(totals2p(:,idx),2));title('Time-GC');xlabel('time (s)');ylabel('GC');hold on;
% plot(tt-0.5, sum(totalp2s(:,idx),2));title('Time-GC');xlabel('time (s)');ylabel('GC');hold off;
% l1=legend('GCs2p','GCp2s');
% l1.Box='off';
box off
% l1.LineWidth=10;
%xlim([-0.5 1.5]);
% ytickformat('%.1e'); 
xticks(-0.2:0.2:1.2);
xtickangle(45)
xlim([-0.2,1.2])
 ylim([0.0004  0.0015])
 % ylim([0  0.002])
 % yticks([0 50 150]);
 % ylim([0  3])
 % yticks([0 3 6]);
set(gca,'FontSize',20)
set(gca,'FontWeight','bold')
% subplot(1,2,2);
% plot(f, sum(totals2p,1));title('Freq-GC');xlabel('Frequency (Hz)');ylabel('GC');hold on;
% plot(f, sum(totalp2s,1));title( 'Freq-GC');xlabel('Frequency (Hz)');ylabel('GC');hold off;
% xlim([-0.5 1.5]);
% xticks([0,50,100]);
% yticks(ylim);
% l2=legend('GCs2p','GCp2s');
% l2.LineWidth=10;
% l2.Box='off';
% box off
% set(gca,'FontSize',15)
%%
idx=f>1&f<4;

totals2p=dtotals2p(:,idx,:);
totalp2s=dtotalp2s(:,idx,:);

totals2p=squeeze(mean(totals2p,2));
totalp2s=squeeze(mean(totalp2s,2));
totals2p=totals2p';
totalp2s=totalp2s';

%%
% 绘制曲线
figure;
set(gcf,'unit','normalized','position',[0,0,0.2,0.25]);

% subplot(1,2,1);
hold on
means2p=mean(totals2p,1);
sem=std(totals2p,1)/sqrt(size(totals2p,1));
shadedErrorBar(tt-0.5,means2p,sem,'lineProps',{'Color','red','LineWidth',2,'LineStyle','-','Marker', 'none'});
meanp2s=mean(totalp2s,1);
sem=std(totalp2s,1)/sqrt(size(totalp2s,1));
shadedErrorBar(tt-0.5,meanp2s,sem,'lineProps',{'Color','blue','LineWidth',2,'LineStyle','-','Marker', 'none'});
% plot(tt-0.5, means2p,'LineWidth',2,'Color','red');
% % title('Time-GC');xlabel('time (s)');ylabel('GC');
% hold on;
% plot(tt-0.5, meanp2s,'LineWidth',2,'Color','blue');
% title('Time-GC');xlabel('time (s)');ylabel('GC');
hold off;
set(gca,'LineWidth',3);
% plot(tt-0.5, sum(totals2p(:,idx),2));title('Time-GC');xlabel('time (s)');ylabel('GC');hold on;
% plot(tt-0.5, sum(totalp2s(:,idx),2));title('Time-GC');xlabel('time (s)');ylabel('GC');hold off;
% l1=legend('GCs2p','GCp2s');
% l1.Box='off';
box off
% l1.LineWidth=10;
%xlim([-0.5 1.5]);
% ytickformat('%.1e'); 
xticks(-0.2:0.2:1.2);
xtickangle(45)
xlim([-0.2,1.2])
 ylim([0.0004  0.0015])
  % ylim([0  0.0015])
 % yticks([0 50 150]);
 % ylim([0  3])
 % yticks([0 3 6]);
set(gca,'FontSize',20)
set(gca,'FontWeight','bold')
% subplot(1,2,2);
% plot(f, sum(totals2p,1));title('Freq-GC');xlabel('Frequency (Hz)');ylabel('GC');hold on;
% plot(f, sum(totalp2s,1));title( 'Freq-GC');xlabel('Frequency (Hz)');ylabel('GC');hold off;
% xlim([-0.5 1.5]);
% xticks([0,50,100]);
% yticks(ylim);
% l2=legend('GCs2p','GCp2s');
% l2.LineWidth=10;
% l2.Box='off';
% box off
% set(gca,'FontSize',15)
