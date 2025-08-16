%%
clear;
clc;

load G:\Global_Results\lfp\multitaper\spec_f.mat
load G:\Global_Results\lfp\multitaper\spec_t.mat

nowpath=pwd;

% workpath='G:\2023_Mn_14325_RightBrain\validData';
workpath='G:\2022_Mn_14365\validData';

cd(workpath);   

weeklist = dir(workpath);
weeklist =weeklist ([weeklist .isdir]);

totalspec=[];

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
        varNames = fieldnames(matData);
        for k = 1:length(varNames)
            varName = varNames{k};
            specData = matData.(varName);
            % totalspec(:,:,num)= normalize(zscore(specData),'range',[-1,1]);
            specData= 10*log10(specData);

            delta(num,:)=mean(specData(frequency>=1&frequency<=3,:),1);
            theta(num,:)=mean(specData(frequency>=4&frequency<=8,:),1);
            alpha(num,:)=mean(specData(frequency>=9&frequency<=12,:),1);
            beta(num,:)=mean(specData(frequency>=12&frequency<=24,:),1);
            gamma(num,:)=mean(specData(frequency>=30&frequency<=50,:),1);
            num=num+1;

        end
        num=1;
        varNames = fieldnames(matData2);
        for k = 1:length(varNames)
            varName = varNames{k};
            frData=matData2.(varName);
            totalfr(:,num)=frData/max(frData);
            num=num+1;

        end
    end
end
totalfr=totalfr';

disp(['combine lfp wave Finished!']);
cd(workpath);  


%%
figure;
set(gcf,'unit','normalized','position',[0,0,0.20,0.25]);
hold on
% 设置长方形区域的坐标和大小
x = 0.2; % 左下角x坐标
y = 0.2; % 左下角y坐标
width = 0.9;  % 长方形的宽度
height = 0.6; % 长方形的高度
% 添加灰色长方形区域
rectangle('Position', [x, y, width, height], 'FaceColor', [0.95, 0.95, 0.95], 'EdgeColor', 'none');

yyaxis left
meanfr=mean(totalfr,1);
sem=std(totalfr,1)/sqrt(size(totalfr,1));
s3=shadedErrorBar(t-1,meanfr,sem,'lineProps',{'Color','black','LineWidth',2});% 绘制发放率变化曲线
% ylabel('Firing Rate (Hz)');
xlim([-0.2 1.2]);
 % ylim([0, 0.6]); 
 ylim([0, 1]); 
 % yticks([0:0.2:0.6])





yyaxis right
meandelta=mean(delta,1);
sem=std(delta,1)/sqrt(size(delta,1));
shadedErrorBar(t-1,meandelta,sem,'lineProps',{'Color',[0.5 0.8 0.82],'LineWidth',2,'LineStyle','-','Marker', 'none'});
meantheta=mean(theta,1);
sem=std(theta,1)/sqrt(size(theta,1));
shadedErrorBar(t-1,meantheta,sem,'lineProps',{'Color',[1 0.3 0.5],'LineWidth',2,'LineStyle','-','Marker', 'none'});
meanalpha=mean(alpha,1);
sem=std(alpha,1)/sqrt(size(alpha,1));
shadedErrorBar(t-1,meanalpha,sem,'lineProps',{'Color',[0.8 0.8 0.2],'LineWidth',2,'LineStyle','-','Marker', 'none'});
meanbeta=mean(beta,1);
sem=std(beta,1)/sqrt(size(beta,1));
shadedErrorBar(t-1,meanbeta,sem,'lineProps',{'Color',[0.5 0.3 0.2],'LineWidth',2,'LineStyle','-','Marker', 'none'});
meangamma=mean(gamma,1);
sem=std(gamma,1)/sqrt(size(gamma,1));
shadedErrorBar(t-1,meangamma,sem,'lineProps',{'Color',[0.1 0.8 0.2],'LineWidth',2,'LineStyle','-','Marker', 'none'});

xlim([-0.5 1.2]);
xticks([-0.2:0.2:1.2]);
yticks([-1:1:3])
ylim([-1 2])

% legend('δ','θ','α','β','γ');
yyaxis right; % 使用右侧的 y 轴
ax = gca;        % 获取当前轴
ax.YColor = 'k'; 
ax.LineWidth=2;
ax.FontWeight="bold"
ax.FontSize=15;
ax.Box='on'








 %%
% 示例数据，实际数据应使用你的Power和FiringRate
Power = meangamma; % 模拟功率数据，假设时长为2500ms
FiringRate = meanfr; % 模拟放电率数据

% 参数设置
startIdx = 36;  % 0.2 秒起始时间点
endIdx = 56;   % 0.8 秒结束时间点
timeWindow = endIdx - startIdx + 1; % 窗口长度 (21 个时间点)

lagRange = -15:1:15; % 时滞范围，对应约 ±450ms

corrValues = zeros(size(lagRange)); % 用来存储每个时滞的相关系数

% 计算在每个时滞下的相关性
for i = 1:length(lagRange)
    L = lagRange(i); % 当前时滞
    PowerWindow = Power(startIdx:endIdx); % 功率在0.2到0.8秒
    RateWindow = FiringRate(startIdx+L:endIdx+L); % 放电率在L到L+0.6秒
    % 计算Spearman相关系数
    corrValues(i) = corr(PowerWindow', RateWindow', 'Type', 'Spearman');
end
corrgamma=corrValues;


% 示例数据，实际数据应使用你的Power和FiringRate
Power = meandelta; % 模拟功率数据，假设时长为2500ms
FiringRate = meanfr; % 模拟放电率数据

% 参数设置
startIdx = 36;  % 0.2 秒起始时间点
endIdx = 56;   % 0.8 秒结束时间点
timeWindow = endIdx - startIdx + 1; % 窗口长度 (21 个时间点)

lagRange = -15:1:15; % 时滞范围，对应约 ±450ms

corrValues = zeros(size(lagRange)); % 用来存储每个时滞的相关系数

% 计算在每个时滞下的相关性
for i = 1:length(lagRange)
    L = lagRange(i); % 当前时滞
    PowerWindow = Power(startIdx:endIdx); % 功率在0.2到0.8秒
    RateWindow = FiringRate(startIdx+L:endIdx+L); % 放电率在L到L+0.6秒
    % 计算Spearman相关系数
    corrValues(i) = corr(PowerWindow', RateWindow', 'Type', 'Spearman');
end
corrdelta=corrValues;

% 示例数据，实际数据应使用你的Power和FiringRatea
Power = meantheta; % 模拟功率数据，假设时长为2500ms
FiringRate = meanfr; % 模拟放电率数据

% 参数设置
startIdx = 36;  % 0.2 秒起始时间点
endIdx = 56;   % 0.8 秒结束时间点
timeWindow = endIdx - startIdx + 1; % 窗口长度 (21 个时间点)

lagRange = -15:1:15; % 时滞范围，对应约 ±450ms

corrValues = zeros(size(lagRange)); % 用来存储每个时滞的相关系数

% 计算在每个时滞下的相关性
for i = 1:length(lagRange)
    L = lagRange(i); % 当前时滞
    PowerWindow = Power(startIdx:endIdx); % 功率在0.2到0.8秒
    RateWindow = FiringRate(startIdx+L:endIdx+L); % 放电率在L到L+0.6秒
    % 计算Spearman相关系数
    corrValues(i) = corr(PowerWindow', RateWindow', 'Type', 'Spearman');
end
corrtheta=corrValues;

% 示例数据，实际数据应使用你的Power和FiringRate
Power = meanalpha; % 模拟功率数据，假设时长为2500ms
FiringRate = meanfr; % 模拟放电率数据

% 参数设置
startIdx = 36;  % 0.2 秒起始时间点
endIdx = 56;   % 0.8 秒结束时间点
timeWindow = endIdx - startIdx + 1; % 窗口长度 (21 个时间点)

lagRange = -15:1:15; % 时滞范围，对应约 ±450ms

corrValues = zeros(size(lagRange)); % 用来存储每个时滞的相关系数

% 计算在每个时滞下的相关性
for i = 1:length(lagRange)
    L = lagRange(i); % 当前时滞
    PowerWindow = Power(startIdx:endIdx); % 功率在0.2到0.8秒
    RateWindow = FiringRate(startIdx+L:endIdx+L); % 放电率在L到L+0.6秒
    % 计算Spearman相关系数
    corrValues(i) = corr(PowerWindow', RateWindow', 'Type', 'Spearman');
end
corralpha=corrValues;

% 示例数据，实际数据应使用你的Power和FiringRate
Power = meanbeta; % 模拟功率数据，假设时长为2500ms
FiringRate = meanfr; % 模拟放电率数据

% 参数设置
startIdx = 36;  % 0.2 秒起始时间点
endIdx = 56;   % 0.8 秒结束时间点
timeWindow = endIdx - startIdx + 1; % 窗口长度 (21 个时间点)

lagRange = -15:1:15; % 时滞范围，对应约 ±450ms

corrValues = zeros(size(lagRange)); % 用来存储每个时滞的相关系数

% 计算在每个时滞下的相关性
for i = 1:length(lagRange)
    L = lagRange(i); % 当前时滞
    PowerWindow = Power(startIdx:endIdx); % 功率在0.2到0.8秒
    RateWindow = FiringRate(startIdx+L:endIdx+L); % 放电率在L到L+0.6秒
    % 计算Spearman相关系数
    corrValues(i) = corr(PowerWindow', RateWindow', 'Type', 'Spearman');
end
corrbeta=corrValues;

% 绘制交叉相关性函数
figure;
hold on
plot(lagRange * 30, corrdelta,'Color',[0.5 0.8 0.82], 'LineWidth', 2); % 将时间点换算成ms (30ms/点)
plot(lagRange * 30, corrtheta,'Color',[1 0.3 0.5], 'LineWidth', 2);
plot(lagRange * 30, corralpha,'Color',[0.8 0.8 0.2], 'LineWidth', 2);
plot(lagRange * 30, corrbeta,'Color',[0.5 0.3 0.2], 'LineWidth', 2);
plot(lagRange * 30, corrgamma,'Color',[0.1 0.8 0.2], 'LineWidth', 2);
xlabel('Time lag (ms)');
ylabel('Spearman correlation');
title('Cross-correlation between Power and Firing Rate');
grid on;

%%
%t=-0.04:0.001:0.04;
t=0:0.001:2.5-1/1000;
stim=totalspec(t>=1&t<=2);
base=totalspec(t>=0&t<=1);
figure;
hold on
plot(0:1/1000:1,1000*(stim-mean(stim)) ,'k', 'LineWidth', 1.5,'Color','black');
    xlabel('Time (s)');ylabel('Voltage (mV)');
plot(0:1/1000:1,1000*(base-mean(base)) ,'k', 'LineWidth', 1.5,'Color',[0.8, 0.8, 0.8]);
%     xlabel('Time (s)');ylabel('Voltage (mV)');
    hold off
lgd=legend('Stimulus','Baseline');
lgd.FontSize = 18;
lgd.Box = 'off';

% 移除x轴轴线，只保留刻度
ax = gca; % 获取当前坐标轴
ax.XTick=[0,0.5,1];

yticks([-4,0,4]);
set(gca,'FontSize',18)
set(gca,'FontWeight','bold')
%%
% [delta,theta,alpha,beta,gamma]=detaild_filter(stim);
%Set the parameters of the MTM.
Fs = 1000;                  %Sampling frequency
TW = 3;						%Choose time-bandwidth product of 3.
ntapers = 2*TW-1;			%Choose the # of tapers.
fpass = [0 100];
pad = 2;

params.Fs = Fs;				%... sampling frequency
params.tapers=[TW,ntapers];	%... time-bandwidth product & tapers.
params.fpass=fpass;       %... fpass
params.pad = pad;             %... padding.
params.trialave = 1;		%... trial average.

%speparam=params;
%speparam.fpass=[0 120];

%时间窗参数length steplength
movingwin=[0.2 0.02];

[spectro,t,frequency]=mtspecgramc(stim,movingwin,params);
%%
delta=mean(spectro(:,frequency>=1&frequency<=3),2);
theta=mean(spectro(:,frequency>=4&frequency<=8),2);
alpha=mean(spectro(:,frequency>=9&frequency<=12),2);
beta=mean(spectro(:,frequency>=12&frequency<=30),2);
gamma=mean(spectro(:,frequency>=30&frequency<=80),2);
%%
%% 
colors = lines(5); 
figure;
%set(gcf,'unit','normalized','position',[0,0,0.3,0.3]);
plot(t,10*log10(delta) ,'Color', colors(1,:));
xlabel('Time (ms)');ylabel('Voltage (mV)');
% % 移除x轴轴线，只保留刻度
% ax = gca; % 获取当前坐标轴
% ax.XTick=[0,0.5,1];
% 
% yticks([-4,0,4]);
% set(gca,'FontSize',18)
% set(gca,'FontWeight','bold')
% title('δ');

hold on

% figure;
% set(gcf,'unit','normalized','position',[0,0,0.3,0.3]);
plot(t,10*log10(theta) ,'k','Color', colors(2,:));
xlabel('Time (ms)');ylabel('Voltage (mV)');
% 移除x轴轴线，只保留刻度
% ax = gca; % 获取当前坐标轴
% ax.XTick=[0,0.5,1];
% 
% yticks([-4,0,4]);
% set(gca,'FontSize',18)
% set(gca,'FontWeight','bold')
% title('θ');


% figure;
% set(gcf,'unit','normalized','position',[0,0,0.3,0.3]);
plot(t,10*log10(alpha) ,'k','Color', colors(3,:));
xlabel('Time (ms)');ylabel('Voltage (mV)');
% 移除x轴轴线，只保留刻度
% ax = gca; % 获取当前坐标轴
% ax.XTick=[0,0.5,1];
% 
% yticks([-4,0,4]);
% set(gca,'FontSize',18)
% set(gca,'FontWeight','bold')
% title('α');
% 
% figure;
% set(gcf,'unit','normalized','position',[0,0,0.3,0.3]);
plot(t,10*log10(beta) ,'k','Color', colors(4,:));
xlabel('Time (ms)');ylabel('Voltage (mV)');
% 移除x轴轴线，只保留刻度
% ax = gca; % 获取当前坐标轴
% ax.XTick=[0,0.5,1];
% 
% yticks([-4,0,4]);
% set(gca,'FontSize',18)
% set(gca,'FontWeight','bold')
% title('β');
% 
% figure;
% set(gcf,'unit','normalized','position',[0,0,0.3,0.3]);
plot(t,10*log10(gamma) ,'k', 'Color', colors(5,:));
xlabel('Time (ms)');ylabel('Voltage (mV)');
% 移除x轴轴线，只保留刻度
% ax = gca; % 获取当前坐标轴
% ax.XTick=[0,0.5,1];
% 
% yticks([-4,0,4]);
% set(gca,'FontSize',18)
% set(gca,'FontWeight','bold')
% title('γ');

legend('delta','theta','alpha','beta','gamma');


%% 
colors = lines(5); 
figure;
set(gcf,'unit','normalized','position',[0,0,0.3,0.3]);
plot(0:1/1000:1,1000*(delta-mean(delta)) ,'k', 'LineWidth', 2,'Color', colors(1,:));
xlabel('Time (ms)');ylabel('Voltage (mV)');
% 移除x轴轴线，只保留刻度
ax = gca; % 获取当前坐标轴
ax.XTick=[0,0.5,1];

yticks([-4,0,4]);
set(gca,'FontSize',18)
set(gca,'FontWeight','bold')
title('δ');



figure;
set(gcf,'unit','normalized','position',[0,0,0.3,0.3]);
plot(0:1/1000:1,1000*(theta-mean(theta)) ,'k', 'LineWidth', 2,'Color', colors(2,:));
xlabel('Time (ms)');ylabel('Voltage (mV)');
% 移除x轴轴线，只保留刻度
ax = gca; % 获取当前坐标轴
ax.XTick=[0,0.5,1];

yticks([-4,0,4]);
set(gca,'FontSize',18)
set(gca,'FontWeight','bold')
title('θ');


figure;
set(gcf,'unit','normalized','position',[0,0,0.3,0.3]);
plot(0:1/1000:1,1000*(alpha-mean(alpha)) ,'k', 'LineWidth', 2,'Color', colors(3,:));
xlabel('Time (ms)');ylabel('Voltage (mV)');
% 移除x轴轴线，只保留刻度
ax = gca; % 获取当前坐标轴
ax.XTick=[0,0.5,1];

yticks([-4,0,4]);
set(gca,'FontSize',18)
set(gca,'FontWeight','bold')
title('α');

figure;
set(gcf,'unit','normalized','position',[0,0,0.3,0.3]);
plot(0:1/1000:1,1000*(beta-mean(beta)) ,'k', 'LineWidth', 2,'Color', colors(4,:));
xlabel('Time (ms)');ylabel('Voltage (mV)');
% 移除x轴轴线，只保留刻度
ax = gca; % 获取当前坐标轴
ax.XTick=[0,0.5,1];

yticks([-4,0,4]);
set(gca,'FontSize',18)
set(gca,'FontWeight','bold')
title('β');

figure;
set(gcf,'unit','normalized','position',[0,0,0.3,0.3]);
plot(0:1/1000:1,1000*(gamma-mean(gamma)) ,'k', 'LineWidth', 2,'Color', colors(5,:));
xlabel('Time (ms)');ylabel('Voltage (mV)');
% 移除x轴轴线，只保留刻度
ax = gca; % 获取当前坐标轴
ax.XTick=[0,0.5,1];

yticks([-4,0,4]);
set(gca,'FontSize',18)
set(gca,'FontWeight','bold')
title('γ');


%%
%%
delta=totalspec(frequency>=1&frequency<=3,:);
theta=totalspec(frequency>=4&frequency<=8,:);
alpha=totalspec(frequency>=9&frequency<=12,:);
beta=totalspec(frequency>=12&frequency<=30,:);
gamma=totalspec(frequency>=30&frequency<=50,:);
%%

meandelta=mean(delta,1);
sem=std(delta,1)/sqrt(size(tmpspec,1));
shadedErrorBar(t-1,meandelta,sem,'lineProps',{'Color',[0.5 0.8 0.82]});

%% 
colors = lines(5); 
figure;
%set(gcf,'unit','normalized','position',[0,0,0.3,0.3]);
plot(t,10*log10(delta) ,'Color', colors(1,:));
xlabel('Time (ms)');ylabel('Voltage (mV)');
% % 移除x轴轴线，只保留刻度
% ax = gca; % 获取当前坐标轴
% ax.XTick=[0,0.5,1];
% 
% yticks([-4,0,4]);
% set(gca,'FontSize',18)
% set(gca,'FontWeight','bold')
% title('δ');

hold on

% figure;
% set(gcf,'unit','normalized','position',[0,0,0.3,0.3]);
plot(t,10*log10(theta) ,'k','Color', colors(2,:));
xlabel('Time (ms)');ylabel('Voltage (mV)');
% 移除x轴轴线，只保留刻度
% ax = gca; % 获取当前坐标轴
% ax.XTick=[0,0.5,1];
% 
% yticks([-4,0,4]);
% set(gca,'FontSize',18)
% set(gca,'FontWeight','bold')
% title('θ');


% figure;
% set(gcf,'unit','normalized','position',[0,0,0.3,0.3]);
plot(t,10*log10(alpha) ,'k','Color', colors(3,:));
xlabel('Time (ms)');ylabel('Voltage (mV)');
% 移除x轴轴线，只保留刻度
% ax = gca; % 获取当前坐标轴
% ax.XTick=[0,0.5,1];
% 
% yticks([-4,0,4]);
% set(gca,'FontSize',18)
% set(gca,'FontWeight','bold')
% title('α');
% 
% figure;
% set(gcf,'unit','normalized','position',[0,0,0.3,0.3]);
plot(t,10*log10(beta) ,'k','Color', colors(4,:));
xlabel('Time (ms)');ylabel('Voltage (mV)');
% 移除x轴轴线，只保留刻度
% ax = gca; % 获取当前坐标轴
% ax.XTick=[0,0.5,1];
% 
% yticks([-4,0,4]);
% set(gca,'FontSize',18)
% set(gca,'FontWeight','bold')
% title('β');
% 
% figure;
% set(gcf,'unit','normalized','position',[0,0,0.3,0.3]);
plot(t,10*log10(gamma) ,'k', 'Color', colors(5,:));
xlabel('Time (ms)');ylabel('Voltage (mV)');
% 移除x轴轴线，只保留刻度
% ax = gca; % 获取当前坐标轴
% ax.XTick=[0,0.5,1];
% 
% yticks([-4,0,4]);
% set(gca,'FontSize',18)
% set(gca,'FontWeight','bold')
% title('γ');

legend('delta','theta','alpha','beta','gamma');


