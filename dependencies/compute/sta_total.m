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
        matData = load(fullfile(sitepath,'results/sta/totallfp_sta.mat'));
        % matData2 = load(fullfile(sitepath,'results/sta/baselinelfp_sta.mat'));
        % if exist(fullfile(sitepath,'combined/p_data_detail2.mat'))
        % load(fullfile(sitepath,'combined/p_data_detail2.mat'));
        % else
        %     load(fullfile(sitepath,'combined/p_data_detail.mat'));
        % end
        varNames = fieldnames(matData);
        for k = 1:length(varNames)
            varName = varNames{k};
            specData = matData.(varName);
            %totalspec(:,num)= normalize(zscore(specData),'range',[-1,1]);
            totalspec(:,num)= specData;
            % specData2 = matData2.(varName);
            % %totalspec2(:,num)= normalize(zscore(specData2),'range',[-1,1]);
            % totalspec2(:,num)= specData2;
            num=num+1;
        end
    end
end
totalspec=mean(totalspec,2);
% baseline=mean(totalspec2,2);
disp(['combine lfp sta Finished!']);
cd(workpath);  
%%
t=-0.04:0.001:0.04;
figure;
hold on
plot(t,1000*(totalspec-mean(totalspec)) ,'k', 'LineWidth', 2,'Color','black');
    xlabel('Time (ms)');ylabel('Voltage (mV)');
plot(t,1000*(baseline-mean(baseline)) ,'k', 'LineWidth', 2,'Color',[0.8, 0.8, 0.8]);
    xlabel('Time (s)');ylabel('Voltage (mV)');
    hold off
lgd=legend('Stimulus','Baseline');
lgd.FontSize = 18;
lgd.Box = 'off';

% 移除x轴轴线，只保留刻度
ax = gca; % 获取当前坐标轴
ax.XTick=[-0.04,0,0.04];
ax.XTickLabel=["40","0","40"];
yticks([-1.5,0,1.5]);
set(gca,'FontSize',18)
%%
[delta,theta,alpha,beta,gamma]=detaild_filter(totalspec);
%% 
colors = lines(5); 
figure;
set(gcf,'unit','normalized','position',[0,0,0.3,0.3]);
plot(t,1000*(delta-mean(delta)) ,'k', 'LineWidth', 2,'Color', colors(1,:));
xlabel('Time (ms)');ylabel('Voltage (mV)');
ax = gca; % 获取当前坐标轴
ax.XTick=[-0.04,0,0.04];
ax.XTickLabel=["40","0","40"];
ylim([-0.02 0.02])
yticks([-0.02,0,0.02]);
set(gca,'FontSize',18)
title('δ');



figure;
set(gcf,'unit','normalized','position',[0,0,0.3,0.3]);
plot(t,1000*(theta-mean(theta)) ,'k', 'LineWidth', 2,'Color', colors(2,:));
xlabel('Time (ms)');ylabel('Voltage (mV)');
ax = gca; % 获取当前坐标轴
ax.XTick=[-0.04,0,0.04];
ax.XTickLabel=["40","0","40"];
ylim([-0.02 0.02])
yticks([-0.02,0,0.02]);
set(gca,'FontSize',18)
title('θ');


figure;
set(gcf,'unit','normalized','position',[0,0,0.3,0.3]);
plot(t,1000*(alpha-mean(alpha)) ,'k', 'LineWidth', 2,'Color', colors(3,:));
xlabel('Time (ms)');ylabel('Voltage (mV)');
ax = gca; % 获取当前坐标轴
ax.XTick=[-0.04,0,0.04];
ax.XTickLabel=["40","0","40"];
ylim([-0.01 0.01])
yticks([-0.01,0,0.01]);
set(gca,'FontSize',18)
title('α');

figure;
set(gcf,'unit','normalized','position',[0,0,0.3,0.3]);
plot(t,1000*(beta-mean(beta)) ,'k', 'LineWidth', 2,'Color', colors(4,:));
xlabel('Time (ms)');ylabel('Voltage (mV)');
ax = gca; % 获取当前坐标轴
ax.XTick=[-0.04,0,0.04];
ax.XTickLabel=["40","0","40"];
ylim([-1 1])
yticks([-1,0,1]);
set(gca,'FontSize',18)
title('β');

figure;
set(gcf,'unit','normalized','position',[0,0,0.3,0.3]);
plot(t,1000*(gamma-mean(gamma)) ,'k', 'LineWidth', 2,'Color', colors(5,:));
xlabel('Time (ms)');ylabel('Voltage (mV)');
ax = gca; % 获取当前坐标轴
ax.XTick=[-0.04,0,0.04];
ax.XTickLabel=["40","0","40"];
ylim([-1 1])
yticks([-1,0,1]);
set(gca,'FontSize',18)
title('γ');



