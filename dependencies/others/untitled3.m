%%
clear;
clc;

neuron0fr=[];
neuron45fr=[];
neuron90fr=[];
neuron135fr=[];
neuron180fr=[];
neuron0di=[];
neuron0oi=[];
neuron45di=[];
neuron45oi=[];
neuron90di=[];
neuron90oi=[];
neuron135di=[];
neuron135oi=[];
neuron180di=[];
neuron180oi=[];
neuron={'neuron0fr','neuron45fr','neuron90fr','neuron135fr','neuron180fr'};
neurondi={'neuron0di','neuron45di','neuron90di','neuron135di','neuron180di'};
neuronoi={'neuron0oi','neuron45oi','neuron90oi','neuron135oi','neuron180oi'};
%%

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
        matData = load(fullfile(sitepath,'results/direction_tun/totaldata.mat'));
        matDatadi = load(fullfile(sitepath,'results/direction_tun/Oidata.mat'));
        matDataoi = load(fullfile(sitepath,'results/orientation4_tun/Oidata.mat'));
        matData=matData.Diresults;
        matDatadi=matDatadi.Dindex;
        matDataoi=matDataoi.Dindex;
        load(fullfile(sitepath,'combined/p_data_detail.mat'));
        varNames = fieldnames(matData);
        for k = 1:length(varNames)
            varName = varNames{k};
            totalfr = matData.(varName);
            [maxv, maxindx] = max(totalfr);
            if maxindx==1
                neuron0fr=[neuron0fr;totalfr];
                neuron0di=[neuron0di;matDatadi.(varName)];
                neuron0oi=[neuron0oi;matDataoi.(varName)];
            elseif maxindx==2
                neuron45fr=[neuron45fr;totalfr];
                neuron45di=[neuron45di;matDatadi.(varName)];
                neuron45oi=[neuron45oi;matDataoi.(varName)];
            elseif maxindx==3
                neuron90fr=[neuron90fr;totalfr];
                neuron90di=[neuron90di;matDatadi.(varName)];
                neuron90oi=[neuron90oi;matDataoi.(varName)];
            elseif maxindx==4
                neuron135fr=[neuron135fr;totalfr];
                neuron135di=[neuron135di;matDatadi.(varName)];
                neuron135oi=[neuron135oi;matDataoi.(varName)];
            elseif maxindx==5
                neuron180fr=[neuron180fr;totalfr];
                neuron180di=[neuron180di;matDatadi.(varName)];
                neuron180oi=[neuron180oi;matDataoi.(varName)];
            end

        end
    end
end


cd(nowpath)

%%

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
        matData = load(fullfile(sitepath,'results/direction_tun/totaldata.mat'));
        matDatadi = load(fullfile(sitepath,'results/direction_tun/Oidata.mat'));
        matDataoi = load(fullfile(sitepath,'results/orientation4_tun/Oidata.mat'));
        matData=matData.Diresults;
        matDatadi=matDatadi.Dindex;
        matDataoi=matDataoi.Dindex;
        load(fullfile(sitepath,'combined/p_data_detail.mat'));
        varNames = fieldnames(matData);
        for k = 1:length(varNames)
            varName = varNames{k};
            totalfr = matData.(varName);
            [maxv, maxindx] = max(totalfr);
            if maxindx==1
                neuron0fr=[neuron0fr;totalfr];
                neuron0di=[neuron0di;matDatadi.(varName)];
                neuron0oi=[neuron0oi;matDataoi.(varName)];
            elseif maxindx==2
                neuron45fr=[neuron45fr;totalfr];
                neuron45di=[neuron45di;matDatadi.(varName)];
                neuron45oi=[neuron45oi;matDataoi.(varName)];
            elseif maxindx==3
                neuron90fr=[neuron90fr;totalfr];
                neuron90di=[neuron90di;matDatadi.(varName)];
                neuron90oi=[neuron90oi;matDataoi.(varName)];
            elseif maxindx==4
                neuron135fr=[neuron135fr;totalfr];
                neuron135di=[neuron135di;matDatadi.(varName)];
                neuron135oi=[neuron135oi;matDataoi.(varName)];
            elseif maxindx==5
                neuron180fr=[neuron180fr;totalfr];
                neuron180di=[neuron180di;matDatadi.(varName)];
                neuron180oi=[neuron180oi;matDataoi.(varName)];
            end

        end
    end
end

%%
close all
uniqAngle=[0,45,90,135,180];
for i=1:5
firing_rate = eval(neuron{i});
firing_di = eval(neurondi{i});
firing_di=mean(firing_di);
firing_oi = eval(neuronoi{i});
firing_oi=mean(firing_oi);

% 计算每个时间点的平均 Firing Rate 和标准误差
mean_firing_rate = mean(firing_rate, 1);  % 计算每个时间点的平均值
stderr_firing_rate = std(firing_rate, 0, 1) / sqrt(size(firing_rate, 1));  % 计算标准误差

% 创建一个新的图形窗口
figure;
set(gcf,'unit','normalized','position',[0,0.2,0.2,0.4]);
% 绘制平均 Firing Rate 及其误差条
errorbar(uniqAngle,mean_firing_rate, stderr_firing_rate,'black', 'LineWidth', 1.5);
xticks(uniqAngle);
set(gca,'FontSize',15)
set(gca,'FontWeight','bold')
xlabel('Angle / °');
ylabel('Firing Rate / Hz');
%title([neuron{i},' Firing Rate']);
ylim([2 8])
yticks([ 2 5 8])



% 标注 trial 数量
num_trials = size(firing_rate, 1);
% 获取当前图形的坐标范围
x_limits = xlim;
y_limits = ylim;

% 在图形右上角标注 trial 数量
text_x_pos = x_limits(2) - 0.4 * (x_limits(2) - x_limits(1));
text_y_pos = y_limits(2) - 0.1 * (y_limits(2) - y_limits(1));
text(text_x_pos, text_y_pos, ['n= ', num2str(num_trials)], 'FontSize', 12, 'Color', 'r');
text(text_x_pos, text_y_pos- 0.1 * (y_limits(2) - y_limits(1)), ['DI=',num2str(firing_di)], 'FontSize', 12, 'Color', 'r');
text(text_x_pos, text_y_pos- 0.2 * (y_limits(2) - y_limits(1)), ['OI=',num2str(firing_oi)], 'FontSize', 12, 'Color', 'r');
% 美化图形
% grid on;
box off;
end

%%
close all
figure;
uniqAngle=[0,45,90,135,180];
for i=1:5
firing_di = eval(neurondi{i});
% 对数据进行排序
sortedData2 = sort(firing_di);
% 计算CDF
n2 = numel(sortedData2); % 数据点的数量
cdf2 = (1:n2) / n2;

% 绘制CDF曲线

set(gcf,'unit','normalized','position',[0,0,0.2,0.4]);
plot(sortedData2, cdf2, 'LineWidth', 2);
hold on
xlabel('DI');
ylabel('CDF');
% title('Cumulative Distribution Function');
%grid on;

end
hold off
% l=legend({'0 degree','45 degree','90 degree','135 degree','180 degree'});
% l.Box='off'
 xlim([0.4 0.6]);
xticks([0.4 0.6]);
yticks([0 0.5 1])
set(gca,'FontSize',20)
set(gca,'LineWidth',2)
set(gca,'FontWeight','bold')

%%
close all
figure;
uniqAngle=[0,45,90,135,180];
for i=1:4
firing_oi = eval(neuronoi{i});
% 对数据进行排序
sortedData2 = sort(firing_oi);
% 计算CDF
n2 = numel(sortedData2); % 数据点的数量
cdf2 = (1:n2) / n2;

% 绘制CDF曲线

set(gcf,'unit','normalized','position',[0,0,0.2,0.4]);
plot(sortedData2, cdf2, 'LineWidth', 2);
hold on
xlabel('OI');
ylabel('CDF');
% title('Cumulative Distribution Function');
%grid on;

end
hold off
% l=legend({'0 degree','45 degree','90 degree','135 degree','180 degree'});
% l.Box='off'
 xlim([0. 0.4]);
 xticks([0 0.4]);
% yticks([0 0.5 1])
set(gca,'FontSize',20)
set(gca,'LineWidth',2)
set(gca,'FontWeight','bold')


%%
close all
figure;
uniqAngle=[0,45,90,135,180];
firing_oi=[];
for i=1:5
firing_oi = [firing_oi;eval(neuronoi{i})];
end
% 对数据进行排序
sortedData2 = sort(firing_oi);
% 计算CDF
n2 = numel(sortedData2); % 数据点的数量
cdf2 = (1:n2) / n2;

% 绘制CDF曲线

set(gcf,'unit','normalized','position',[0,0,0.2,0.4]);
plot(sortedData2, cdf2, 'LineWidth', 3,'Color','b');
hold on
xlabel('OI');
ylabel('CDF');
% title('Cumulative Distribution Function');
%grid on;


% hold off
% l=legend({'0 degree','45 degree','90 degree','135 degree','180 degree'});
% l.Box='off'
%  xlim([0.4 0.6]);
% xticks([0.4 0.6]);
% yticks([0 0.5 1])
set(gca,'FontSize',20)
set(gca,'LineWidth',2)
set(gca,'FontWeight','bold')


uniqAngle=[0,45,90,135,180];
firing_di=[];
for i=1:5
firing_di = [firing_di;eval(neurondi{i})];
end
% 对数据进行排序
sortedData2 = sort(firing_di);
% 计算CDF
n2 = numel(sortedData2); % 数据点的数量
cdf2 = (1:n2) / n2;

% 绘制CDF曲线

% set(gcf,'unit','normalized','position',[0,0,0.2,0.4]);
plot(sortedData2, cdf2, 'LineWidth', 3,'Color','r');
% hold on
xlabel('OI');
ylabel('CDF');
% title('Cumulative Distribution Function');
%grid on;


hold off
% l=legend({'0 degree','45 degree','90 degree','135 degree','180 degree'});
% l.Box='off'
%  xlim([0.4 0.6]);
% xticks([0.4 0.6]);
% yticks([0 0.5 1])
xlim([0 1])
set(gca,'FontSize',20)
set(gca,'LineWidth',2)
set(gca,'FontWeight','bold')


