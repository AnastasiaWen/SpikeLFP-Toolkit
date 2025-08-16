%%
clear;
clc;

neuron0fr=[];
neuron45fr=[];
neuron90fr=[];
neuron135fr=[];
neuron180fr=[];
neuron0di=[];
neuron45di=[];
neuron90di=[];
neuron135di=[];
neuron180di=[];
neuron={'neuron0fr','neuron45fr','neuron90fr','neuron135fr','neuron180fr'};
neurondi={'neuron0di','neuron45di','neuron90di','neuron135di','neuron180di'};
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
        matDatadi = load(fullfile(sitepath,'results/direction_tun/DIdata.mat'));
        matData=matData.Diresults;
        matDatadi=matDatadi.Dindex;
        load(fullfile(sitepath,'combined/p_data_detail.mat'));
        varNames = fieldnames(matData);
        for k = 1:length(varNames)
            varName = varNames{k};
            totalfr = matData.(varName);
            [maxv, maxindx] = max(totalfr);
            if maxindx==1
                neuron0fr=[neuron0fr;totalfr];
                neuron0di=[neuron0di;matDatadi.(varName)];
            elseif maxindx==2
                neuron45fr=[neuron45fr;totalfr];
                neuron45di=[neuron45di;matDatadi.(varName)];
            elseif maxindx==3
                neuron90fr=[neuron90fr;totalfr];
                neuron90di=[neuron90di;matDatadi.(varName)];
            elseif maxindx==4
                neuron135fr=[neuron135fr;totalfr];
                neuron135di=[neuron135di;matDatadi.(varName)];
            elseif maxindx==5
                neuron180fr=[neuron180fr;totalfr];
                neuron180di=[neuron180di;matDatadi.(varName)];
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
        matDatadi = load(fullfile(sitepath,'results/direction_tun/DIdata.mat'));
        matData=matData.Diresults;
        matDatadi=matDatadi.Dindex;
        load(fullfile(sitepath,'combined/p_data_detail.mat'));
        varNames = fieldnames(matData);
        for k = 1:length(varNames)
            varName = varNames{k};
            totalfr = matData.(varName);
            [maxv, maxindx] = max(totalfr);
            if maxindx==1
                neuron0fr=[neuron0fr;totalfr];
                neuron0di=[neuron0di;matDatadi.(varName)];
            elseif maxindx==2
                neuron45fr=[neuron45fr;totalfr];
                neuron45di=[neuron45di;matDatadi.(varName)];
            elseif maxindx==3
                neuron90fr=[neuron90fr;totalfr];
                neuron90di=[neuron90di;matDatadi.(varName)];
            elseif maxindx==4
                neuron135fr=[neuron135fr;totalfr];
                neuron135di=[neuron135di;matDatadi.(varName)];
            elseif maxindx==5
                neuron180fr=[neuron180fr;totalfr];
                neuron180di=[neuron180di;matDatadi.(varName)];
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

% 计算每个时间点的平均 Firing Rate 和标准误差
mean_firing_rate = mean(firing_rate, 1);  % 计算每个时间点的平均值
stderr_firing_rate = std(firing_rate, 0, 1) / sqrt(size(firing_rate, 1));  % 计算标准误差

% 创建一个新的图形窗口
figure;
set(gcf,'unit','normalized','position',[0,0.2,0.2,0.4]);
% 绘制平均 Firing Rate 及其误差条
errorbar(uniqAngle,mean_firing_rate, stderr_firing_rate, 'b', 'LineWidth', 1.5);
xticks(uniqAngle);
set(gca,'FontSize',15)
xlabel('Angle / °');
ylabel('Firing Rate / Hz');
%title([neuron{i},' Firing Rate']);
ylim([1.5 8])
yticks([ 8])


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
figure;
% set(gcf,'unit','normalized','position',[0,0,0.25,0.3]);
plot(sortedData2, cdf2, 'LineWidth', 2);
hold on
xlabel('DI');
ylabel('CDF');
title('Cumulative Distribution Function');
%grid on;

end
hold off
l=legend({'0 degree','45 degree','90 degree','135 degree','180 degree'});
l.Box='off'
xticks([0 0.5 1]);
yticks([0 0.5 1])
set(gca,'FontSize',15)

