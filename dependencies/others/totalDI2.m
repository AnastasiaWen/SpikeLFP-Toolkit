%%
clear;
clc;

nowpath=pwd;

workpath='G:\2022_Mn_14365\validData';
cd(workpath);   

weeklist = dir(workpath);
weeklist =weeklist ([weeklist .isdir]);

totalDI1=[];

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
        matData = load(fullfile(sitepath,'results\direction_tun\DIdata.mat'));
        matData=matData.Dindex;
        if exist(fullfile(sitepath,'combined/p_data_detail2.mat'))
        load(fullfile(sitepath,'combined/p_data_detail2.mat'));
        else
            load(fullfile(sitepath,'combined/p_data_detail.mat'));
        end
        varNames = fieldnames(matData);
        for k = 1:length(varNames)
            varName = varNames{k};
            specData = matData.(varName);
            totalDI1= [totalDI1,specData];
            num=num+1;
        end
    end
end



%%

nowpath=pwd;

workpath='G:\2023_Mn_14325_RightBrain\validData';
cd(workpath);   

weeklist = dir(workpath);
weeklist =weeklist ([weeklist .isdir]);

totalDI2=[];

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
        matData = load(fullfile(sitepath,'results\direction_tun\DIdata.mat'));
        matData=matData.Dindex;
        if exist(fullfile(sitepath,'combined/p_data_detail2.mat'))
        load(fullfile(sitepath,'combined/p_data_detail2.mat'));
        else
            load(fullfile(sitepath,'combined/p_data_detail.mat'));
        end
        varNames = fieldnames(matData);
        for k = 1:length(varNames)
            varName = varNames{k};
            specData = matData.(varName);
            totalDI2= [totalDI2,specData];
            num=num+1;
        end
    end
end

figure;

% 对数据进行排序
sortedData1 = sort(totalDI1);
% 计算CDF
n1 = numel(sortedData1); % 数据点的数量
cdf1 = (1:n1) / n1;

% 对数据进行排序
sortedData2 = sort(totalDI2);
% 计算CDF
n2 = numel(sortedData2); % 数据点的数量
cdf2 = (1:n2) / n2;

% 绘制CDF曲线
figure;
% set(gcf,'unit','normalized','position',[0,0,0.25,0.3]);
plot(sortedData1, cdf1, 'LineWidth', 2);
hold on
plot(sortedData2, cdf2, 'LineWidth', 2);
hold off
xlabel('DI');
ylabel('CDF');
title('Cumulative Distribution Function');
%grid on;
l=legend({'Mn1','Mn2'});
l.Box='off'
xticks([0 0.5 1]);
yticks([0 0.5 1])
set(gca,'FontSize',15)

% saveas(gcf,fullfile('G:\Global_Results\lfp\multitaper\mn2','incspec_cut.png'));
% save(fullfile('G:\Global_Results\lfp\multitaper\mn2','inc_lfp_spec_cut.mat'),"totalspec");
