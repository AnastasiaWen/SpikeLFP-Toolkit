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
        matDatas2p=matData.totalGCs2p;
        % load(fullfile(sitepath,'combined/p_data_detail.mat'));
        varNames = fieldnames(matDatas2p);
        for k = 1:length(varNames)
            varName = varNames{k};
            
            if ~isempty(matDatap2s.(varName))
                if ~exist([varName(1:5),'p2s'])
                eval([varName(1:5),'p2s(:,:,1)= matDatap2s.(varName);'])
                eval([varName(1:5),'s2p(:,:,1)= matDatas2p.(varName);'])
                else
                eval([varName(1:5),'p2s(:,:,end+1)= matDatap2s.(varName);'])
                eval([varName(1:5),'s2p(:,:,end+1)= matDatas2p.(varName);'])
                end
            end
        end
    end
end
% 
% dtotalp2s=totalp2s;
% dtotals2p=totals2p;
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
        matDatas2p=matData.totalGCs2p;
        % load(fullfile(sitepath,'combined/p_data_detail.mat'));
        varNames = fieldnames(matDatas2p);
        for k = 1:length(varNames)
            varName = varNames{k};
            
            if ~isempty(matDatap2s.(varName))
                if ~exist([varName(1:5),'p2s'])
                                    tmp = matDatap2s.(varName);
                tmp=tmp./mean(tmp(6:12,:));

                eval([varName(1:5),'p2s(:,:,1)= tmp;'])
                            tmp = matDatas2p.(varName);
               tmp=tmp./mean(tmp(6:12,:));
                eval([varName(1:5),'s2p(:,:,1)= tmp;'])
                else
                                    tmp = matDatap2s.(varName);
                tmp=tmp./mean(tmp(6:12,:));

                eval([varName(1:5),'p2s(:,:,end+1)= tmp;'])
                            tmp = matDatas2p.(varName);
               tmp=tmp./mean(tmp(6:12,:));
                eval([varName(1:5),'s2p(:,:,end+1)= tmp;'])
                end
            end
        end
    end
end
% dtotalp2s=totalp2s;
% dtotals2p=totals2p;
cd(nowpath)


%%
close all
probeDic=[3,7,11,15,14,10,6,2,1,5,9,13,16,12,8,4];
vars=who('-regexp','^LFP');
% 绘制曲线
% figure;
for i=1:numel(vars)/2
    chname=vars{2*i};
    chname=str2num(chname(4:5));
    figure;
set(gcf,'unit','normalized','position',[0,0,0.3,0.25]);
 % subplot(8,2,i);
% hold on
totals2p=eval(vars{2*i});
totalp2s=eval(vars{2*i-1});

idx=f>0&f<50;

totals2p=totals2p(:,idx,:);
totalp2s=totalp2s(:,idx,:);

totals2p=squeeze(mean(totals2p,2));
totalp2s=squeeze(mean(totalp2s,2));
totals2p=totals2p';
totalp2s=totalp2s';

hold on
means2p=mean(totals2p,1);
sem=std(totals2p,1)/sqrt(size(totals2p,1));
shadedErrorBar(tt-0.5,means2p,sem,'lineProps',{'Color','red','LineWidth',2,'LineStyle','-','Marker', 'none'});
meanp2s=mean(totalp2s,1);
sem=std(totalp2s,1)/sqrt(size(totalp2s,1));
shadedErrorBar(tt-0.5,meanp2s,sem,'lineProps',{'Color','blue','LineWidth',2,'LineStyle','-','Marker', 'none'});
hold off
% 
% plot(tt-0.5, means2p,'LineWidth',2,'Color','red');
% % title('Time-GC');xlabel('time (s)');ylabel('GC');
% hold on;
% plot(tt-0.5, meanp2s,'LineWidth',2,'Color','blue');
% % title('Time-GC');xlabel('time (s)');ylabel('GC');
% hold off;
% plot(tt-0.5, sum(totals2p(:,idx),2));title('Time-GC');xlabel('time (s)');ylabel('GC');hold on;
% plot(tt-0.5, sum(totalp2s(:,idx),2));title('Time-GC');xlabel('time (s)');ylabel('GC');hold off;
% l1=legend('GCs2p','GCp2s');
% l1.Box='off';
% box off
% l1.LineWidth=10;
% xlim([-0.2 1.2]);
xticks(-0.2:0.2:1.2);
xlim([-0.2,1.2])
ylim([0  0.004])
% yticks([0 0.25]);
set(gca,'FontSize',15)
set(gca,'FontWeight','bold')
set(gca,'LineWidth',3);
title(['LFP',num2str(probeDic(chname))]);
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
end
