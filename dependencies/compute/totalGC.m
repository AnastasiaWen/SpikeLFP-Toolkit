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
% workpath='G:\2022_Mn_14365\validData';
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
            totalp2s(:,:,num) = matDatap2s.(varName);
            totals2p(:,:,num) = matDatas2p.(varName);
            num=num+1;
            end
        end
    end
end

dtotalp2s=totalp2s;
dtotals2p=totals2p;
cd(nowpath)

%%

nowpath=pwd;

workpath='G:\2022_Mn_14365\validData';
% workpath='G:\2023_Mn_14325_RightBrain\validData';
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
                tmp = matDatap2s.(varName);
                tmp=tmp./mean(tmp(2:9,:));
            totalp2s(:,:,num) = tmp;
            tmp = matDatas2p.(varName);
               tmp=tmp./mean(tmp(2:9,:));
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
f=matData.f;
%%
idx=f>30&f<50;

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
 % ylim([0.0004  0.0015])
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
 % ylim([0.0004  0.0015])
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
 % ylim([0.0004  0.0015])
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
 % ylim([0.0004  0.0015])
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
