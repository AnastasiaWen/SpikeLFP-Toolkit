%%
clear;
clc;

load('G:\2022_Mn_14365\script\GC_f_tt.mat')
totalC=[];
totalphi=[];
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
        if isfile(fullfile(sitepath,'results/newGC/totalGC.mat'))
            if isfile(fullfile(sitepath,'results/channel_cohe_ch/totalCoh2.mat'))
                matData = load(fullfile(sitepath,'results/channel_cohe_ch/totalCoh2.mat'));
            else
                matData = load(fullfile(sitepath,'results/channel_cohe_ch/totalCoh.mat')); 
            end
        else
            continue
        end
        matDataC=matData.totalC;
        matDataphi=matData.totalphi;
        % load(fullfile(sitepath,'combined/p_data_detail.mat'));
        varNames = fieldnames(matDataC);
        for k = 1:length(varNames)
            varName = varNames{k};
            
            if ~isempty(matDataC.(varName))
                if ~exist([varName(1:5),'C'])
                eval([varName(1:5),'C(:,:,1)= matDataC.(varName);'])
                eval([varName(1:5),'phi(:,:,1)= matDataphi.(varName);'])
                else
                eval([varName(1:5),'C(:,:,end+1)= matDataC.(varName);'])
                eval([varName(1:5),'phi(:,:,end+1)= matDataphi.(varName);'])
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

workpath='G:\2022_Mn_14365\Test';
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
            if isfile(fullfile(sitepath,'results/channel_cohe_ch/totalCoh2.mat'))
                matData = load(fullfile(sitepath,'results/channel_cohe_ch/totalCoh2.mat'));
            else
                matData = load(fullfile(sitepath,'results/channel_cohe_ch/totalCoh.mat')); 
            end
        else
            continue
        end
        matDataC=matData.totalC;
        matDataphi=matData.totalphi;
        % load(fullfile(sitepath,'combined/p_data_detail.mat'));
        varNames = fieldnames(matDataC);
        for k = 1:length(varNames)
            varName = varNames{k};
            
            if ~isempty(matDataC.(varName))
                if ~exist([varName(1:5),'C'])
                eval([varName(1:5),'C(:,:,1)= matDataC.(varName);'])
                eval([varName(1:5),'phi(:,:,1)= matDataphi.(varName);'])
                else
                eval([varName(1:5),'C(:,:,end+1)= matDataC.(varName);'])
                eval([varName(1:5),'phi(:,:,end+1)= matDataphi.(varName);'])
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
set(gcf,'unit','normalized','position',[0,0,0.20,0.25]);
 % subplot(8,2,i);
% hold on
totalphi=eval(vars{2*i});
totalC=eval(vars{2*i-1});
% if size(totalC,3)==1
%     continue;
% end

idx=f>0&f<50;

totalC=totalC(:,idx,:);
totalphi=totalphi(:,idx,:);

totalC=squeeze(mean(totalC,2));
totalphi=squeeze(mean(totalphi,2));
totalC=totalC';
totalphi=totalphi';

% hold on
means2p=mean(totalC,1);
% sem=std(totalC,1)/sqrt(size(totalC,1));
% shadedErrorBar(tt-0.5,means2p,sem,'lineProps',{'Color','red','LineWidth',2,'LineStyle','-','Marker', 'none'});
% % meanp2s=mean(totalphi,1);
% % sem=std(totalphi,1)/sqrt(size(totalphi,1));
% % shadedErrorBar(tt-0.5,meanp2s,sem,'lineProps',{'Color','blue','LineWidth',2,'LineStyle','-','Marker', 'none'});
% hold off
% 
plot(tt-0.5, means2p,'LineWidth',2,'Color','red');
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
ylim([0  0.2])
box off
% yticks(0:0.5:1);
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
