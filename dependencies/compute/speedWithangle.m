function speedWithangle(sitepath)
%speed angle一起计算 细节要再看看 尤其是加载数据是哪个
close all;
conpath=fullfile(sitepath,'results\speed_angle_con\');
cd(conpath);
load results.mat %这个要变

savePath = fullfile(sitepath,'\results\speed_angle_sep');
mkdir(savePath);

%%
uniqSpeed=[5,10,20,25,40,50,80,100,160,200,240,320];
lowSpeed=uniqSpeed(1:6);% speed 5,10,20,25,40,50 per 10 times
highSpeed=uniqSpeed(7:12);% speed 80,100,160,200,240,320 per 10 times

uniqAngle=[0,45,90,135,180];% 0,45,90,135,180 per 24times

Fs=1000;
recordtime=-1:1/Fs:1.5-1/Fs; % 实际变为0-2.5s
timeindx=recordtime>=0.2 & recordtime <=1 | abs(recordtime-0.2)<eps | abs(recordtime-1)<eps;
needtime=recordtime(timeindx);
%%
vars=who;
indxre=startsWith(vars, 'result');
namesre=vars(find(indxre==1));

%%
for i=1:numel(namesre)
    tmpname=namesre{i};

    figure;
    %set(gcf,'unit','normalized','position',[0,0,0.17,0.62]);

    set(gcf,'unit','normalized','position',[0,0,1,1]);
    for j=1:12
        plotfire=[];
        plotlfp=[];
        plotspk=[];
        %speed为变量研究angle
        for k=1:5
            eval(['plotfire=[plotfire,mean(',tmpname,'(k,j).fr)];']);
            eval(['plotlfp=[plotlfp,',tmpname,'(k,j).lfp];']);
            eval(['plotspk=[plotspk,',tmpname,'(k,j).spkts];']);
        end
    
    ax=subplot(4,6,fix(j/7)*13+mod(j,7));
    plotlfp=lfp_norm(plotlfp);
    [dN,t]=binspikes(plotspk,50,[1.2 2]);

    plot(needtime,rescale(mean(plotlfp,2),-1,1));
    hold on;
    bar(t-1,rescale(sum(dN,2)),'r');
    alpha(0.1)
    hold off;
    title(['Speed ',num2str(uniqSpeed(j)),'mm/s']);
    %xlim([-0.5,1.5]);
    xlabel('Time/s');
    %xticks([0.2,0.4,0.6,0.8,1.0]);
    ylabel('Volt / Count');


    subplot(4,6,fix(j/7)*13+mod(j,7)+6);
    plot(uniqAngle,plotfire,'-o','Color','black','LineWidth',1);
    xlabel('Angle/°');
    ylabel('Firng Rate');
    xticks([0,45,90,135,180]);
    ylim([0,max(plotfire)+0.5]);
    %xticklabels({'0','45','90','135','180'});
    end
    sgtitle(['Neuron ',tmpname(end-2:end)]);
    % windowsize=3;
    % smoothfr=smoothdata(fire,'gaussian',windowsize);
    % %imagesc(smoothfr);
    % imagesc(fire);
    % xticks(1:12);
    % xticklabels({'5','10','20','25','40','50','80','100','160','200','240','320'});
    % yticks(1:5);
    % yticklabels({'0','45','90','135','180'});
    % %colormap("jet");
    % colorbar;
    % title('Speed-Angle FR');
    % xlabel('Speed mm/s');
    % ylabel('Angle');
    saveas(gcf,fullfile(savePath,['Neuron ',tmpname(end-2:end),'_Angle.png']))
end


%%
for i=1:numel(namesre)
    tmpname=namesre{i};

    figure;
    set(gcf,'unit','normalized','position',[0,0,1,0.5]);
    for j=1:5
        plotfire=[];
        plotlfp=[];
        plotspk=[];
        %speed为变量研究angle
        for k=1:12
            eval(['plotfire=[plotfire,mean(',tmpname,'(j,k).fr)];']);
            eval(['plotlfp=[plotlfp,',tmpname,'(j,k).lfp];']);
            eval(['plotspk=[plotspk,',tmpname,'(j,k).spkts];']);
        end
    
    ax=subplot(2,5,j);
    plotlfp=lfp_norm(plotlfp);
    [dN,t]=binspikes(plotspk,50,[1.2 2]);

    plot(needtime,rescale(mean(plotlfp,2),-1,1));
    hold on;
    bar(t-1,rescale(sum(dN,2)),'r');
    alpha(0.1)
    hold off;
    title(['Angle ',num2str(uniqAngle(j))]);
    %xlim([-0.5,1.5]);
    xlabel('Time/s');
    %xticks([0.2,0.4,0.6,0.8,1.0]);
    ylabel('Volt / Count');

    subplot(2,5,j+5);
    plot(uniqSpeed,plotfire,'-o','Color','black','LineWidth',1);
    xlabel('Speed mm/s');
    ylabel('Firng Rate');
    xticks([5,10,20,25,40,50,80,100,160,200,240,320]);
    xlim([5,320]);
    ylim([0,max(plotfire)+0.5]);
    %xticklabels({'0','45','90','135','180'});
    end
    sgtitle(['Neuron ',tmpname(end-2:end)]);
    % windowsize=3;
    % smoothfr=smoothdata(fire,'gaussian',windowsize);
    % %imagesc(smoothfr);
    % imagesc(fire);
    % xticks(1:12);
    % xticklabels({'5','10','20','25','40','50','80','100','160','200','240','320'});
    % yticks(1:5);
    % yticklabels({'0','45','90','135','180'});
    % %colormap("jet");
    % colorbar;
    % title('Speed-Angle FR');
    % xlabel('Speed mm/s');
    % ylabel('Angle');
    saveas(gcf,fullfile(savePath,['Neuron ',tmpname(end-2:end),'_Speed.png']))
end


