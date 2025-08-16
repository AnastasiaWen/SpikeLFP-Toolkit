clc;clear;
load StimSpeed.mat
load data_nocut.mat
TotalSpeed=reshape(speedRec',120,1);
sortspeed=sort(TotalSpeed);
uniqSpeed=unique(sortspeed);
lowSpeed=uniqSpeed(1:6);% speed 5,10,20,25,40,50 per 10 times
highSpeed=uniqSpeed(7:12);% speed 80,100,160,200,240,320 per 10 times



vars = who;
indxL=startsWith(vars, 'LFP');
indxS=startsWith(vars, 'SPKC');
nameL=vars(find(indxL==1));
nameS=vars(find(indxS==1));


time=1-0.2;
for i=1:length(nameS)
    eval(['spk=',nameS{i},';']);
    spkname=nameS{i};
    spkname=spkname(5:end);
    for i = 1:6
    lowFR(i).Speed = lowSpeed(i);
    highFR(i).Speed = highSpeed(i);
    lowFR(i).FR=[];
    highFR(i).FR=[];
    end

    for j=1:length(spk)
        tmpspk=spk(j).times;
        tmpspeed=TotalSpeed(j);
        if isempty(tmpspk)
            fr=0;
        else
            count=tmpspk>=0.2 &tmpspk<=1;
            count=sum(count);
            fr=count/time;
        end
        
        if tmpspeed<60
            indx=find(lowSpeed==tmpspeed);
            tmpfr=lowFR(indx).FR;
            tmpfr=[tmpfr,fr];
            lowFR(indx).FR=tmpfr;
        else
            indx=find(highSpeed==tmpspeed);
            tmpfr=highFR(indx).FR;
            tmpfr=[tmpfr,fr];
            highFR(indx).FR=tmpfr;
        end
    end
    
    highfire=[];lowfire=[];
    for i = 1:6
    % 将当前结构体的 'FR' 字段的值垂直堆叠到结果数组中
    highfire = vertcat(highfire, highFR(i).FR);
    lowfire = vertcat(lowfire, lowFR(i).FR);
    end

    frname=['speFR',spkname];
    eval([frname,'.low=lowfire;']);
    eval([frname,'.high=highfire;']);


end

%%
clc;
load StimAngle.mat
%load data_nocut.mat
TotalAngle=reshape(angleRec',120,1);
sortangle=sort(TotalAngle);
uniqAngle=unique(sortangle);% 0,45,90,135,180 per 24times


% vars = who;
% indxL=startsWith(vars, 'LFP');
% indxS=startsWith(vars, 'SPKC');
% nameL=vars(find(indxL==1));
% nameS=vars(find(indxS==1));


time=1-0.2;
for i=1:length(nameS)
    eval(['spk=',nameS{i},';']);
    spkname=nameS{i};
    spkname=spkname(5:end);
    for i = 1:5
    aFR(i).Angle = uniqAngle(i);
    aFR(i).FR=[];
    end

    for j=1:length(spk)
        tmpspk=spk(j).times;
        tmpangle=TotalAngle(j);
        if isempty(tmpspk)
            fr=0;
        else
            count=tmpspk>=0.2 &tmpspk<=1;
            count=sum(count);
            fr=count/time;
        end
        
        indx=find(uniqAngle==tmpangle);
        tmpfr=aFR(indx).FR;
        tmpfr=[tmpfr,fr];
        aFR(indx).FR=tmpfr;
    end
    
    angleFR=[];
    for i = 1:5
    % 将当前结构体的 'FR' 字段的值垂直堆叠到结果数组中
        angleFR = vertcat(angleFR, aFR(i).FR);
    end
    frname=['angFR',spkname];
    eval([frname,'=angleFR;']);

end
%%
close all;
hist = readtable('spk_hist.csv');
time=hist{:,1};

namehist = hist.Properties.VariableNames;
rindx=time>=0 & time<=1.5;
rtime=time(rindx);
binWidth=0.01;%bin
% 计算 bin 的中心
bin_centers = rtime(1:end) + binWidth / 2;

vars=who;
indxFR=startsWith(vars, 'speFR');
namesFR=vars(find(indxFR==1));
indxFR=startsWith(vars, 'angFR');
nameaFR=vars(find(indxFR==1));
for i=1:numel(namesFR)
    column_data = hist{:, i+1};
    rdata=column_data(rindx);

    eval(['lowfire=',namesFR{i},'.low;']);
    eval(['highfire=',namesFR{i},'.high;']);
    eval(['anglefire=',nameaFR{i},';']);
    figure;
    set(gcf,'unit','normalized','position',[0,0,0.15,0.84]);
    ax=subplot(4,1,1);
    bar(ax,bin_centers, rdata, 'hist');
    line(ax,[0.2,0.2],ylim,'Color','red','LineWidth',2);
    line(ax,[1,1],ylim,'Color','Blue','LineWidth',2);
    xlim([0 1.5]);
    psthname=namesFR{i};
    psthname=psthname(6:end);
    title(['PSTH ',psthname]);
    xlabel('time/s');
    ylabel('count');
    subplot(4,1,2);
    errorbar(lowSpeed,mean(lowfire,2),std(lowfire,0,2)/sqrt(6),'-o','LineWidth',2);
    if max(mean(lowfire,2))<max(mean(highfire,2))
        ymax=max(mean(highfire,2))+max(std(highfire,0,2))/sqrt(6)+0.5;
    else
        ymax=max(mean(lowfire,2))+max(std(lowfire,0,2))/sqrt(6)+0.5;
    end

    if min(mean(lowfire,2))<min(mean(highfire,2))
        ymin=min(mean(lowfire,2))-max(std(lowfire,0,2))/sqrt(6)-0.5;
    else
        ymin=min(mean(highfire,2))-max(std(highfire,0,2))/sqrt(6)-0.5;
    end
    title('LowSpeed ');
    xlabel('Speed mm/s');
    ylabel('Fring Rate count/s');
    ylim([ymin,ymax]);
    xlim([5,50]);
    subplot(4,1,3);
    errorbar(highSpeed,mean(highfire,2),std(highfire,0,2)/sqrt(6),'-o','LineWidth',2,'Color','red');
    title(['highSpeed ']);
    xlabel('Speed mm/s');
    ylabel('Fring Rate count/s');
    ylim([ymin,ymax]);
    
    subplot(4,1,4);
    ymax=max(mean(anglefire,2))+max(std(anglefire,0,2))/sqrt(5)+0.5;
    ymin=min(mean(anglefire,2))-max(std(anglefire,0,2))/sqrt(5)-0.5;
    errorbar(uniqAngle,mean(anglefire,2),std(anglefire,0,2)/sqrt(5),'-o','LineWidth',2,'Color','black');
    title(['Direction ']);
    xlabel('Angle ');
    ylabel('Fring Rate count/s');
    ylim([ymin,ymax]);
    xlim([0,180]);
end

