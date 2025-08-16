function combine_speed_angle(sitepath)
close all;
combinepath=fullfile(sitepath,'combined');
cd(combinepath);
load Combined.mat

savePath = fullfile(sitepath,'results\speed_angle_con');
mkdir(savePath);
%%
Fs = 1000;                  %Sampling frequency
TW = 3;						%Choose time-bandwidth product of 3.
ntapers = 2*TW-1;			%Choose the # of tapers.
fpass = [0 30];
pad = 2;

params.Fs = Fs;				%... sampling frequency
params.tapers=[TW,ntapers];	%... time-bandwidth product & tapers.
params.fpass=fpass;       %... fpass
params.pad = pad;             %... padding.
params.trialave = 1;		%... trial average.
%%
uniqSpeed=[5,10,20,25,40,50,80,100,160,200,240,320];
lowSpeed=uniqSpeed(1:6);% speed 5,10,20,25,40,50 per 10 times
highSpeed=uniqSpeed(7:12);% speed 80,100,160,200,240,320 per 10 times

uniqAngle=[0,45,90,135,180];% 0,45,90,135,180 per 24times

%%
vars = who;
numbers_channel = {};
channels_name={};
for i = 1:length(vars)
    matches = regexp(vars{i}, '\d+', 'match');
    if ~isempty(matches)
        numbers_channel = [numbers_channel, matches];
        channels_name=[channels_name,['channel',matches{1}]];
    end
end
numbers_channel = unique(numbers_channel);
channels_name=unique(channels_name);
%%
% aquire 0.2s-1s 实际是1.2s-2s 
recordtime=-1:1/Fs:1.5-1/Fs; % 实际变为0-2.5s
timeindx=recordtime>=0.2 & recordtime <=1 | abs(recordtime-0.2)<eps | abs(recordtime-1)<eps;
needtime=recordtime(timeindx);
vars_lfp = who('-regexp', '^LFP');
vars_spk = who('-regexp', '^SPKC');
for k=1:numel(vars_spk)
    eval(['lfp=',vars_lfp{k},';']);
    eval(['spk=',vars_spk{k},';']);
    for kk=1:numel(spk)
        tmpts=spk(kk).times;
        tmpindx=tmpts>=1.2 & tmpts<=2 | abs(tmpts-1.2)<eps | abs(tmpts-2)<eps;
        tmpts=tmpts(tmpindx);
        % extratime=1.2;%onset变为0
        % spk(kk).times=tmpts-extratime;   
        spk(kk).times=tmpts;   
    end
    lfp=lfp(timeindx,:);
    eval([vars_lfp{k},'=lfp;']);
    eval([vars_spk{k},'=spk;']);
end
disp('now ts is 0.2-1s!(actually 1.2-2s) ')
%disp('now ts is 0.2-1s! and onset is 0 point!!')
%%
indxL=startsWith(vars, 'LFP');
indxS=startsWith(vars, 'SPKC');
indxAngle=startsWith(vars, 'StimAngle');
indxSpeed=startsWith(vars, 'StimSpeed');
nameL=vars(find(indxL==1));
nameS=vars(find(indxS==1));
nameAngle=vars(find(indxAngle==1));
nameSpeed=vars(find(indxSpeed==1));

time=[1.2 2];
%%
totalresults=struct();
for n=1:length(nameS)
    eval(['spk=',nameS{n},';']);
    spkname=nameS{n};
    spkname=spkname(5:end);
    eval(['lfp=',nameL{n},';']);
    eval(['TotalAngle=',nameAngle{n},';']);
    eval(['TotalSpeed=',nameSpeed{n},';']);

    rename=['result',spkname];
    eval([rename,'(5,12).fr=[];']);
    %
    for i=1:12
        for j=1:5
            rc=find(TotalAngle==uniqAngle(j)&TotalSpeed==uniqSpeed(i));
            refr=[];
            relfp=[];
            respk=[];
            respks=[];
            for k=1:length(rc)
                tmpspk=spk(rc(k)).times;
                tmpspks.times=tmpspk;
                spkcount=binspikes(spk(rc(k)),Fs,time);
                tmplfp=lfp(:,rc(k));
                if isempty(tmpspk)
                    fr=0;
                else
                    count=sum(spkcount(:));
                    fr=count/(time(end)-time(1));
                end
                refr=[refr,fr];
                %respk=[respk;spkcount];
                relfp=[relfp,tmplfp];
                respks=[respks,tmpspks];
            end
            %[C,phi,~,Sf,Ss,f]=coherencycpb(relfp,respk,params);
            eval([rename,'(j,i).fr=refr;'])
            %eval([rename,'(j,i).spkcount=respk;'])
            eval([rename,'(j,i).spkts=respks;'])
            eval([rename,'(j,i).lfp=relfp;'])
            %eval([rename,'(j,i).coherence=C;'])
            %eval([rename,'(j,i).fcohe=f;'])
            eval(['totalresults.',rename,'=',rename,';'])

        end
    end

end
save(fullfile(savePath,'results.mat'),'-struct',"totalresults");
% %%
% close all;
% %画sin函数拟合
% for n=1:length(nameS)
%     eval(['spk=',nameS{n},';']);
%     spkname=nameS{n};
%     spkname=spkname(5:end);
%     eval(['lfp=',nameL{n},';']);
%     eval(['TotalAngle=',nameAngle{n},';']);
%     eval(['TotalSpeed=',nameSpeed{n},';']);
% 
%     rename=['result',spkname];
%     eval([rename,'(5,12).fr=[];']);
%     %
%     figure;
%     sgtitle(nameS{n});
%     for i=1:3
%             rc=find(TotalSpeed==uniqSpeed(i));
%             respk=[];
%             for k=1:length(rc)
%                 tmpspk=spk(rc(k)).times;
%                 tmpspks.times=tmpspk;
%                 spkcount=binspikes(spk(rc(k)),Fs,time);
%                 respk=[respk,spkcount];
%             end
%             respk=sum(respk,2)';
%             sinspk=sin(2*pi*needtime*uniqSpeed(i))+1;
%             subplot(6,1,i);
%             plot(needtime,respk);
%             hold on
%             plot(needtime,sinspk);
%             hold off
%             title(['speed',num2str(uniqSpeed(i))]);
%     end
% 
% end
% %%
% close all;
% 
% vars=who;
% indxre=startsWith(vars, 'result');
% namesre=vars(find(indxre==1));
% for i=1:numel(namesre)
%     tmpname=namesre{i};
%     eval(['tmpre=',tmpname,';']);
%     for j=1:12
%         %speed为变量研究angle
%         for k=1:5
%             eval(['plotfire(k,j)=mean(',tmpname,'(k,j).fr);']);
%         end
%     end
% 
%     plotfire=log(plotfire+1);
%     %subplot(1,2,1);
%     z=plotfire;
%     x=uniqSpeed;
%     y=uniqAngle;
%     [Xq, Yq] = meshgrid(min(x):5:max(x), min(y):5:max(y));
%     % 二维插值
%     Zq = interp2(x, y, z, Xq, Yq, 'spline'); % 使用样条插值
%     % subplot(1,2,1);
%     [xqs,yqs]=meshgrid(x,y);
%     plotfire=[];
% 
%     %寻找peak
%     spot=Zq;
%     spot(spot<0.7*max(spot(:)))=0;
%     PeaksMap = imregionalmax(spot);
%     [maxRow,maxCol]=find(PeaksMap==1);
%     peaks=zeros(1,length(maxCol));
% 
%     figure;
%     set(gcf,'unit','normalized','position',[0,0,1,0.62]);
%     subplot(1,2,1);
%     surf(Xq, Yq, Zq);
%     xticks(uniqSpeed);
%     yticks([0,45,90,135,180]);
%     title('Firing Rate');
%     xlabel('Speed(mm/s)')
%     title(['Neuron ',extractAfter(tmpname,'result')]);
%     ylabel('Angle(°)');
%     %shading interp;
% 
%     spksCohe=[];
%     lfpsCohe=[];
%     hold on
%     for ii=1:length(maxRow)
%         peaks(ii)=Zq(maxRow(ii),maxCol(ii));
%         peakx(ii)=Xq(1,maxCol(ii));
%         peaky(ii)=Yq(maxRow(ii),1);
%         distances = abs(peakx(ii) - uniqSpeed);
%         sindex=find(distances==min(distances));
%         distances = abs(peaky(ii) - uniqAngle);%找最接近peak的angle以及speed
%         aindex=find(distances==min(distances));
% 
%         for jj=1:numel(sindex)
%             for kk=1:numel(aindex)
%                 spkts=tmpre(aindex(kk),sindex(jj)).spkts;
%                 dN=binspikes(spkts,Fs,[1.2 2]);
%                 spksCohe=[spksCohe,dN];
%                 lfpsCohe=[lfpsCohe,tmpre(aindex(kk),sindex(jj)).lfp];
%             end
%         end
%         plot3(Xq(1,maxCol(ii)),Yq(maxRow(ii),1),peaks(ii),'k.','MarkerSize',20);
%         text(Xq(1,maxCol(ii)),Yq(maxRow(ii),1),peaks(ii),['Peak ',num2str(ii),' (',num2str(Xq(1,maxCol(ii))),',',num2str(Yq(maxRow(ii),1)),',',num2str(peaks(ii)),')']);
%     end
%     hold off
% 
%     lfpsCohe=lfp_norm(lfpsCohe);
%     [Cohe,phi,~,Sf,Ss,f]=coherencycpb(lfpsCohe,spksCohe,params);
% 
%     subplot(1,2,2);
%     plot(f, Cohe, 'k', 'LineWidth', 2);xlim([0 30]);xlabel('Frequency/Hz');ylabel('Coherence');title('Coherence');set(gca, 'FontSize', 14);
%     xticks([0,4,8,12,16,21,30]);
% 
%     fileName = ['neuron',extractAfter(namesre{i},'result'),'.png'];
%     % 使用saveas保存图形到指定路径
%     saveas(gcf, fullfile(savePath, fileName));
% end

    %%
% %%
% for i=1:numel(namesre)
%     column_data = hist{:, i+1};
%     rdata=column_data(rindx);
%     tmpname=namesre{i};
%     psthname=namesre{i};
%     psthname=psthname(3:end);
%     lfpname=['LFP',psthname(1:2)];
% 
%     % figure;
%     % %set(gcf,'unit','normalized','position',[0,0,0.17,0.62]);
%     % 
%     % wave = readtable([trialname,'_SPKC',tmpname(end-2:end),'.csv']);
%     % %subplot(3,1,1);
%     % shadedErrorBar(wave{:,1},wave{:,2},wave{:,3},'lineprops','-r');
%     % title(['Waveform',psthname]);
%     % xlabel('time/s');
%     % ylabel('v');
% 
%     figure;
%     set(gcf,'unit','normalized','position',[0,0,1,0.5]);
%     for j=1:5
%         %speed为变量研究angle
%         for k=1:12
%             eval(['plotfire(k)=',tmpname,'(j,k).fr;']);
%             eval(['plotlfp(:,k)=',tmpname,'(j,k).lfp;']);
%             eval(['plotspk(k,:)=',tmpname,'(j,k).spkcount;']);
%         end
% 
%     ax=subplot(2,5,j);
%     meanlfp=mean(plotlfp,2);
% 
%     plot(recordtime,meanlfp-mean(meanlfp));
%     hold on;
%     plot(recordtime,mean(plotspk,1));
% 
%     hold off;
%     ylim([min(meanlfp-mean(meanlfp))-0.01,max(meanlfp-mean(meanlfp))+0.01]);
%     line(ax,[0.2,0.2],ylim,'Color','red','LineWidth',1.5);
%     line(ax,[1,1],ylim,'Color','blue','LineWidth',1.5);
%     title(['Angle ',num2str(uniqAngle(j))]);
%     xlim([-0.5,1.5]);
%     xlabel('Time/s');
%     %xticks([0.2,0.4,0.6,0.8,1.0]);
%     ylabel('Volt / Count');
%     % ax=subplot(4,1,3);
%     % bar(ax,bin_centers, rdata, 'hist');
%     % line(ax,[0.2,0.2],ylim,'Color','red','LineWidth',2);
%     % line(ax,[1,1],ylim,'Color','Blue','LineWidth',2);
%     % xlim([0 1.5]);
%     % 
%     % title(['PSTH ',psthname]);
%     % xlabel('time/s');
%     % ylabel('count');
% 
%     subplot(2,5,j+5);
%     plot(uniqSpeed,plotfire,'-o');
%     xlabel('Speed mm/s');
%     ylabel('Firng Rate');
%     xticks([5,10,20,25,40,50,80,100,160,200,240,320]);
%     xlim([5,320]);
%     ylim([0,max(plotfire)+0.5]);
%     %xticklabels({'0','45','90','135','180'});
%     end
%     sgtitle(['Neuron ',tmpname(end-2:end)]);
%     % windowsize=3;
%     % smoothfr=smoothdata(fire,'gaussian',windowsize);
%     % %imagesc(smoothfr);
%     % imagesc(fire);
%     % xticks(1:12);
%     % xticklabels({'5','10','20','25','40','50','80','100','160','200','240','320'});
%     % yticks(1:5);
%     % yticklabels({'0','45','90','135','180'});
%     % %colormap("jet");
%     % colorbar;
%     % title('Speed-Angle FR');
%     % xlabel('Speed mm/s');
%     % ylabel('Angle');
% 
% end
%     cd(currentDir);
% 
