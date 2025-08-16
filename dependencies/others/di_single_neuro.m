%%


clear;clc


sig_oi=[];
sig_di=[];
sig_oiy=[];
sig_dix=[];
%%
nowpath=pwd;

 % workpath='G:\2023_Mn_14325_RightBrain\validData';
workpath='G:\2022_Mn_14365\validData';
cd(workpath);   

weeklist = dir(workpath);
weeklist =weeklist ([weeklist .isdir]);



% num=1;

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
        matData=load(fullfile(sitepath,'results\direction_tun\totaldata.mat'));
        matDatap=load(fullfile(sitepath,'results\direction_tun\Pdata.mat'));
        matDataDi=load(fullfile(sitepath,'results\direction_tun\Oidata.mat'));
        % matData=load(fullfile(sitepath,'results\orientation4_tun\totaldata.mat'));
        % matDatap=load(fullfile(sitepath,'results\orientation4_tun\Pdata.mat'));
        % matDataOi=load(fullfile(sitepath,'results\orientation4_tun\Oidata.mat'));
        %         matData=load(fullfile(sitepath,'results\orientationC_tun\totaldata.mat'));
        % matDatap=load(fullfile(sitepath,'results\orientationC_tun\Pdata.mat'));
        % matDataOi=load(fullfile(sitepath,'results\orientationC_tun\Oidata.mat'));
        matDatap=matDatap.Presults;
        matData=matData.Diresults;
        % matDataOi=matDataOi.Dindex;
        matDataDi=matDataDi.Dindex;
        % load(fullfile(sitepath,'combined/p_data_detail.mat'));
        varNames = fieldnames(matDatap);
        for k = 1:length(varNames)
            varName = varNames{k};
            
            if matDatap.(varName)<0.05
                disp([weekname,',',sitename,',',varName,',','p:',num2str(matDatap.(varName)),'OI:',num2str(matDataOi.(varName))])
                % preferdir(matData.(varName));
                sig_oi=[sig_oi,matDataOi.(varName)];
                sig_oiy=[sig_oiy,matDataDi.(varName)];
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
        matData=load(fullfile(sitepath,'results\direction_tun\totaldata.mat'));
        matDatap=load(fullfile(sitepath,'results\direction_tun\Pdata.mat'));
        matDataDi=load(fullfile(sitepath,'results\direction_tun\Oidata.mat'));
        % matData=load(fullfile(sitepath,'results\orientation_tun\totaldata.mat'));
        % matDatap=load(fullfile(sitepath,'results\orientation_tun\Pdata.mat'));
        matDataOi=load(fullfile(sitepath,'results\orientation4_tun\Oidata.mat'));
        %         matData=load(fullfile(sitepath,'results\orientationC_tun\totaldata.mat'));
        % matDatap=load(fullfile(sitepath,'results\orientationC_tun\Pdata.mat'));
        % matDataOi=load(fullfile(sitepath,'results\orientationC_tun\Oidata.mat'));
        matDatap=matDatap.Presults;
        matData=matData.Diresults;
        matDataOi=matDataOi.Dindex;
        matDataDi=matDataDi.Dindex;
        % load(fullfile(sitepath,'combined/p_data_detail.mat'));
        varNames = fieldnames(matDatap);
        for k = 1:length(varNames)
            varName = varNames{k};
            
            if matDatap.(varName)<0.05
                disp([weekname,',',sitename,',',varName,',','p:',num2str(matDatap.(varName)),'DI:',num2str(matDataDi.(varName))])
                % preferdir(matData.(varName));
                sig_di=[sig_di,matDataDi.(varName)];
                sig_dix=[sig_dix,matDataOi.(varName)];
            end

        end
    end
end
% 
% dtotalp2s=totalp2s;
% dtotals2p=totals2p;
cd(nowpath)

%%
figure;
set(gcf,'unit','normalized','position',[0,0,0.2,0.4]);
hold on;
    
    % 绘制 sig_di 的数据点，标记为 'o'
    plot(sig_dix, sig_di, 'ro', 'MarkerSize', 8, 'DisplayName', 'DI Points');
    
    % 绘制 sig_oi 的数据点，标记为 'x'
    plot(sig_oi,sig_oiy, 'bx', 'MarkerSize', 8, 'DisplayName', 'OI Points');
    
    % 添加 45 度虚线
    plot([0 1], [0 1], 'k--', 'LineWidth', 1);
    
    % 设置坐标轴范围为 0 到 1
    xlim([0 1]);
    ylim([0 1]);
    yticks([0,0.5,1])
    % 添加标签和标题
    xlabel('OI (Orientation Index)');
    ylabel('DI (Direction Index)');
    title('DI vs OI Scatter Plot');
    set(gca,'FontSize',15)
set(gca,'LineWidth',2)
set(gca,'FontWeight','bold')
    % 添加图例
    legend('show');
%%
% uniqSpeed=[5,10,20,25,40,50,80,100,160,200,240,320];
Fs=1000;
uniqAngle=[0,45,90,135,180];% 0,45,90,135,180 per 24times

res=zeros(5,12);
for i=1:numel(uniqAngle)
    for j=1:numel(uniqSpeed)
        rc=find(StimSpeed10d==uniqSpeed(j)&StimAngle10d==uniqAngle(i));
        res(i,j)=numel(rc);
    end
end
res

%%
clear
clc
close all

matData=load("p_data_cut_inc.mat");

uniqSpeed=[5,10,20,25,40,50,80,100,160,200,240,320];
Fs=1000;
uniqAngle=[0,45,90,135,180];% 0,45,90,135,180 per 24times

lowspeed=uniqSpeed(1:6);
highspeed=uniqSpeed(7:12);

varName='SPKC09b';
totalspkc = matData.(varName);
totalangle=matData.(['StimAngle',varName(end-2:end)]);
totalspeed=matData.(['StimSpeed',varName(end-2:end)]);
totallfp=matData.(['LFP',varName(end-2:end-1)]);
lfporder=matData.(['Order',varName(end-2:end-1)]);
[totallfp,totalspkc] = match_trials(totallfp, totalspkc, lfporder, matData.(['Order',varName(end-2:end)]));
totalfr=cell(5,12);
totalsd=cell(5,12);
totalact=cell(12,5);
for i=1:5
    for j=1:12
        rc=find(totalangle==uniqAngle(i)&totalspeed==uniqSpeed(j));
        refr=[];
        act=[];
        if ~isempty(rc)
            for k=1:length(rc)
                act=[act;totalspkc(rc(k)).times];
                spkcount_after=binspikes(totalspkc(rc(k)),Fs,[1,2]);
                fr_after=sum(spkcount_after)/1;
                % spkcount_before=binspikes(totalspkc(rc(k)),Fs,[0.5,0.75]);
                % fr_before=sum(spkcount_before)/0.25;
                %refr=[refr,(fr_after-fr_before)/(fr_before+eps)];
                %refr=[refr,(fr_after)/(fr_before+eps)];
                refr=[refr,fr_after];
            end
            totalact(j,i)={act};
            totalfr(i,j)={refr};
        else
            totalfr(i,j)={0};
            totalact(j,i)={[]};
        end
    end
end

y_base = 1.2;  % 起始的 y 值
y_spacing = 1.2;  % 每个 trial 之间的间隔
    figure;
    set(gcf,'unit','normalized','position',[0,0.2,0.18,0.2]);
hold on;
for j=1:5
    spike_times=[];

    for i=1:6
        spike_times=totalact{i,j};
        spike_times=spike_times(spike_times>=1 &spike_times<=2);
        y_offset = y_base + (i-1) * y_spacing;
    
    for h = 1:length(spike_times)
        % 使用 plot 函数画竖条
        plot([spike_times(h)+(j-1)*1.2, spike_times(h)+(j-1)*1.2], [y_offset, y_offset + 1], 'k', 'LineWidth', 1.5,'Color',[1 0.8 0.25]);
    end

    end

    %hold off
end
% hold off
xlim([1 8]);
% xticks([1,2,2.5,3.5,4,5,5.5,6.5,7,8])
ylim([1 8.5])

y_base = 1.2;  % 起始的 y 值
y_spacing = 1.2;  % 每个 trial 之间的间隔
%     figure;
% hold on;
for j=1:5
    spike_times=[];

    for i=7:12
        spike_times=totalact{i,j};
        spike_times=spike_times(spike_times>=1 &spike_times<=2);
        y_offset = y_base + (i-1) * y_spacing;
    
    for h = 1:length(spike_times)
        % 使用 plot 函数画竖条
        plot([spike_times(h)+(j-1)*1.2, spike_times(h)+(j-1)*1.2], [y_offset, y_offset + 1], 'k', 'LineWidth', 1.5,'Color',[0 0.5 0.8]);
    end

    end
xlim([1 7]);
xticks([1,2,2.2,3.2,3.4,4.4,4.6,5.6,5.8,7])
xticks([])
yticks([])
box off
 ylim([1 16])

    %hold off
end

xx=[2.1,3.3,4.5,5.7]
    for h = 1:length(xx)
        % 使用 plot 函数画竖条
        plot([xx(h), xx(h)], [1,16], 'k', 'LineWidth', 1.2,'Color',[0.5 0.5 0.5]);
    end
    axis off


fr=[];
sdd=[];
for j=1:5
    spike_times=[];
        rc=find(totalangle==uniqAngle(j)&totalspeed<uniqSpeed(7));
                refr=[];
                if ~isempty(rc)
                    for k=1:length(rc)
                        spkcount_after=binspikes(totalspkc(rc(k)),Fs,[1,2]);
                        fr_after=sum(spkcount_after)/1;
                        refr=[refr,fr_after];
                    end
                    fr=[fr,mean(refr)];
                    sdd=[sdd,std(refr) / sqrt(size(refr,2))];
                else
				%实际上已经不可能是0，我已经剔除掉了在一个方向上没有任何trial的神经元
                    fr=[fr,0];
                    sdd=[sdd,0];
                end
end
figure;
set(gcf,'unit','normalized','position',[0,0.2,0.15,0.25]);
hold on
errorbar(uniqAngle,fr,sdd, 'LineWidth', 1.5,'Color',[1 0.8 0.25]);

 sin_component = sum(fr .* sin(deg2rad(uniqAngle)));
        cos_component = sum(fr .* cos(deg2rad(uniqAngle)));
        numerator = sqrt(sin_component^2 + cos_component^2);
        
        % 计算分母部分
        denominator = sum(fr);
        
        % 计算方向指数 (DI)
        DI = numerator / denominator
% xlim([1 8]);
% xticks([1,2,2.5,3.5,4,5,5.5,6.5,7,8])
% ylim([1 8.5])

fr=[];
sdd=[];
for j=1:5
    spike_times=[];
        rc=find(totalangle==uniqAngle(j)&totalspeed>uniqSpeed(6));
                refr=[];
                if ~isempty(rc)
                    for k=1:length(rc)
                        spkcount_after=binspikes(totalspkc(rc(k)),Fs,[1,2]);
                        fr_after=sum(spkcount_after)/1;
                        refr=[refr,fr_after];
                    end
                    fr=[fr,mean(refr)];
                    sdd=[sdd,std(refr) / sqrt(size(refr,2))];
                else
				%实际上已经不可能是0，我已经剔除掉了在一个方向上没有任何trial的神经元
                    fr=[fr,0];
                    sdd=[sdd,0];
                end
end
errorbar(uniqAngle,fr,sdd, 'LineWidth', 1.5,'Color',[0 0.5 0.8]);

 sin_component = sum(fr .* sin(deg2rad(uniqAngle)));
        cos_component = sum(fr .* cos(deg2rad(uniqAngle)));
        numerator = sqrt(sin_component^2 + cos_component^2);
        
        % 计算分母部分
        denominator = sum(fr);
        
        % 计算方向指数 (DI)
        DI = numerator / denominator
ylim([1 10]);
% xticks([1,2,2.5,3.5,4,5,5.5,6.5,7,8])
% ylim([1 8.5])

fr=[];
sdd=[];
for j=1:5
    spike_times=[];
        rc=find(totalangle==uniqAngle(j));
                refr=[];
                if ~isempty(rc)
                    for k=1:length(rc)
                        spkcount_after=binspikes(totalspkc(rc(k)),Fs,[1,2]);
                        fr_after=sum(spkcount_after)/1;
                        refr=[refr,fr_after];
                    end
                    fr=[fr,mean(refr)];
                    sdd=[sdd,std(refr) / sqrt(size(refr,2))];
                else
				%实际上已经不可能是0，我已经剔除掉了在一个方向上没有任何trial的神经元
                    fr=[fr,0];
                    sdd=[sdd,0];
                end
end
errorbar(uniqAngle,fr,sdd, 'LineWidth', 1.5,'Color',[0.5 0.5 0.5]);

 sin_component = sum(fr .* sin(deg2rad(uniqAngle)));
        cos_component = sum(fr .* cos(deg2rad(uniqAngle)));
        numerator = sqrt(sin_component^2 + cos_component^2);
        
        % 计算分母部分
        denominator = sum(fr);
        
        % 计算方向指数 (DI)
        DI = numerator / denominator
ylim([0 4]);
xticks(uniqAngle)
% yticks([0,5,10])
yticks([0 2 4])
set(gca,'FontSize',15)
set(gca,'FontWeight','bold')
set(gca,'LineWidth',2)
xlim([0 180])
%%
for j=1:5
    spike_times=[];
        rc=find(totalangle==uniqAngle(j));
                refr=[];
                if ~isempty(rc)
                    anglelfp=totallfp(:,rc);
                    lfp=lfp_norm(anglelfp); %norm lfp
    
                    spectro=[];
                    % multi taper method:
                    [spectro,t,frequency]=mtspecgramc(lfp,movingwin,params);
                    spectro=spectro';
                    meanidx=t>=0.5 &t<=0.75;
                    meanspec=spectro(:,meanidx);
                    meanspec=mean(meanspec,2);
                    spectro=spectro./meanspec;% normalize
                    specData= 10*log10(sepctro);

                    delta(num,:)=mean(specData(frequency>=1&frequency<=3,:),1);
                    theta(num,:)=mean(specData(frequency>=4&frequency<=8,:),1);
                    alpha(num,:)=mean(specData(frequency>=9&frequency<=12,:),1);
                    beta(num,:)=mean(specData(frequency>=12&frequency<=24,:),1);
                    gamma(num,:)=mean(specData(frequency>=30&frequency<=50,:),1);
                    

                    figure;
                    set(gcf,'unit','normalized','position',[0,0,0.20,0.25]);
                    hold on
                    % 设置长方形区域的坐标和大小
                    x = 0.2; % 左下角x坐标
                    y = 0.2; % 左下角y坐标
                    width = 0.9;  % 长方形的宽度
                    height = 0.6; % 长方形的高度
                    % 添加灰色长方形区域
                    rectangle('Position', [x, y, width, height], 'FaceColor', [0.95, 0.95, 0.95], 'EdgeColor', 'none');
                    
                    yyaxis left
                    meanfr=mean(totalfr,1);
                    sem=std(totalfr,1)/sqrt(size(totalfr,1));
                    s3=shadedErrorBar(t-1,meanfr,sem,'lineProps',{'Color','black','LineWidth',2});% 绘制发放率变化曲线
                    % ylabel('Firing Rate (Hz)');
                    xlim([-0.2 1.2]);
                     % ylim([0, 0.6]); 
                     ylim([0, 1]); 
                     % yticks([0:0.2:0.6])
                    
                    
                    
                    
                    
                    yyaxis right
                    meandelta=mean(delta,1);
                    sem=std(delta,1)/sqrt(size(delta,1));
                    shadedErrorBar(t-1,meandelta,sem,'lineProps',{'Color',[0.5 0.8 0.82],'LineWidth',2,'LineStyle','-','Marker', 'none'});
                    meantheta=mean(theta,1);
                    sem=std(theta,1)/sqrt(size(theta,1));
                    shadedErrorBar(t-1,meantheta,sem,'lineProps',{'Color',[1 0.3 0.5],'LineWidth',2,'LineStyle','-','Marker', 'none'});
                    meanalpha=mean(alpha,1);
                    sem=std(alpha,1)/sqrt(size(alpha,1));
                    shadedErrorBar(t-1,meanalpha,sem,'lineProps',{'Color',[0.8 0.8 0.2],'LineWidth',2,'LineStyle','-','Marker', 'none'});
                    meanbeta=mean(beta,1);
                    sem=std(beta,1)/sqrt(size(beta,1));
                    shadedErrorBar(t-1,meanbeta,sem,'lineProps',{'Color',[0.5 0.3 0.2],'LineWidth',2,'LineStyle','-','Marker', 'none'});
                    meangamma=mean(gamma,1);
                    sem=std(gamma,1)/sqrt(size(gamma,1));
                    shadedErrorBar(t-1,meangamma,sem,'lineProps',{'Color',[0.1 0.8 0.2],'LineWidth',2,'LineStyle','-','Marker', 'none'});
                    
                    xlim([-0.5 1.2]);
                    xticks([-0.2:0.2:1.2]);
                    yticks([-1:1:3])
                    ylim([-1 2])
                    
                    % legend('δ','θ','α','β','γ');
                    yyaxis right; % 使用右侧的 y 轴
                    ax = gca;        % 获取当前轴
                    ax.YColor = 'k'; 
                    ax.LineWidth=2;
                    ax.FontWeight="bold"
                    ax.FontSize=15;
                    ax.Box='on'
                else
				%实际上已经不可能是0，我已经剔除掉了在一个方向上没有任何trial的神经元
                end
end