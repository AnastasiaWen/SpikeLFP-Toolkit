%clear;
clc;
trialpath=pwd;
% load(fullfile(trialpath,'data_car.mat'));
savepath=fullfile(trialpath,'ica_re');
mkdir(savepath);

% 使用ICA分解信号
close all
clc
nComponents=15;
Fs=1000;
vars=who;
vars_LFP=who('-regexp','^LFP');
vars_before=who('-regexp','^before');
vars_spk=who('-regexp','^SPKC');

for trr=1:size(LFP05,2)
    
    lim=0.1;
    orig_sig=[];
    re_ica=1;
    lfp_data=[];
    close all;
    clearvars hw sw h1 h2 c1 c2 af1 af2
    hw=figure;
    set(gcf,'unit','normalized','position',[0,0,0.4,1]);
    ha=tight_subplot(8,2,0.04);
    % sgtitle('waveform');
    sw=figure;
    set(gcf,'unit','normalized','position',[0.4,0,0.4,1]);
    fa=tight_subplot(8,2,0.04);
    % sgtitle('spectrum');
    for ii = 1:numel(vars_LFP)
        cname=vars_LFP{ii};
        cname=['channel',cname(end-1:end)];
        tmp_LFP=eval(vars_LFP{ii});
        tmp_LFP=tmp_LFP(:,trr);

                %% 带通滤波 0.5Hz-300Hz
d = designfilt('bandpassfir', 'FilterOrder', 300, ...
       'CutoffFrequency1', 1, 'CutoffFrequency2', 300,...
       'SampleRate', Fs,'Window','hamming');
% d=designfilt('bandpassiir', ...
%     'FilterOrder',4,'HalfPowerFrequency1',0.5, ...
%     'HalfPowerFrequency2',300,'SampleRate',Fs);
tmp_LFP=filtfilt(d,tmp_LFP);

        axes(ha(ii));
        plot(-1:1/1000:1.4999,tmp_LFP);
        title(cname);
        ylim([-lim lim]);

        [faxis,Sxx]=fftspecc(tmp_LFP);
        axes(fa(ii));
        plot(faxis,10*log10(Sxx));
        xlim([0 100]);
        title(cname);



        lfp_data=[lfp_data,tmp_LFP];
        ylim([-80 0]);
    end
    
    prompt="回车跳过";
    jump= input(prompt,'s');

    if isempty(jump)
     continue;
    end

    f_order=[];
    prompt="额外的filter:";
    f_or = input(prompt,'s');

    f_order=[50 150];
    while ~isempty(f_or)
        f_order = str2num(f_or);
        %lfp_data=[];
        close all
        hw=figure;
        set(gcf,'unit','normalized','position',[0,0,0.4,1]);
        ha=tight_subplot(8,2,0.04);
        % sgtitle('waveform');
        sw=figure;
        set(gcf,'unit','normalized','position',[0.4,0,0.4,1]);
        fa=tight_subplot(8,2,0.04);
        % sgtitle('spectrum');


        for ii = 1:numel(vars_LFP)
        cname=vars_LFP{ii};
        cname=['channel',cname(end-1:end)];
        tmp_LFP=lfp_data(:,ii);

        %optional filt
        for ff=1:numel(f_order)
        ffilt = f_order(ff); % 欲消除的频率为50Hz
        Q = 100; % 带宽参数
        w = ffilt/(Fs/2);
        bw = w/Q;
        [b, a] = iirnotch(w, bw);
        tmp_LFP=filtfilt(b,a,tmp_LFP);
        end

        axes(ha(ii));
        plot(-1:1/1000:1.4999,tmp_LFP);
        title(cname);
        ylim([-0.08 0.08]);

        [faxis,Sxx]=fftspecc(tmp_LFP);
        axes(fa(ii));
        plot(faxis,10*log10(Sxx));
        ylim([-80 0])
        xlim([0 100]);
        title(cname);
        lfp_data(:,ii)=tmp_LFP;
        end
        prompt="额外的filter:";
    f_or = input(prompt,'s');
    end
prompt="调整ylim:";
li = input(prompt,'s');

while ~isempty(li)
lim = str2num(li);
axesHandles = findall(hw, 'Type', 'axes');  % 获取当前 figure 中所有的轴
set(axesHandles, 'YLim', [-lim lim]);
prompt="调整ylim:";
li = input(prompt,'s');
end
orig_sig=lfp_data;
saveas(hw,fullfile(savepath,['tr',num2str(trr),'orig.jpg']));
saveas(sw,fullfile(savepath,['tr',num2str(trr),'orig_s.jpg']));

cmrre=lfp_data;
meancmr=mean(cmrre,2);
c1=figure;
set(gcf,'unit','normalized','position',[0,0,0.4,1]);
ca1=tight_subplot(8,2,0.04);
% sgtitle('waveform');
c2=figure;
set(gcf,'unit','normalized','position',[0.4,0,0.4,1]);
ca2=tight_subplot(8,2,0.04);
for ii = 1:size(cmrre,2)
    cmrre(:,ii)=cmrre(:,ii)-meancmr;
    cname=vars_LFP{ii};
    cname=['CMR',cname(end-1:end)];

     axes(ca1(ii));
        plot(-1:1/1000:1.4999,cmrre(:,ii));
        title(cname);
        ylim([-lim  lim]);

        [faxis,Sxx]=fftspecc(cmrre(:,ii));
        axes(ca2(ii));
        plot(faxis,10*log10(Sxx));
        xlim([0 100]);
        title(cname);
        ylim([-80 0]);
end


while ~isempty(re_ica)
    clean_lfp=[];
if re_ica==0
    lfp_data=orig_sig;
    close(hw);close(sw);close(af1);close(af2);
    hw=figure;
    set(gcf,'unit','normalized','position',[0,0,0.4,1]);
    ha=tight_subplot(8,2,0.04);
    sw=figure;
    set(gcf,'unit','normalized','position',[0.4,0,0.4,1]);
    fa=tight_subplot(8,2,0.04);
    for ii = 1:size(lfp_data,2)
        cname=vars_LFP{ii};
        cname=['channel',cname(end-1:end)];
        tmp_LFP=lfp_data(:,ii);
        axes(ha(ii));
        plot(-1:1/1000:1.4999,tmp_LFP);
        title(cname);
        ylim([-0.08 0.08]);

        [faxis,Sxx]=fftspecc(tmp_LFP);
        axes(fa(ii));
        plot(faxis,10*log10(Sxx));
        xlim([0 100]);
        ylim([-80 0])
        title(cname);
    end
    prompt="调整ylim:";
    li = input(prompt,'s');

    while ~isempty(li)
    lim = str2num(li);
    axesHandles = findall(hw, 'Type', 'axes');  % 获取当前 figure 中所有的轴
    set(axesHandles, 'YLim', [-lim lim]);
    prompt="调整ylim:";
    li = input(prompt,'s');
    end
end
if exist('h1')
    close(h1);
close(h2);
end
% rica
% close all
[coeff, score, latent] = pca(zscore(lfp_data));
nComponents=size(lfp_data,2);
icasig=[];
%while size(icasig,2)<nComponents
% ICA
%rica_models= rica(score, nComponents); % MATLAB's reconstructionICA
rica_models= rica(lfp_data, nComponents);
% Get independent components
icasig = transform(rica_models, lfp_data);
%end



% %% FastICA
% % Goal of ICA is to minimise Gaussianity of the sources
% % Note - FastICA takes matrix dimensions (signals x samples), opposite to
% %   the general dimensions used for MATLAB functions (samples x signals)
% close all
% data =lfp_data';
% [coeff, score, latent] = pca(zscore(lfp_data));
% 
% [icasig, A, W]= fastica(data);
% 
% % icasig=[];
% % %每行是观测值
% % % while size(icasig,1)<size(data,1)
% % while size(icasig,1)<nComponents
% % [icasig, A, W]= fastica(score','maxNumIterations',1000, 'numOfIC', nComponents);
% % end
% 
% icasig=icasig';
%% % 3. 可视化独立成分的频谱
h1=figure;
set(gcf,'unit','normalized','position',[0,0,0.4,1]);
cw=tight_subplot(8,2,0.04);
% sgtitle('ICA wave');
h2=figure;
set(gcf,'unit','normalized','position',[0.4,0,0.4,1]);
cs=tight_subplot(8,2,0.04);
low_power=[];
% sgtitle('ICA spec');
for i = 1:size(icasig,2)
    axes(cw(i));
    % plot(icasig(:,i));
   plot(-1:1/1000:1.4999,icasig(:,i))
   ylim([-lim lim]);
    title(['Independent Component ' num2str(i)]);
   % figure;
    [faxis,Sxx]=fftspecc(icasig(:,i));
    axes(cs(i));
    plot(faxis,10*log10(Sxx));
    low_freq_idx = faxis <= 30;
    low_power = [low_power,sum(Sxx(low_freq_idx))];
   
    xlim([0 100]);
    ylim([-80 0]);
    title(['spectrum ' num2str(i)]);
end
indx=find(low_power>=(mean(low_power)+std(low_power)));
[~,sind]=sort(low_power);
sind=sind';
for hh=1:numel(indx)
    axes(cw(indx(hh)));
    title(['**Independent Component ' num2str(indx(hh))],'Color','r');
    axes(cs(indx(hh)));
    title(['**spectrum ' num2str(indx(hh))],'Color','r');
end
% 
% prompt="调整ylim:";
% lim = input(prompt,'s');
% lim = str2num(lim);
% while ~isempty(lim)
% axesHandles = findall(h1, 'Type', 'axes');  % 获取当前 figure 中所有的轴
% set(axesHandles, 'YLim', [-lim lim]);
% prompt="调整ylim:";
% lim = input(prompt,'s');
% lim = str2num(lim);
% end
%%
%close all
disp(strcat("10Hz低频power参考值排名：",num2str(sind)));
prompt=strcat("请输入去除的成分序号:(推荐值",num2str(indx),")   ");
% % 4. 识别噪声成分，假设第 5 个成分是噪声
noise_components = input(prompt,'s');
if noise_components=='y'
    noise_components=indx;
else
noise_components = str2num(noise_components);
end
if ~isempty(noise_components)
    % % 5. 将噪声成分置零
    % icasig(noise_components, :) = 0;
    % 
    % % 6. 重构去噪后的 LFP 信号
    % clean_data = A * icasig;  % [16 × 时间点]
    % clean_lfp = clean_data';  % [时间点 × 16]
    
    % 使用 rica 模型中的系数重组信号
    icasig(:, noise_components) = 0;  % 将第 k 个成分设为零
    clean_lfp = icasig * rica_models.TransformWeights';  % 重组信号
    
    %%
    %close(hw);close(sw);
    af1=figure;
    set(gcf,'unit','normalized','position',[0,0,0.4,1]);
    ha=tight_subplot(8,2,0.04);
    af2=figure;
    set(gcf,'unit','normalized','position',[0.4,0,0.4,1]);
    fa=tight_subplot(8,2,0.04);
    for ii = 1:size(clean_lfp,2)
        cname=vars_LFP{ii};
        cname=['channel',cname(end-1:end)];
        tmp_LFP=clean_lfp(:,ii);
        axes(ha(ii));
        plot(-1:1/1000:1.4999,tmp_LFP);
        title(cname);
        ylim([-lim lim]);

        [faxis,Sxx]=fftspecc(tmp_LFP);
        axes(fa(ii));
        plot(faxis,10*log10(Sxx));
        xlim([0 100]);
        ylim([-80 0])
        title(cname);
    end
    lfp_data=clean_lfp;
    prompt="调整ylim:";
    li = input(prompt,'s');
    lim=0.08;
    while ~isempty(li)
        lim = str2num(li);
    axesHandles = findall(hw, 'Type', 'axes');  % 获取当前 figure 中所有的轴
    set(axesHandles, 'YLim', [-lim lim]);
    prompt="调整ylim:";
    li = input(prompt,'s');
    end
end
prompt="继续ica(回车下一个，1继续，0重置):";
in = input(prompt);
re_ica=in;
end

if ~isempty(clean_lfp)
for ii = 1:numel(vars_LFP)
    eval([vars_LFP{ii},'(:,trr)=clean_lfp(:,ii);']);
end

saveas(h1,fullfile(savepath,['tr',num2str(trr),'comp.jpg']));
saveas(h2,fullfile(savepath,['tr',num2str(trr),'comp_s.jpg']));
saveas(af1,fullfile(savepath,['tr',num2str(trr),'clean.jpg']));
saveas(af2,fullfile(savepath,['tr',num2str(trr),'clean_s.jpg']));
saveas(c1,fullfile(savepath,['tr',num2str(trr),'car.jpg']));
saveas(c2,fullfile(savepath,['tr',num2str(trr),'car.jpg']));
end
disp(['进度：',num2str(trr),'/120']);
end
