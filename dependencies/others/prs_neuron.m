clear;clc;
hist = readtable('spk_hist.csv');
time=hist{:,1};
names = hist.Properties.VariableNames;% 获取表格中所有的列名

%% 显著性检验neuron放电
% numneuron=1;
% totalp=[];
% pindx=time>=0.12 & time <=1;
% nindx=time<=0.12 & time >=-0.76;
% for i = 2:size(hist,2)
%         column_data = hist{:, i};
%         pdata=column_data(pindx);
%         ndata=column_data(nindx);
%         [h, p, ci, stats] = ttest(ndata, pdata);
%         t=stats.tstat;
%         totalp=[totalp,p];
%         if p<0.05
%             pneuron{numneuron}=names{i};
%             numneuron=numneuron+1;
%         end
% end

rindx=time>=0 & time<=1.5;
rtime=time(rindx);
binWidth=0.01;%bin
% 计算 bin 的中心
bin_centers = rtime(1:end) + binWidth / 2;

SAtemplate=zeros(length(rtime),1);
indx=rtime>=0.2 & rtime <=1;
SAtemplate(indx)=1;

numneuron=1;
for i = 2:size(hist,2)
        column_data = hist{:, i};
        column_data=normalize(column_data,"range");
        rdata=column_data(rindx);
        
        name=names{i};name=name(5:end);
        
        RAcoef=0;
        RAwin=0.01;
        for win=0.01:0.01:0.10
            RAtemplate=zeros(length(rtime),1);
            indx=rtime>=0.20 & rtime <=(0.20+win);
            RAtemplate(indx)=1;
            indx=rtime>=1 & rtime <=(1.0+win);
            RAtemplate(indx)=1;
            [r,lags]=xcorr(rdata, RAtemplate);
            r=max(r);
            if r>RAcoef
                RAcoef=r;
                RAwin=win;
            end
        end
           
        [r,lags]=xcorr(rdata, SAtemplate);
        SAcoef=max(r);
        totalcoef(numneuron,1)=RAcoef;
        totalcoef(numneuron,2)=SAcoef;
        if abs(RAcoef)<abs(SAcoef)
            pneuron{1,numneuron}=name;
            pneuron{2,numneuron}='SA';
            % 绘制直方图
            figure;
            hold on;
            bar(bin_centers, rdata, 'hist'); 
            plot(rtime,SAtemplate,'Color','red','LineWidth',2);
            % 添加标签和标题
            xlabel('Time');
            ylabel('Count');
            title(name);
            hold off;  
        else
            pneuron{1,numneuron}=name;
            pneuron{2,numneuron}='RA';
            % 绘制直方图
            RAtemplate=zeros(length(rtime),1);
            indx=rtime>=0.20 & rtime <=(0.20+RAwin);
            RAtemplate(indx)=1;
            indx=rtime>=1 & rtime <=(1.0+RAwin);
            RAtemplate(indx)=1;

            figure;
            hold on;
            bar(bin_centers, rdata, 'hist'); 
            plot(rtime,RAtemplate,'Color','red','LineWidth',2);
            % 添加标签和标题
            xlabel('Time');
            ylabel('Count');
            name=names(i);
            title(name);
            hold off; 
        end
        numneuron=numneuron+1;
end

