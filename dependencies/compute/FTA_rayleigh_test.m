function FTA_rayleigh_test(sitepath)
% %coorperate with lfp_wavelet_spectro.m
% close all;
% combinepath=fullfile(sitepath,'combined');
% cd(combinepath);
% load p_data.mat
% 
% savePath = fullfile(sitepath,'results\spectro\wavelet');
% rmdir(savePath,'s'); % 选择
% savePath = fullfile(sitepath,'results\spectro\wavelet\nocut');
% mkdir(savePath);

% coorperate with lfp_wavelet_spectro.m
close all;
combinepath=fullfile(sitepath,'combined');
% cd(combinepath);
% load p_data_cut_sd.mat


savePath = fullfile(sitepath,'results\fta');
cd(savePath);
matData=load('totallfp_fta0to30.mat');



%%
Fs=1000;
pdata={};
num=0;
vars = who;
vars_lfp = who('-regexp', '^LFP');
vars_spk = who('-regexp', '^SPKC');

varnames=fieldnames(matData);
for i=1:numel(varnames)
    data=matData.(varnames{i});
    values=mean(data,1);
    
    % 假设你的相位数据为 angles，范围从 -π 到 π
    angles = linspace(-pi, pi, 801);
    % 将相位转换为单位向量
    unit_vectors = exp(1i * angles);
    % 计算加权的单位向量和
    R = sum(values .* unit_vectors) / sum(values);
    % 计算 Rayleigh 统计量
    R_stat = abs(R) * sqrt(sum(values));
    % 计算对应的 P 值
    p_value = exp(-R_stat^2);
    % 判断显著性
    alpha = 0.05; % 显著性水平
    if p_value < alpha
        disp('Rayleigh test: Significant unimodal distribution (P < 0.05)');
        num=num+1;
        pdata(num)={varnames{i}};
        figure;
        polarplot(angles, values, '-o');
        title('Field Triggered Average Phase Distribution');
        saveas(gcf, fullfile(savePath, ['p_',varnames{i},'.png']));
    else
        % disp('Rayleigh test: Not significant (P >= 0.05)');
    end
end
% 可视化分布

    save(fullfile(savePath,'pname'),'pdata');
    % save(fullfile(savePath,'lowlfp_spec.mat'),'-struct','lowlfp_spec');
