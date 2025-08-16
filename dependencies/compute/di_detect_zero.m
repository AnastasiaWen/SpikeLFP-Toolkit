function di_detect_zero(sitepath)
% coorperate with lfp_wavelet_spectro.m
% close all;
% combinepath=fullfile(sitepath,'combined');
% cd(combinepath);
% load Combined.mat
% 
% % savePath = fullfile(sitepath,'results\spectro\wavelet');
% % rmdir(savePath,'s'); % 选择
% savePath = fullfile(sitepath,'results\spectro\wavelet\nocut');
% mkdir(savePath);

% coorperate with lfp_wavelet_spectro.m
close all;
diPath = fullfile(sitepath,'results\direction_tun');

load(fullfile(diPath,'totaldata.mat'))
% 使用 structfun 来检查每个字段是否包含 0
has_zero = structfun(@(x) any(x(:) == 0), Diresults);

% 找到包含 0 值的字段名称
fields_with_zero = fieldnames(Diresults);
fields_with_zero = fields_with_zero(has_zero);

% 计算每个字段中 0 值的数量
zero_counts = structfun(@(x) sum(x(:) == 0), Diresults);

% 仅保留包含 0 值的字段和对应的数量
zero_counts = zero_counts(has_zero);

% 显示包含 0 值的字段名称及其数量
disp('Fields containing zero values and their counts:');
for i = 1:length(fields_with_zero)
    fprintf('%s: %d zeros\n', fields_with_zero{i}, zero_counts(i));
end


% 使用 structfun 来检查每个字段是否包含 0
has_zero = structfun(@(x) any(isnan(x(:) )), Diresults);

% 找到包含 0 值的字段名称
fields_with_zero = fieldnames(Diresults);
fields_with_zero = fields_with_zero(has_zero);

% 计算每个字段中 0 值的数量
zero_counts = structfun(@(x) sum(isnan(x(:) )), Diresults);

% 仅保留包含 0 值的字段和对应的数量
zero_counts = zero_counts(has_zero);

% 显示包含 0 值的字段名称及其数量
disp('Fields containing NaN values and their counts:');
for i = 1:length(fields_with_zero)
    fprintf('%s: %d zeros\n', fields_with_zero{i}, zero_counts(i));
end