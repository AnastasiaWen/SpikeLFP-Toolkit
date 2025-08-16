function num=neuro_num(sitepath)
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
combinepath=fullfile(sitepath,'combined');
cd(combinepath);
if isfile('Combined_cut_new.mat')
    load('Combined_cut_new.mat');
else
    load('Combined_cut.mat');
end
% matData = load(fullfile(sitepath,'combined/combined_.mat'));
% load(fullfile(sitepath,'combined/p_data_detail.mat'));
% load(fullfile(sitepath,'results\direction_tun\totaldata.mat'));
% load(fullfile(sitepath,'results\speed_tun\totaldata.mat'));

vars_spk=who('-regexp','^SPKC');
num=numel(vars_spk);
