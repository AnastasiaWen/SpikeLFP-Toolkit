function p_neuron(sitepath)

%sitepath='G:\2022_Mn_14365\validData\week08\D2_site02';
% coorperate with lfp_wavelet_spectro.m
close all;
combinepath=fullfile(sitepath,'combined');
cd(combinepath);
load p_data_detail.mat

increaseFolder = fullfile(sitepath,'results\spectro\multitaper\cut/increase');
mkdir(increaseFolder);
decreaseFolder = fullfile(sitepath,'results\spectro\multitaper\cut/decrease');
mkdir(decreaseFolder);
sourceFolder = fullfile(sitepath,'results\spectro\multitaper\cut');

%%

if ~isempty(increaseName)
    for k=1:numel(increaseName)
        tmpname=increaseName{k};
        sourceFile = fullfile(sourceFolder, ['total_LFP',tmpname(end-2:end),'.png']);
        destinationFile = fullfile(increaseFolder, ['total_LFP',tmpname(end-2:end),'.png']);
        movefile(sourceFile, destinationFile);
        %disp(['0']);
        
    end
end
if ~isempty(decreaseName)
    for k=1:numel(decreaseName)
        tmpname=decreaseName{k};
        sourceFile = fullfile(sourceFolder, ['total_LFP',tmpname(end-2:end),'.png']);
        destinationFile = fullfile(decreaseFolder, ['total_LFP',tmpname(end-2:end),'.png']);
        movefile(sourceFile, destinationFile);
        %disp(['0']);
        
    end
end

%%
increaseFolder = fullfile(sitepath,'results\spectro\multitaper\nocut/increase');
mkdir(increaseFolder);
decreaseFolder = fullfile(sitepath,'results\spectro\multitaper\nocut/decrease');
mkdir(decreaseFolder);
sourceFolder = fullfile(sitepath,'results\spectro\multitaper\nocut');

%%

if ~isempty(increaseName)
    for k=1:numel(increaseName)
        tmpname=increaseName{k};
        sourceFile = fullfile(sourceFolder, ['total_LFP',tmpname(end-2:end),'.png']);
        destinationFile = fullfile(increaseFolder, ['total_LFP',tmpname(end-2:end),'.png']);
        movefile(sourceFile, destinationFile);
        %disp(['0']);
        
    end
end
if ~isempty(decreaseName)
    for k=1:numel(decreaseName)
        tmpname=decreaseName{k};
        sourceFile = fullfile(sourceFolder, ['total_LFP',tmpname(end-2:end),'.png']);
        destinationFile = fullfile(decreaseFolder, ['total_LFP',tmpname(end-2:end),'.png']);
        movefile(sourceFile, destinationFile);
        %disp(['0']);
        
    end
end

