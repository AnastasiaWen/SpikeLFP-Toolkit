function lfp_extract_only1(trialpath,opn)
% used to aquaire one lfp signal
cd(trialpath);
%%
%clear;
nex = actxserver('NeuroExplorer.Application');
folders=pwd;
files = dir(fullfile(folders, '*.pl2'));
doc = nex.OpenDocument(fullfile(folders,files.name));

neuroncount=doc.NeuronCount;


EVT02=doc.Variable('EVT02');
EVT02=EVT02.Timestamps();

%only WB01
WBnames={'WB01'};
%% lfp extract
mkdir('Signal_Plot')
folderpath=pwd;
cd Signal_Plot\

Fs=1000;
for i=1:length(WBnames)
    WBnex=doc.Variable(WBnames{i});
    WB=WBnex.ContinuousValues();
    lfpname=['LFP',WBnames{i}(3:end)];
    if exist('opn')
        lfpdata=lfpextract(WB,WBnames{i},opn);
    else
        lfpdata=lfpextract(WB,WBnames{i}); 
    end
    eval([lfpname,'=lfpdata;']);
    eval(['before_',lfpname,'=lfpdata;']);
    disp(['Extracting ',lfpname,' finished']);
    close all
end

cd(folderpath);
close all
vars_before_lfp=who('-regexp', '^before_LFP');
%% raster lfp
vars_lfp=who('-regexp', '^LFP');
win=[1 1.5]; %event triggerd window (related to event02 as 0)
for i=1:length(WBnames)
    eval(['tmplfp=',vars_lfp{i},';'])
    trialslfp=createdatamatc(tmplfp',EVT02,Fs,win);
    eval([vars_lfp{i},'=trialslfp;']);
    disp(['Cut trials ',vars_lfp{i},' finished']);
end


%% save data

% 存储以 "LFP" 开头的变量
for i = 1:numel(vars_lfp)
    varName = vars_lfp{i};
    save(fullfile(trialpath,'data.mat'), varName, '-append');  % 将变量 LFP01 添加到 data.mat 文件中
end


% 存储以 "before_LFP" 开头的变量
for i = 1:numel(vars_before_lfp)
    varName = vars_before_lfp{i};
    save(fullfile(trialpath,'data.mat'), varName, '-append');  % 将变量 LFP01 添加到 data.mat 文件中
end

disp('Save data finished!');

%%
doc.Close();

end