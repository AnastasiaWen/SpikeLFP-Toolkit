function spk_lfp_extract(trialpath,opn)
% this script is used to extract lfp&spk in a pl2 to generate a data.mat
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

%for 22 week03 no wb01
%WBnames={'WB02','WB03','WB04','WB05','WB06','WB07','WB08','WB09','WB10','WB11','WB12','WB13','WB14','WB15','WB16'};
WBnames={'WB01','WB02','WB03','WB04','WB05','WB06','WB07','WB08','WB09','WB10','WB11','WB12','WB13','WB14','WB15','WB16'};
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

%% raster spike
win=[1 1.5];
for i=1:neuroncount
    neuron=doc.Neuron(i);
    neuronts=neuron.Timestamps();
    neuronname=neuron.Name;
    if length(neuronname)~=11
        continue
    end
    trialsspk=createdatamatpt(neuronts,EVT02,win)
    eval([neuronname(5:end),'=trialsspk;']);
    disp(['Cut trials ',neuronname(5:end),' finished']);
end
vars_spk=who('-regexp', '^SPKC');

%% save data
savedata = struct;

% 存储以 "LFP" 开头的变量
for i = 1:numel(vars_lfp)
    varName = vars_lfp{i};
    savedata.(varName) = eval(varName); % 将变量存入结构体中
end

% 存储以 "SPKC" 开头的变量
for i = 1:numel(vars_spk)
    varName = vars_spk{i};
    savedata.(varName) = eval(varName); % 将变量存入结构体中
end

% 存储以 "before_LFP" 开头的变量
for i = 1:numel(vars_before_lfp)
    varName = vars_before_lfp{i};
    savedata.(varName) = eval(varName); % 将变量存入结构体中
end

save('data.mat', '-struct', 'savedata');
disp('Save data finished!');

%%
%doc.Close();

end