% 该脚本用于删除实验中多余的EVT02 (均为120)

clear;clc;
nex = actxserver('NeuroExplorer.Application');
folders=pwd;
files = dir(fullfile(folders, '*.pl2'));
doc = nex.OpenDocument(fullfile(folders,files.name));
%doc = nex.ActiveDocument;

neuroncount=doc.NeuronCount;

EVT02=doc.Variable('EVT02');
EVT02=EVT02.Timestamps();

detrial=[];
cleaned_EVT02 = EVT02(1); % 保留第一个元素
for i=2:length(EVT02)
    if EVT02(i)-cleaned_EVT02<2.5
        detrial=[detrial,i];
    elseif EVT02(i)-cleaned_EVT02>5
            if  i==length(EVT02)
                detrial=[detrial,i];
            elseif i==2
                detrial=[detrial,1];
            elseif EVT02(i+1)-EVT02(i)>=2.5 && EVT02(i+1)-EVT02(i)<=5 && i+1~=length(EVT02)
                cleaned_EVT02=EVT02(i);
            else
                detrial=[detrial,i];
            end
    else
        cleaned_EVT02=EVT02(i);
    end
end
clearvars -except detrial EVT02

%%

load data.mat
vars_lfp=who('-regexp', '^LFP');
vars_spk=who('-regexp', '^SPKC');
vars_blfp=who('-regexp', '^before');
for i=1:length(vars_lfp)
    eval([vars_lfp{i},'(:,detrial)=[];'])
end
for i=1:length(vars_spk)
    eval([vars_spk{i},'(detrial)=[];'])
end

savevars=[vars_spk;vars_lfp,;vars_blfp];
save('data.mat',savevars{:});