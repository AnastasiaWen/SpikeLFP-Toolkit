function waveform_extract(trialpath)
%   this function is used to extract all waveform in 
%   sorted spikes and generate a .mat
cd(trialpath);
%%
%clear;
nex = actxserver('NeuroExplorer.Application');
folders=pwd;
[wfpath,trialname]=fileparts(folders);
wfpath=fullfile(wfpath,'waveform');
%mkdir(wfpath);
files = dir(fullfile(folders, '*.pl2'));
doc = nex.OpenDocument(fullfile(folders,files.name));
%doc = nex.ActiveDocument;

wavecount=doc.WaveCount;

for i=1:wavecount
    wave=doc.Wave(i);
    wfs=wave.WaveformValues();
    wavename=wave.Name();
    if length(wavename)~=14
        continue
    end
    wfname=[trialname(11:end),'_',wavename(9:end),'.mat'];
    save(fullfile(wfpath,wfname),"wfs");
end


doc.Close();

