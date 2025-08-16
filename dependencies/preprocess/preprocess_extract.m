% this script is used to batch process .pl2 in a week
% to aquaire data.mat(lfp&spk) and waveform ( need the
%  nex .py script to coorperate)

clear;
clc;

nowpath=pwd;

weekpath='G:\2022_Mn_14365\20221013_Expt_Week03\fromPlexon_RecordPL2';
cd(weekpath);   

sitelist = dir(weekpath);
sitelist =sitelist ([sitelist .isdir]);9999


for i=3:numel(sitelist)
    sitename=sitelist(i).name;
    sitepath=fullfile(weekpath,sitename);
    disp(['Pre-processing...',sitename]);
    triallist=dir(sitepath);
    % 筛选出文件夹名
    triallist=triallist([triallist.isdir]);

    if numel(triallist)==2
        continue;
    end

    mkdir(fullfile(sitepath,'waveform'));
    for j=3:numel(triallist)
        trialname=triallist(j).name;
        if ~strcmp(trialname,'waveform')
            disp(['Processing...',trialname]);
            trialpath=fullfile(sitepath,trialname);
            
            %lfp_extract_only1(trialpath);
            %spk_lfp_extract(trialpath);
            waveform_extract(trialpath); 
            
        end
    end
end

cd(weekpath);   

%system('shutdown.exe -s -t 14400')