
%%
clear;
clc;

nowpath=pwd;

weekpath='G:\2022_Mn_14365\validData\week03';
cd(weekpath);   

sitelist = dir(weekpath);
sitelist =sitelist ([sitelist .isdir]);


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
    for j=3:numel(triallist)-3
        trialname=triallist(j).name;
        if ~strcmp(trialname,'waveform')
            disp(['Processing...',trialname]);
            trialpath=fullfile(sitepath,trialname);
            % LFP_bad_channel(trialpath);
            LFP_sd_trial(trialpath);
            % lfp_carspec(trialpath,'data_sd.mat');
            %lfp_extract_only1(trialpath);
            %spk_lfp_extract(trialpath);
            % waveform_extract(trialpath); 
            
        end
    end
end

cd(weekpath);   

%system('shutdown.exe -s -t 14400')

% this script is used to batch process .pl2 in a week
% to aquaire data.mat(lfp&spk) and waveform ( need the
%  nex .py script to coorperate)
cd(weekpath);  


%%
clear;
clc;

nowpath=pwd;

weekpath='G:\2022_Mn_14365\validData\week04';
cd(weekpath);   

sitelist = dir(weekpath);
sitelist =sitelist ([sitelist .isdir]);


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
    for j=3:numel(triallist)-3
        trialname=triallist(j).name;
        if ~strcmp(trialname,'waveform')
            disp(['Processing...',trialname]);
            trialpath=fullfile(sitepath,trialname);
            % LFP_bad_channel(trialpath);
            LFP_sd_trial(trialpath);
            % lfp_carspec(trialpath,'data_sd.mat');
            %lfp_extract_only1(trialpath);
            %spk_lfp_extract(trialpath);
            % waveform_extract(trialpath); 
            
        end
    end
end

cd(weekpath);   

%system('shutdown.exe -s -t 14400')

% this script is used to batch process .pl2 in a week
% to aquaire data.mat(lfp&spk) and waveform ( need the
%  nex .py script to coorperate)
cd(weekpath);  


%%
clear;
clc;

nowpath=pwd;

weekpath='G:\2022_Mn_14365\validData\week06';
cd(weekpath);   

sitelist = dir(weekpath);
sitelist =sitelist ([sitelist .isdir]);


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
    for j=3:numel(triallist)-3
        trialname=triallist(j).name;
        if ~strcmp(trialname,'waveform')
            disp(['Processing...',trialname]);
            trialpath=fullfile(sitepath,trialname);
            % LFP_bad_channel(trialpath);
            LFP_sd_trial(trialpath);
            % lfp_carspec(trialpath,'data_sd.mat');
            %lfp_extract_only1(trialpath);
            %spk_lfp_extract(trialpath);
            % waveform_extract(trialpath); 
            
        end
    end
end

cd(weekpath);   

%system('shutdown.exe -s -t 14400')

% this script is used to batch process .pl2 in a week
% to aquaire data.mat(lfp&spk) and waveform ( need the
%  nex .py script to coorperate)
cd(weekpath);  


%%
clear;
clc;

nowpath=pwd;

weekpath='G:\2022_Mn_14365\validData\week07';
cd(weekpath);   

sitelist = dir(weekpath);
sitelist =sitelist ([sitelist .isdir]);


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
    for j=3:numel(triallist)-3
        trialname=triallist(j).name;
        if ~strcmp(trialname,'waveform')
            disp(['Processing...',trialname]);
            trialpath=fullfile(sitepath,trialname);
            % LFP_bad_channel(trialpath);
            LFP_sd_trial(trialpath);
            % lfp_carspec(trialpath,'data_sd.mat');
            %lfp_extract_only1(trialpath);
            %spk_lfp_extract(trialpath);
            % waveform_extract(trialpath); 
            
        end
    end
end

cd(weekpath);   

%system('shutdown.exe -s -t 14400')

% this script is used to batch process .pl2 in a week
% to aquaire data.mat(lfp&spk) and waveform ( need the
%  nex .py script to coorperate)
cd(weekpath);  

%%
clear;
clc;

nowpath=pwd;

weekpath='G:\2022_Mn_14365\validData\week08';
cd(weekpath);   

sitelist = dir(weekpath);
sitelist =sitelist ([sitelist .isdir]);


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
    for j=3:numel(triallist)-3
        trialname=triallist(j).name;
        if ~strcmp(trialname,'waveform')
            disp(['Processing...',trialname]);
            trialpath=fullfile(sitepath,trialname);
            % LFP_bad_channel(trialpath);
            LFP_sd_trial(trialpath);
            % lfp_carspec(trialpath,'data_sd.mat');
            %lfp_extract_only1(trialpath);
            %spk_lfp_extract(trialpath);
            % waveform_extract(trialpath); 
            
        end
    end
end

cd(weekpath);   

%system('shutdown.exe -s -t 14400')

% this script is used to batch process .pl2 in a week
% to aquaire data.mat(lfp&spk) and waveform ( need the
%  nex .py script to coorperate)
cd(weekpath);  


%%
clear;
clc;

nowpath=pwd;

weekpath='G:\2022_Mn_14365\validData\week10';
cd(weekpath);   

sitelist = dir(weekpath);
sitelist =sitelist ([sitelist .isdir]);


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
    for j=3:numel(triallist)-3
        trialname=triallist(j).name;
        if ~strcmp(trialname,'waveform')
            disp(['Processing...',trialname]);
            trialpath=fullfile(sitepath,trialname);
            % LFP_bad_channel(trialpath);
            LFP_sd_trial(trialpath);
            % lfp_carspec(trialpath,'data_sd.mat');
            %lfp_extract_only1(trialpath);
            %spk_lfp_extract(trialpath);
            % waveform_extract(trialpath); 
            
        end
    end
end

cd(weekpath);   

%system('shutdown.exe -s -t 14400')

% this script is used to batch process .pl2 in a week
% to aquaire data.mat(lfp&spk) and waveform ( need the
%  nex .py script to coorperate)
cd(weekpath);  


%%
clear;
clc;

nowpath=pwd;

weekpath='G:\2023_Mn_14325_RightBrain\validData\Week04';
cd(weekpath);   

sitelist = dir(weekpath);
sitelist =sitelist ([sitelist .isdir]);


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
    for j=3:numel(triallist)-3
        trialname=triallist(j).name;
        if ~strcmp(trialname,'waveform')
            disp(['Processing...',trialname]);
            trialpath=fullfile(sitepath,trialname);
            % LFP_bad_channel(trialpath);
            LFP_sd_trial(trialpath);
            % lfp_carspec(trialpath,'data_sd.mat');
            %lfp_extract_only1(trialpath);
            %spk_lfp_extract(trialpath);
            % waveform_extract(trialpath); 
            
        end
    end
end

cd(weekpath);   

%system('shutdown.exe -s -t 14400')

% this script is used to batch process .pl2 in a week
% to aquaire data.mat(lfp&spk) and waveform ( need the
%  nex .py script to coorperate)
cd(weekpath);  
%%
clear;
clc;

nowpath=pwd;

weekpath='G:\2023_Mn_14325_RightBrain\validData\Week05';
cd(weekpath);   

sitelist = dir(weekpath);
sitelist =sitelist ([sitelist .isdir]);


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
    for j=3:numel(triallist)-3
        trialname=triallist(j).name;
        if ~strcmp(trialname,'waveform')
            disp(['Processing...',trialname]);
            trialpath=fullfile(sitepath,trialname);
            % LFP_bad_channel(trialpath);
            LFP_sd_trial(trialpath);
            % lfp_carspec(trialpath,'data_sd.mat');
            %lfp_extract_only1(trialpath);
            %spk_lfp_extract(trialpath);
            % waveform_extract(trialpath); 
            
        end
    end
end

cd(weekpath);   

%system('shutdown.exe -s -t 14400')

% this script is used to batch process .pl2 in a week
% to aquaire data.mat(lfp&spk) and waveform ( need the
%  nex .py script to coorperate)
cd(weekpath);  


