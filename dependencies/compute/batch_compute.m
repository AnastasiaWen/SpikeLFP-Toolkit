% this script is used to batch process .pl2 in a week
% to aquaire data.mat(lfp&spk) and waveform ( need the
%  nex .py script to coorperate)

clear;
clc;

nowpath=pwd;

workpath='G:\2023_Mn_14325_RightBrain\validData';
cd(workpath);   

weeklist = dir(workpath);
weeklist =weeklist ([weeklist .isdir]);

dir_tun_num={};
sitenum=0;
for j=3:numel(weeklist)
    weekname=weeklist(j).name;
    weekpath=fullfile(workpath,weekname);
    disp(['Pre-processing...',weekname]);
    cd(weekpath);   
    dir_tun_num(sitenum+1,1)={weekname};
    sitelist = dir(weekpath);
    sitelist =sitelist ([sitelist .isdir]);
    
    for i=3:numel(sitelist)
        sitename=sitelist(i).name;
        sitepath=fullfile(weekpath,sitename);
        disp(['Pre-processing...',sitename]);
        % combinedetail_order(sitepath);
        % CutTrial_data(sitepath) 
        % Combined_cut_ne(sitepath) 
        % % % 
        %   [onum,to,in,de]=p_neuron(sitepath);
        %  sitenum=sitenum+1;
        % dir_tun_num(sitenum,3:6)=num2cell([onum,to,in,de]);
        % dir_tun_num(sitenum,2)={sitename};
        % p_neuron_inc(sitepath)
         % lfp_spectro(sitepath);
         % GC_Other_channel(sitepath);
        % FTA_rayleigh_test(sitepath);
        % FTA_site(sitepath);
        % GC_site(sitepath);
        % % 
        % spk_feature(sitepath)
        % sitenum=sitenum+1;
        % anova2sa(sitepath)
        oi4_sp_neuro_tune(sitepath)
        % dir_tun_num(sitenum,3:8)=num2cell(di_sp_neuro_tune(sitepath));
        % dir_tun_num(sitenum,2)={sitename};
        % clean_up_GCre(sitepath)
        % clean_up_Core(sitepath)

        % 
        % sitenum=sitenum+1;
        % dir_tun_num(sitenum,3)=num2cell(neuro_num(sitepath));
        % dir_tun_num(sitenum,2)={sitename};

% lfp_sta_detail(sitepath);

        % LFP_sd_delete22(sitepath);

        %di_detect_zero(sitepath);
        % di_sp_neuro_tune(sitepath)
        %p_neuron_classify(sitepath)
         %di_sp_neuro_tune(sitepath);
        %channel_lfp_depth(sitepath);
         
        % if exist(fullfile(sitepath,'results/spectro/wavelet/nocut'), 'dir') ~= 7
          % lfp_spectro(sitepath);
        % end
         % CutTrial_data(sitepath) 
         % lfp_spectro(sitepath);
         % spk_feature(sitepath);
        %channel_lfp(sitepath);
        %channel_cohe(sitepath);
        % combine_speed_angle(sitepath);
        % speedWithangle(sitepath);
    
    end
end
disp(['Finished!']);
cd(workpath);   
% 
% filename = 'neuro_oi.xlsx';
% % filename = 'neuro_num_p.xlsx';
% xlswrite(filename, dir_tun_num);
% system('shutdown.exe -s -t 60')
%%
clear;
clc;


nowpath=pwd;

workpath='G:\2022_Mn_14365\validData';
cd(workpath);   

weeklist = dir(workpath);
weeklist =weeklist ([weeklist .isdir]);

dir_tun_num={};
sitenum=0;
for j=3:numel(weeklist)
    weekname=weeklist(j).name;
    weekpath=fullfile(workpath,weekname);
    disp(['Pre-processing...',weekname]);
    cd(weekpath);   
    dir_tun_num(sitenum+1,1)={weekname};
    sitelist = dir(weekpath);
    sitelist =sitelist ([sitelist .isdir]);
    
    for i=3:numel(sitelist)
        sitename=sitelist(i).name;
        sitepath=fullfile(weekpath,sitename);
        disp(['Pre-processing...',sitename]);
        % GC_Other_channel(sitepath)
        % combinedetail_order(sitepath);
        %   CutTrial_data(sitepath) 
        % Combined_cut_ne(sitepath) 
        % % 
        %   [onum,to,in,de]=p_neuron(sitepath);
        %  sitenum=sitenum+1;
        % dir_tun_num(sitenum,3:6)=num2cell([onum,to,in,de]);
        % dir_tun_num(sitenum,2)={sitename};
        % p_neuron_inc(sitepath)
         % lfp_spectro(sitepath);
         % GC_Other_channel(sitepath);
        % FTA_site(sitepath);
        % FTA_rayleigh_test(sitepath);
        % GC_site(sitepath);
        % sitenum=sitenum+1;
        % dir_tun_num(sitenum,3)=num2cell(neuro_num(sitepath));
        % dir_tun_num(sitenum,2)={sitename};
        % 
% lfp_sta_detail(sitepath);

        % LFP_sd_delete22(sitepath);
        % lfp_spectro(sitepath);
        %di_detect_zero(sitepath);
        %di_sp_neuro_tune(sitepath)
        %p_neuron_classify(sitepath)
         % di_sp_neuro_tune(sitepath);
        %channel_lfp_depth(sitepath);
        %          spk_feature(sitepath)
        % sitenum=sitenum+1;
        % dir_tun_num(sitenum,3:8)=num2cell(di_sp_neuro_tune(sitepath));
        % dir_tun_num(sitenum,2)={sitename};
        % clean_up_GCre(sitepath)
        % if exist(fullfile(sitepath,'results/spectro/wavelet/nocut'), 'dir') ~= 7
          %lfp_spectro(sitepath);
        % end
         % CutTrial_data(sitepath) 
         % lfp_spectro(sitepath);
         % spk_feature(sitepath);
        % channel_lfp(sitepath);
        %channel_cohe(sitepath);
        % combine_speed_angle(sitepath);
        % speedWithangle(sitepath);
        % clean_up_Core(sitepath)
        oi4_sp_neuro_tune(sitepath)
        % anova2sa(sitepath)
    
    end
end
disp(['Finished!']);
cd(workpath);    

% filename = 'neuro_oi.xlsx';
% % filename = 'neuro_num_p.xlsx';
% xlswrite(filename, dir_tun_num);