% %%d
% clear;
% clc;
% 
% nowpath=pwd;
% 
% workpath='G:\2023_Mn_14325_RightBrain\validData';
% cd(workpath);   
% 
% weeklist = dir(workpath);
% weeklist =weeklist ([weeklist .isdir]);
% 
% totalspkc=[];
% totalangle=[];
% totalspeed=[];
% 
% for j=3:numel(weeklist)
%     weekname=weeklist(j).name;
%     weekpath=fullfile(workpath,weekname);
%     disp(['Pre-processing...',weekname]);
%     cd(weekpath);   
% 
%     sitelist = dir(weekpath);
%     sitelist =sitelist ([sitelist .isdir]);
% 
%     for i=3:numel(sitelist)
%         sitename=sitelist(i).name;
%         sitepath=fullfile(weekpath,sitename);
%         disp(['Pre-processing...',sitename]);
%         matData = load(fullfile(sitepath,'combined/Combined_cut.mat'));
%         varNames = fieldnames(matData);
%         for k = 1:length(varNames)
%             varName = varNames{k};
%             if startsWith(varName, 'SPKC')
%                 spkcData = matData.(varName);
%                 totalspkc= [totalspkc,spkcData];
%             elseif startsWith(varName, 'StimAngle')
%                 spkcData = matData.(varName);
%                 totalangle= [totalangle;spkcData];
%             elseif startsWith(varName, 'StimSpeed')
%                 spkcData = matData.(varName);
%                 totalspeed= [totalspeed;spkcData];
%             end
%         end
%     end
% end
% disp(['combine spkc Finished!']);
% cd(workpath);   
% 
% clearvars -except totalspkc totalangle totalspeed


%%
%uniqSpeed=[5,10,20,25,40,50,80,100,160,200,240,320];
Fs=1000;
uniqAngle=[0,45,90,135,180];% 0,45,90,135,180 per 24times

%%
totalresults=struct();
totalfr=[];
for i=1:5
        rc=find(totalangle==uniqAngle(i));
        refr=[];
        for k=1:length(rc)
            spkcount_after=binspikes(totalspkc(rc(k)),Fs,[1.2,2]);
            fr_after=sum(spkcount_after)/0.8;
            spkcount_before=binspikes(totalspkc(rc(k)),Fs,[0.5,0.75]);
            fr_before=sum(spkcount_before)/0.25;
            %refr=[refr,(fr_after-fr_before)/(fr_before+eps)];
            refr=[refr,(fr_after)/(fr_before+eps)];
        end
        totalfr=[totalfr,mean(refr)];
end
totalfr=zscore(totalfr);
figure;plot(uniqAngle,(totalfr));
xlabel('angle/°'); ylabel('FR(z-score)');
%save(fullfile(savePath,'results.mat'),'-struct',"totalresults");


