function anova2sa(sitepath)
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
    matData=load('p_data_cut_inc.mat');
% matData = load(fullfile(sitepath,'combined/combined_.mat'));
% load(fullfile(sitepath,'combined/p_data_detail.mat'));
% load(fullfile(sitepath,'results\direction_tun\totaldata.mat'));
% load(fullfile(sitepath,'results\speed_tun\totaldata.mat'));



diPath = fullfile(sitepath,'results\anova');
rmdir(diPath,"s")
mkdir(diPath);


%%
uniqSpeed=[5,10,20,25,40,50,80,100,160,200,240,320];
Fs=1000;
uniqAngle=[0,45,90,135,180];% 0,45,90,135,180 per 24times

maxnum=zeros(6,1);
neuron0fr={};
neuron45fr={};
neuron90fr={};
neuron135fr={};
neuron180fr={};
neuronmix={};
neurondetail={neuron0fr,neuron45fr,neuron90fr,neuron135fr,neuron180fr,neuronmix};
%%
anovaresults=struct();
Dindex=struct();
Presults=struct();
% incPath=fullfile(diPath,'increase');
% mkdir(incPath);
% decPath=fullfile(diPath,'decrease');
% mkdir(decPath);

varNames = fieldnames(matData);
for h = 1:length(varNames)
    flag=true;
    varName = varNames{h};
    if startsWith(varName, 'SPKC')
        totalspkc = matData.(varName);
        totalangle=matData.(['StimAngle',varName(end-2:end)]);
        totalspeed=matData.(['StimSpeed',varName(end-2:end)]);
        totalfr=cell(5,12);
        totalsd=cell(5,12);
        for i=1:5
            for j=1:12
                rc=find(totalangle==uniqAngle(i)&totalspeed==uniqSpeed(j));
                refr=[];
                if numel(rc)>=2
                    for k=1:length(rc)
                        spkcount_after=binspikes(totalspkc(rc(k)),Fs,[1.2,2]);
                        fr_after=sum(spkcount_after)/0.8;
                        % spkcount_before=binspikes(totalspkc(rc(k)),Fs,[0.5,0.75]);
                        % fr_before=sum(spkcount_before)/0.25;
                        %refr=[refr,(fr_after-fr_before)/(fr_before+eps)];
                        %refr=[refr,(fr_after)/(fr_before+eps)];
                        refr=[refr,fr_after];
                    end
                    totalfr(i,j)={refr};
                    % totalsd(i,j)={std(refr) / sqrt(size(refr,2));
                else
                    flag=false;
                    break;
				%实际上已经不可能是0，我已经剔除掉了在一个方向上没有任何trial的神经元
 % totalfr(i,j)={0};
                    % totalsd(i,j)=0;
                end
                if flag==false
                    break;
                end
            end
            if flag==false
                    break;
                end
        end
        if flag==false
                    continue;
        end
        p=anova2_single_neuro(totalfr);
      if any(p<0.05)
        anovaresults.(varName)=p;
        disp(sitepath);
      end

    end
end
save(fullfile(diPath,'totaldata.mat'),"anovaresults");
close all
