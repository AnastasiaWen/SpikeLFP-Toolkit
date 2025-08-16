function maxnum=oi_sp_neuro_tune(sitepath)
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



diPath = fullfile(sitepath,'results\orientation_tun');
% rmdir(diPath,"s")
mkdir(diPath);

spPath = fullfile(sitepath,'results\speed_tun');
% rmdir(spPath,"s")
% mkdir(spPath);

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
Diresults=struct();
Dindex=struct();
Presults=struct();
% incPath=fullfile(diPath,'increase');
% mkdir(incPath);
% decPath=fullfile(diPath,'decrease');
% mkdir(decPath);

varNames = fieldnames(matData);
for h = 1:length(varNames)
    varName = varNames{h};
    if startsWith(varName, 'SPKC')
        totalspkc = matData.(varName);
        totalangle=matData.(['StimAngle',varName(end-2:end)]);
        totalfr=[];
        totalsd=[];
        for i=1:5
                rc=find(totalangle==uniqAngle(i));
                refr=[];
                if ~isempty(rc)
                    for k=1:length(rc)
                        spkcount_after=binspikes(totalspkc(rc(k)),Fs,[1.2,2]);
                        fr_after=sum(spkcount_after)/0.8;
                        spkcount_before=binspikes(totalspkc(rc(k)),Fs,[0.5,0.75]);
                        fr_before=sum(spkcount_before)/0.25;
                        %refr=[refr,(fr_after-fr_before)/(fr_before+eps)];
                        %refr=[refr,(fr_after)/(fr_before+eps)];
                        refr=[refr,fr_after];
                    end
                    totalfr=[totalfr,mean(refr)];
                    totalsd=[totalsd,std(refr) / sqrt(size(refr,2))];
                else
				%实际上已经不可能是0，我已经剔除掉了在一个方向上没有任何trial的神经元
                    totalfr=[totalfr,0];
                    totalsd=[totalsd,0];
                end
        end
        Diresults.(varName)=totalfr;
        
        %totalfr=zscore(totalfr);
        % totalfr = normalize(totalfr,"range");
        sin_component = sum(totalfr .* sin(2*deg2rad(uniqAngle)));
        cos_component = sum(totalfr .* cos(2*deg2rad(uniqAngle)));
        numerator = sqrt(sin_component^2 + cos_component^2);
        
        % 计算分母部分
        denominator = sum(totalfr);
        
        % 计算方向指数 (DI)
        DI = numerator / denominator;
        Dindex.(varName)=DI;
        % min_fr = min(totalfr);
        % max_fr = max(totalfr);
        % scaled_mean_fr = (totalfr- min_fr) / (max_fr - min_fr);
        % scaled_stderr_fr = totalsd / (max_fr - min_fr); % 注意标准误也要缩放
        
		 % 进行随机化测试
		 num_permutations=1000;
        DI_perm = zeros(num_permutations, 1);
        for p = 1:num_permutations
            perm_angle = uniqAngle(randperm(length(uniqAngle)));  % 随机打乱角度
            sin_component_perm = sum(totalfr .* sin(2*deg2rad(perm_angle)));
            cos_component_perm = sum(totalfr .* cos(2*deg2rad(perm_angle)));
            numerator_perm = sqrt(sin_component_perm^2 + cos_component_perm^2);
            denominator_perm = sum(totalfr);
            DI_perm(p) = numerator_perm / denominator_perm;  % 计算每次打乱后的 OI 值
        end
        
        % 计算 p 值（实际 OI 值与随机化分布的比较）
        p_value = mean(DI_perm >= DI);
		Presults.(varName)=p_value;
        
        % 打印结果
        fprintf('Neuron %s: OI = %.3f, p-value = %.3f\n', varName, DI, p_value);
		
		
        [maxv, maxindx] = max(totalfr);
        
        if numel(find(totalfr==maxv))>1
            maxum(6)=maxnum(6)+1;
        else
           maxnum(maxindx)=maxnum(maxindx)+1;
        end
        figure;
        
        errorbar(uniqAngle,totalfr,totalsd,'LineWidth',2);
        xlabel('angle(°)'); ylabel('FR(Hz)');
        xticks(uniqAngle);
        % ylim([0 1]);
        xlim([0 180]);
        title(['OI=',num2str(DI),',','p:',num2str(p_value)]);
        
        saveas(gcf, fullfile(diPath, ['OiTune_',varName,'.png']));

    end
end

save(fullfile(diPath,'totaldata.mat'),"Diresults");
save(fullfile(diPath,'Oidata.mat'),"Dindex");
save(fullfile(diPath,'Pdata.mat'),"Presults");
close all
