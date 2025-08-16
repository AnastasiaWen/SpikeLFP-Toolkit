function [onum,to,in,de]=p_neuron(sitepath)

%sitepath='G:\2022_Mn_14365\validData\week08\D2_site02';
% coorperate with lfp_wavelet_spectro.m
close all;
combinepath=fullfile(sitepath,'combined');
cd(combinepath);
load Combined_cut_new.mat
if isfile(fullfile('p_data.mat'))
    delete('p_data.mat');
    delete('p_data_cut.mat');
    delete('p_data_detail.mat');
end
increaseFolder = fullfile(sitepath,'results/spk_feature/increase');
mkdir(increaseFolder);
decreaseFolder = fullfile(sitepath,'results/spk_feature/decrease');
mkdir(decreaseFolder);
sourceFolder = fullfile(sitepath,'results/spk_feature');

increaseName={};
decreaseName={};
% destinationFolder = fullfile(sitepath,'results/spk_feature/nosig');
% if isfolder(destinationFolder)
%     rmdir(destinationFolder,'s');
%     mkdir(destinationFolder);
% else
%     mkdir(destinationFolder);
% end

%%
Fs=1000;

vars = who;
vars_lfp = who('-regexp', '^LFP');
vars_spk = who('-regexp', '^SPKC');
vars_speed = who('-regexp', '^StimSpeed');
vars_angle = who('-regexp', '^StimAngle');
vars_wave = who('-regexp', '^Waveform');
vars_wavetime = who('-regexp', '^Wavetime');

num=0;
savedata=[];
decrease=0;
increase=0;
onum=numel(vars_spk);
for k=1:numel(vars_spk)
    
    %disp(vars_spk{k})
    eval(['spk=',vars_spk{k},';']);
if ~isempty(spk)
    spkc1=binspikes(spk,100,[0.5,0.9]); 
    fr1=sum(spkc1,2)/0.01;

    spkc2=binspikes(spk,100,[1,1.4]);
    fr2=sum(spkc2,2)/0.01;

    spkc3=binspikes(spk,100,[1.4,1.8]);
    fr3=sum(spkc3,2)/0.01;

    spkc4=binspikes(spk,100,[1.8,2.2]);
    fr4=sum(spkc4,2)/0.01;
    %fr2=fr2(1:end-1);
    spkc5=binspikes(spk,1/0.03,[1,2.2]);
    fr5=sum(spkc5,2)/0.03;


    [p2, h] = signrank(fr1, fr2);
    [p3, h] = signrank(fr1, fr3);
    [p4, h] = signrank(fr1, fr4);
    [p5, h] = signrank(fr1, fr5);
    if p5<0.05
        num=num+1;
        p=5;
        savedata=[savedata,k];
        if mean(fr5)-mean(fr1)>0
            % sourceFile = fullfile(sourceFolder, [vars_spk{k},'.png']);
            % destinationFile = fullfile(increaseFolder, [vars_spk{k},'.png']);
            % movefile(sourceFile, destinationFile);
            increaseName=[increaseName,vars_spk{k}];
            increase=increase+1;
        elseif mean(fr5)-mean(fr1)<0
            % sourceFile = fullfile(sourceFolder, [vars_spk{k},'.png']);
            % destinationFile = fullfile(decreaseFolder, [vars_spk{k},'.png']);
            % movefile(sourceFile, destinationFile);
            decreaseName=[decreaseName,vars_spk{k}];
            decrease=decrease+1;
        end
      
    elseif p4<0.05
        num=num+1;
        p=4;
        savedata=[savedata,k];
        if mean(fr4)-mean(fr1)>0
            % sourceFile = fullfile(sourceFolder, [vars_spk{k},'.png']);
            % destinationFile = fullfile(increaseFolder, [vars_spk{k},'.png']);
            % movefile(sourceFile, destinationFile);
            increaseName=[increaseName,vars_spk{k}];
            increase=increase+1;
        elseif mean(fr4)-mean(fr1)<0
            % sourceFile = fullfile(sourceFolder, [vars_spk{k},'.png']);
            % destinationFile = fullfile(decreaseFolder, [vars_spk{k},'.png']);
            % movefile(sourceFile, destinationFile);
            decreaseName=[decreaseName,vars_spk{k}];
            decrease=decrease+1;
        end

    elseif p3<0.05
        num=num+1;
        p=3;
        savedata=[savedata,k];
        if mean(fr3)-mean(fr1)>0
            % sourceFile = fullfile(sourceFolder, [vars_spk{k},'.png']);
            % destinationFile = fullfile(increaseFolder, [vars_spk{k},'.png']);
            % movefile(sourceFile, destinationFile);
            increaseName=[increaseName,vars_spk{k}];
            increase=increase+1;
        elseif mean(fr3)-mean(fr1)<0
            % sourceFile = fullfile(sourceFolder, [vars_spk{k},'.png']);
            % destinationFile = fullfile(decreaseFolder, [vars_spk{k},'.png']);
            % movefile(sourceFile, destinationFile);
            decreaseName=[decreaseName,vars_spk{k}];
             decrease=decrease+1;
        end
        %disp(['1']);
    elseif  p2<0.05
        num=num+1;
        p=2;
        savedata=[savedata,k];
        if mean(fr2)-mean(fr1)>0
            % sourceFile = fullfile(sourceFolder, [vars_spk{k},'.png']);
            % destinationFile = fullfile(increaseFolder, [vars_spk{k},'.png']);
            % movefile(sourceFile, destinationFile);
            increaseName=[increaseName,vars_spk{k}];
            increase=increase+1;
        elseif mean(fr2)-mean(fr1)<0
            % sourceFile = fullfile(sourceFolder, [vars_spk{k},'.png']);
            % destinationFile = fullfile(decreaseFolder, [vars_spk{k},'.png']);
            % movefile(sourceFile, destinationFile);
            decreaseName=[decreaseName,vars_spk{k}];
             decrease=decrease+1;
        end
        %disp(['1']);
        %disp(['1']);
    else
        % sourceFile = fullfile(sourceFolder, [vars_spk{k},'.png']);
        % destinationFile = fullfile(destinationFolder, [vars_spk{k},'.png']);
        % movefile(sourceFile, destinationFile);
        % %disp(['0']);
        
    end
end
end
disp(sitepath);
disp(['Total:',num2str(numel(vars_spk)),';','Sig:',num2str(num)]);
disp(['inc:',num2str(increase),';','dec:',num2str(decrease),';','totalsig:',num2str(increase+decrease)]);
in=increase;de=decrease;to=increase+decrease;

save(fullfile(combinepath,'p_data_detail.mat'),"increaseName","decreaseName");
%%
disp(['Sig:',num2str(num)]);

vars_spk=vars_spk(savedata);
% vars_lfp = vars_lfp(savedata);
vars_speed =vars_speed(savedata);
vars_angle = vars_angle(savedata);
vars_wave = vars_wave(savedata);
vars_wavetime = vars_wavetime (savedata);
% 使用 cellfun 和 strrep 替换 'SPKC' 为 'Order'
vars_orderS = cellfun(@(s) strrep(s, 'SPKC', 'Order'), vars_spk, 'UniformOutput', false);
vars_orderL = cellfun(@(s) strrep(s, 'LFP', 'Order'), vars_lfp, 'UniformOutput', false);

% if ~isempty(vars_spk)
savename=[vars_lfp;vars_spk;vars_angle;vars_speed;vars_wave;vars_wavetime;vars_orderS;vars_orderL];
save(fullfile(combinepath,'p_data.mat'),savename{:});

% load Combined.mat
% save(fullfile(combinepath,'p_data.mat'),savename{:});
% end

%%

dellfp=[];
for i = 1:numel(vars_lfp)
    deltrial=[];
    ldeltrial=[];
    eval([ 'tmplfp=',vars_lfp{i},';']);
    lfporder_name=vars_lfp{i};
    lfporder_name=['Order',lfporder_name(end-1:end)];
    lfporder=eval(lfporder_name);

    tmpname=vars_lfp{i};tmpname=tmpname(end-1:end);tmpname=['^SPKC',tmpname];
    vars_sp = vars_spk(~cellfun('isempty', regexp(vars_spk, tmpname)));
    if isempty(vars_sp)
        dellfp=[dellfp,i];
        continue;
    else

        sp_indx=cell(1,numel(vars_sp));
        sp_orders=[];
        for ss=1:numel(vars_sp)
            tmpname=vars_sp{ss};tmpname=tmpname(end-2:end);
            orderspk=eval(['Order',tmpname]);
            sp_orders=[sp_orders,num2cell(orderspk,1)];
        end
            
        combined = sp_orders;
        
        % 初始化一个 cell 用于存储 unique 的结果
        unique_elements = {};
        
        % 遍历 combined 数组，比较每个元素是否已经存在于 unique_elements 中
        for ii = 1:length(combined)
            current_element = combined{ii};
            is_unique = true;  % 标记当前元素是否唯一
            
            % 遍历 unique_elements，检查当前元素是否已经存在
            for j = 1:length(unique_elements)
                if isequal(current_element, unique_elements{j})
                    is_unique = false;  % 如果已经存在，将标记设置为 false
                    break;
                end
            end
            
            % 如果是唯一的元素，添加到 unique_elements 中
            if is_unique
                unique_elements{end+1} = current_element;
            end
        end
    
    
        a1 = num2cell(lfporder,1);
        a2 = unique_elements;
        
        % 初始化一个数组用于存储 a1 中不在 a2 中的元素的索引
        idx_not_in_a2 = [];
        
        % 遍历 a1 的每个元素，并检查是否存在于 a2 中
        for ii = 1:length(a1)
            current_element = a1{ii};
            is_in_a2 = false;  % 标记当前元素是否在 a2 中
            
            % 遍历 a2 并比较
            for j = 1:length(a2)
                if isequal(current_element, a2{j})
                    is_in_a2 = true;  % 如果存在，标记为 true
                    break;  % 退出内部循环
                end
            end
            
            % 如果当前元素不在 a2 中，记录它的索引
            if ~is_in_a2
                idx_not_in_a2 = [idx_not_in_a2, ii];
            end
        end
        
        eval([vars_lfp{i},'(:,idx_not_in_a2)=[];']);
        eval([lfporder_name,'(:,idx_not_in_a2)=[];']);
    end
end

vars_lfp(dellfp)=[];
vars_orderL = cellfun(@(s) strrep(s, 'LFP', 'Order'), vars_lfp, 'UniformOutput', false);

savename=[vars_lfp;vars_spk;vars_angle;vars_speed;vars_wave;vars_wavetime;vars_orderS;vars_orderL];
save(fullfile(combinepath,'p_data_cut.mat'),savename{:});


