function p_neuron_inc(sitepath)

%sitepath='G:\2022_Mn_14365\validData\week08\D2_site02';
% coorperate with lfp_wavelet_spectro.m
close all;
combinepath=fullfile(sitepath,'combined');
cd(combinepath);
load p_data_cut.mat
load p_data_detail.mat
if isfile(fullfile('p_data_cut_inc.mat'))
    delete('p_data_cut_inc.mat');
end
% destinationFolder = fullfile(sitepath,'results/spk_feature/nosig');
% if isfolder(destinationFolder)
%     rmdir(destinationFolder,'s');
%     mkdir(destinationFolder);
% else
%     mkdir(destinationFolder);
% end

%%
Fs=1000;

vars_lfp = who('-regexp', '^LFP');
vars_spk = increaseName';


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
    vars_sp = increaseName(~cellfun('isempty', regexp(increaseName, tmpname)));
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
vars_orderS = cellfun(@(s) strrep(s, 'SPKC', 'Order'), vars_spk, 'UniformOutput', false);
vars_orderL = cellfun(@(s) strrep(s, 'LFP', 'Order'), vars_lfp, 'UniformOutput', false);
vars_speed = cellfun(@(s) strrep(s, 'SPKC', 'StimSpeed'), vars_spk, 'UniformOutput', false);
vars_angle = cellfun(@(s) strrep(s, 'SPKC', 'StimAngle'), vars_spk, 'UniformOutput', false);
vars_wave =cellfun(@(s) strrep(s, 'SPKC', 'Waveform'), vars_spk, 'UniformOutput', false);
vars_wavetime = cellfun(@(s) strrep(s, 'SPKC', 'Wavetime'), vars_spk, 'UniformOutput', false);

savename=[vars_lfp;vars_spk;vars_angle;vars_speed;vars_wave;vars_wavetime;vars_orderS;vars_orderL];
save(fullfile(combinepath,'p_data_cut_inc.mat'),savename{:});


