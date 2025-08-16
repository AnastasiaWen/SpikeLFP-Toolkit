Fs = 1000;                  %Sampling frequency
TW = 3;						%Choose time-bandwidth product of 3.
ntapers = 2*TW-1;			%Choose the # of tapers.
fpass = [0 100];
pad = 2;

params.Fs = Fs;				%... sampling frequency
params.tapers=[TW,ntapers];	%... time-bandwidth product & tapers.
params.fpass=fpass;       %... fpass
params.pad = pad;             %... padding.
params.trialave = 0;		%... trial average.
%%

% 初始化输出数据结构
combinedLFPData = {};      % 用于存储合并后的LFP矩阵
combinedSPKCData = {};     % 用于存储合并后的SPKC struct
combinedStimAngle = {};    % 用于存储StimAngle
combinedStimSpeed = {};    % 用于存储StimSpeed

% 获取所有LFP和SPK通道的变量名
LFPvars = who('LFP*');   % 所有以LFP开头的变量名
SPKCvars = who('SPKC*'); % 所有以SPKC开头的变量名

% 原始Order结构
totalOrder = [repmat(5, 1, 120), repmat(6, 1, 120), repmat(7, 1, 120), repmat(8, 1, 120);
              1:120, 1:120, 1:120, 1:120];  % 4个实验，每个实验120个trial

% 计数器，用于记录保留的有效trial数量
validTrialCount = 0;

% 遍历所有trial并进行匹配
for trialIdx = 1:size(totalOrder, 2)  % 遍历总Order中的所有trial
    LFPMatrix = [];  % 每个trial的合并LFP矩阵 (numLFPx2500)
    SPKCStruct = struct();  % 每个trial的合并SPKC struct
    trialStimAngle = [];    % 用于存储该trial的StimAngle
    trialStimSpeed = [];    % 用于存储该trial的StimSpeed
    
    % 获取当前trial的order值，根据totalOrder的顺序进行处理
    orderVal = totalOrder(:, trialIdx);  % 取出当前trial的Order

    % 初始化有效通道列表
    validChannels = {};  % 有效的通道列表
    
    % 遍历所有LFP通道，筛选出当前trial有数据的通道
    for i = 1:length(LFPvars)
        LFPvarName = LFPvars{i};
        LFPdata = eval(LFPvarName);  % 获取LFP数据
        LFPOrderVar = ['Order', LFPvarName(4:end)];  % 找到对应的Order变量
        LFPOrder = eval(LFPOrderVar);  % 获取LFP的Order值
        
        % 找到与当前order值匹配的trial索引
        matchIdx = find(ismember(LFPOrder', orderVal', 'rows'));  % 查找匹配的trial索引
        if ~isempty(matchIdx)
            % 如果有数据，将该通道加入有效通道列表
            validChannels = [validChannels, LFPvarName(4:end)];  % 记录通道名
        end
    end
    
    % 只保留有效通道数大于等于 4 的trial
    if length(validChannels) >= 4  % 至少保证LFP和SPKC通道各有2个
        % 随机选择一半的有效通道作为LFP通道
        numValidChannels = length(validChannels);
        numLFPSelected = floor(numValidChannels / 2);
        selectedLFPChannels = validChannels(randperm(numValidChannels, numLFPSelected));  % 选择的LFP通道名
        
        % 剩下的通道用于SPKC
        remainingChannels = setdiff(validChannels, selectedLFPChannels);
        selectedSPKCChannels = remainingChannels(randperm(length(remainingChannels), numLFPSelected));  % 选择的SPKC通道名
        
        % 遍历选中的LFP通道并收集数据
        for i = 1:length(selectedLFPChannels)
            LFPvarName = ['LFP', selectedLFPChannels{i}];  % 重新构造LFP变量名
            LFPdata = eval(LFPvarName);  % 获取对应的LFP数据矩阵 (2500xTrials)
            LFPOrderVar = ['Order', selectedLFPChannels{i}];  % 找到对应的Order变量
            LFPOrder = eval(LFPOrderVar);  % 获取LFP的Order值
            
            % 找到与当前order值匹配的trial索引
            matchIdx = find(ismember(LFPOrder', orderVal', 'rows'));  % 查找匹配的trial索引
            if ~isempty(matchIdx)
                % 如果匹配，将该LFP数据加入LFPMatrix中
                LFPMatrix = [LFPMatrix; LFPdata(:, matchIdx)'];  % 转置后横向拼接
            end
        end
        
        % 遍历选中的SPKC通道并收集对应neuron的SPKC子通道数据，同时记录Stim参数
        for j = 1:length(selectedSPKCChannels)
            baseSPKCChannel = ['SPKC', selectedSPKCChannels{j}];  % SPKC通道的基本名称
            
            % 找到所有对应的子通道，例如 SPKC01a, SPKC01b
            spkcSubChannels = SPKCvars(startsWith(SPKCvars, baseSPKCChannel));
            
            % 遍历每个SPKC子通道
            for k = 1:length(spkcSubChannels)
                SPKCvarName = spkcSubChannels{k};  % 获取SPKC子通道名称
                SPKCdata = eval(SPKCvarName);  % 获取SPKC struct数组
                SPKCOrderVar = ['Order', SPKCvarName(5:end)];  % 找到对应的Order变量
                SPKCOrder = eval(SPKCOrderVar);  % 获取SPKC的Order值
                
                % 找到与当前order值匹配的trial索引
                matchIdx = find(ismember(SPKCOrder', orderVal', 'rows'));  % 查找匹配的trial索引
                if ~isempty(matchIdx)
                    % 记录StimAngle和StimSpeed
                    stimAngleVar = ['StimAngle', SPKCvarName(5:end)];
                    stimSpeedVar = ['StimSpeed', SPKCvarName(5:end)];
                    stimAngleData = eval(stimAngleVar);  % 获取StimAngle数据
                    stimSpeedData = eval(stimSpeedVar);  % 获取StimSpeed数据
                    
                    trialStimAngle = [trialStimAngle; stimAngleData(matchIdx)];
                    trialStimSpeed = [trialStimSpeed; stimSpeedData(matchIdx)];
                    
                    % 如果有多个matchIdx，逐个处理
                    for m = 1:length(matchIdx)
                        % 如果 SPKCStruct 是空的，直接赋值
                        if isempty(fieldnames(SPKCStruct))
                            SPKCStruct = SPKCdata(matchIdx(m));
                        else
                            % 使用cat(1, ...) 垂直拼接结构体，避免字段数不匹配
                            SPKCStruct = cat(1, SPKCStruct, SPKCdata(matchIdx(m)));
                        end
                    end
                end
            end
        end
        % 只保留有数据的trial
        if ~isempty(LFPMatrix) || ~isempty(fieldnames(SPKCStruct))  % 检查SPKCStruct是否有数据
            % 增加有效trial计数
            validTrialCount = validTrialCount + 1;
            [dN,t]=binspikes(SPKCStruct,100,[0.8 2.2]);
            dN=transpose(dN(1:end-1,:));
            dN(dN~=0)=1;
            
            LFPvol=transpose(LFPMatrix(:,602:2201));
            [S,t,f]=mtspecgramc(LFPvol,[0.2 0.05],params);
            findx=(f<=30 & f<=50);
            S=S(:,findx,:);
            S=squeeze(mean(S,2));
            S=S';
            numNaNs=4;
            newS=[];
            for nn=1:size(S,1)
                originalArray=S(nn,:);
                filledArray = reshape([originalArray; NaN(numNaNs, length(originalArray))], 1, []);
                newS(nn,:)=filledArray;
            end
            % 保存当前trial的LFP和SPKC数据
            combinedLFPData{validTrialCount} = newS(:,2:end-numNaNs);  % 存储合并的LFP矩阵
            combinedSPKCData{validTrialCount} = dN;  % 存储合并的SPKC结构体
            
            % 保存当前trial的StimAngle和StimSpeed数据
            combinedStimAngle{validTrialCount} = unique(trialStimAngle);
            combinedStimSpeed{validTrialCount} = unique(trialStimSpeed);
        end
    end
end




% % 显示结果
% disp('保留的有效LFP数据:');
% disp(combinedLFPData);
% 
% disp('保留的有效SPKC数据:');
% disp(combinedSPKCData);
% 
% disp('StimAngle数据:');
% disp(combinedStimAngle);
% 
% disp('StimSpeed数据:');
% disp(combinedStimSpeed);

%%
Y_Obs=combinedLFPData{1};
N_Obs=combinedSPKCData{1};
clearvars -except N_Obs Y_Obs