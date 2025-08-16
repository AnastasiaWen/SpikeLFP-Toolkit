function combinedetail_order(trialpath)
cd(fullfile(trialpath,'combined'));
if isfile('combinedDetail2.mat')
matdata=load('CombineDetail2.mat');
else
matdata=load('CombineDetail.mat');
end
vars=fieldnames(matdata);
orderData=struct();

for i=1:numel(vars)
tmpname=vars{i};
C=matdata.(vars{i});
% 使用正则表达式提取 'TRIAL' 后面的数字
trialNumbers = cellfun(@(x) str2double(regexp(x, '(?<=TRIAL)\d+', 'match')), C);

% 根据提取的数字对 cell 数组进行排序
[~, sortedIdx] = sort(trialNumbers);

% 获取排序后的结果
sortedC = C(sortedIdx);
matdata.(vars{i})=sortedC;

% 显示排序后的结果
% disp('排序后的数组：');
% disp(sortedC);
% 
% firstRow = repelem(trialNumbers, 120);
% secondRow = repmat(1:120, 1, numel(C));
% resultArray = [firstRow; secondRow];
end

save('CombineDetail.mat','-struct',"matdata");
end