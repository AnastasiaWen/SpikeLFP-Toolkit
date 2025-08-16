function p=anova2_single_neuro(fr)
% [numTextures, numSpeeds] = size(fr); % 获取纹理和速度的数量
% 
% % 将 firing rate 矩阵展开为向量
% fr_values = fr(:);  % 将5x12的矩阵变成一个60x1的向量
% 
% % 创建对应的纹理和速度因子
% textures = repmat((1:numTextures)', numSpeeds, 1);  % 纹理标签，重复每个纹理12次
% speeds = repmat({'s1','s2','s3','s4','s5','s6','s7','s8','s9','s10','s11','s12'}, numTextures, 1);  % 速度标签，重复5次
% speeds = speeds(:);  % 将速度标签也转成列向量
% 
% % 2. 执行双因素方差分析
% [p, tbl, stats] = anovan(fr_values, {textures, speeds}, 'model', 2, ...
%     'varnames', {'Texture', 'Speed'});
% 
% % p 值是0.05为阈值时判断显著性的标准，如果 p 值小于0.05，则认为该因素显著
% 
% % 3. 查看结果
% disp('双因素ANOVA结果：');
% disp(tbl);
% 
% % 如果你想绘制可视化结果
% figure;
% boxplot(fr_values, {textures, speeds}, 'Labels', {'Texture', 'Speed'});
% xlabel('Conditions');
% ylabel('Firing Rate');
% title('Firing Rate by Texture and Speed');

% 4. 如果你想进一步查看 stats 的信息，可以使用multcompare来进一步比较
% multcompare(stats, 'Dimension', 1); % 针对纹理的多重比较
% multcompare(stats, 'Dimension', 2); % 针对速度的多重比较


% 假设 fr 是 5x12 的 cell 数组，每个 cell 中有 1×n 的观测值
% 5 表示5种纹理，12表示12种速度

% 初始化cell数据
% 你的 fr 是 5x12 的 cell 数组，每个 cell 中包含 1xn 的观测值向量
% 这里假设 fr 已经加载或初始化

% 初始化保存展开数据的向量
data = [];  % 用于保存所有观测值
texture_factor = [];  % 用于保存纹理因子
speed_factor = [];  % 用于保存速度因子

% 遍历 fr 数组，将数据展开为一个向量，同时记录对应的因子水平
for texture = 1:size(fr, 1)  % 遍历纹理
    for speed = 1:size(fr, 2)  % 遍历速度
        current_data = fr{texture, speed};  % 当前 cell 的数据
        
        % 添加当前数据到 data 向量中
        data = [data; current_data(:)];
        
        % 根据当前 texture 和 speed 的索引，生成相应的因子水平
        texture_factor = [texture_factor; repmat(texture, length(current_data), 1)];
        speed_factor = [speed_factor; repmat(speed, length(current_data), 1)];
    end
end

% 使用 anovan 进行双因素方差分析，考虑交互作用
[p] = anovan(data, {texture_factor, speed_factor}, ...
                         'model', 'interaction', ...
                         'varnames', {'Texture', 'Speed'},'display','off');
% 
% % 显示分析结果
% disp('ANOVA table:');
% disp(tbl);

