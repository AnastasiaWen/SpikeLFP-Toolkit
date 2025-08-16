% Fisher合并p值函数
function overall_p_value = fisher_combined_p(p_values)
    % Fisher合并p值方法
    % 输入: p_values 是一个时间维度的p值数组
    % 输出: 合并后的总体p值
    
    % 计算Fisher检验统计量
    chi_square_stat = -2 * sum(log(p_values(:)));  % 将所有p值取对数并乘以-2
    degrees_of_freedom = 2 * numel(p_values);  % 自由度为2倍的p值数量
    
    % 通过卡方分布计算总体p值
    overall_p_value = 1 - chi2cdf(chi_square_stat, degrees_of_freedom);
end