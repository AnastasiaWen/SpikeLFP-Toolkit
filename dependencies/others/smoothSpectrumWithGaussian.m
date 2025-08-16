function smoothedMatrix = smoothSpectrumWithGaussian(magnitudeMatrix, sigma)
    % magnitudeMatrix: 输入的频谱矩阵
    % sigma: 高斯平滑的标准差
    
    % 创建高斯核
    size = 2*round(3*sigma) + 1; % 根据 sigma 创建适当大小的核
    x = linspace(-size / 2, size / 2, size);
    gaussKernel = exp(-x .^ 2 / (2 * sigma ^ 2));
    gaussKernel = gaussKernel / sum(gaussKernel); % 归一化
    
    % 对频谱矩阵进行卷积平滑
    smoothedMatrix = conv2(magnitudeMatrix, gaussKernel, 'same');
end