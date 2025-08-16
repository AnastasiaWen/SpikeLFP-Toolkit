function smoothedMatrix = smoothSpectrum(magnitudeMatrix, windowSize)
    % magnitudeMatrix: 输入的频谱矩阵
    % windowSize: 平滑窗口大小（移动平均窗口）
    
    % 对每一列（即每个频率）进行移动平均平滑
    smoothedMatrix = transpose(movmean(magnitudeMatrix', windowSize, 2));
end