function preferdir(meanFR)
% 定义方向角度（弧度制）
angles = deg2rad([0, 45, 90, 135, 180]);
% angles = deg2rad([0, 45, 90, 135]);
% angles = deg2rad([ 45, 90, 135, 180]);

% 定义 von Mises 拟合函数
vonMises = @(b, theta) b(1) * exp(b(2) * cos(2*(theta - b(3))));  % b(1)=A, b(2)=kappa, b(3)=mu

% 初始参数估计：幅值 A、集中参数 kappa 和优选方向 mu
initialParams = [max(meanFR), 1, pi/2];  % A = max firing rate, kappa = 1, mu = 90度

% 使用最小二乘法拟合 von Mises 函数
optimOptions = optimset('MaxFunEvals', 1000, 'MaxIter', 1000);  % 设定最大迭代次数
fitParams = lsqcurvefit(vonMises, initialParams, angles, meanFR, [], [], optimOptions);

% 提取拟合参数
A_fit = fitParams(1);      % 幅值
kappa_fit = fitParams(2);  % 调制强度
mu_fit = fitParams(3);     % 优选方向 (弧度)

% 将 mu_fit 转换为度数
preferred_orientation = rad2deg(mu_fit);

% 绘制拟合结果与原始数据的比较
theta_fit = linspace(0, pi, 100);  % 生成拟合的角度范围
fittedFR = vonMises(fitParams, theta_fit);  % 计算拟合曲线

figure;
plot(rad2deg(angles), meanFR, 'bo', 'MarkerFaceColor', 'b');  % 原始数据
hold on;
plot(rad2deg(theta_fit), fittedFR, 'r-', 'LineWidth', 2);  % 拟合曲线
xlabel('Orientation (degrees)');
ylabel('Firing rate (Hz)');
title(['Preferred Orientation: ', num2str(preferred_orientation), '°']);
legend('Data', 'von Mises Fit');
grid on;

% 显示结果
fprintf('Preferred orientation (degrees): %.2f\n', preferred_orientation);
fprintf('Kappa (concentration parameter): %.2f\n', kappa_fit);