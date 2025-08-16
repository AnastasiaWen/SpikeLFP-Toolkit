function [faxis,Sxx]=fftspecc(data,p)
Fs = 1000; % 采样频率

pl=0;
if nargin>1
    if p=='p'
        pl=1;
    end
end


% 计算频谱
L = length(data); % 数据长度

t = 1:L; % 时间向量
dt=1/Fs;
t=t/Fs;
T=t(end);
%data = data+sin(2*pi*30*t) + 0.5*sin(2*pi*120*t); % 包含两个正弦波的示例数据

xf = fft(data); % 进行快速傅里叶变换
Sxx=2*dt^2/T*(xf.*conj(xf));
Sxx = Sxx(1:L/2+1); % 获取正频率范围
df=1/max(T);
fNQ=1/dt/2;
faxis=(0:df:fNQ);

if pl==1
%绘制频谱
figure;
subplot(1,2,1);
plot(faxis,10*log10(Sxx)) 
xlim([0 200]);ylim([-150 0])
title('单侧振幅谱')
xlabel('频率 (Hz)')
ylabel('振幅')

subplot(1,2,2);
plot(t,data);
%xlim([0 1]);
title('原始信号')
xlabel('time')
ylabel('振幅')
end

end