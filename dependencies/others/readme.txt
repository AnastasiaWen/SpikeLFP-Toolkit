1.中间调用的函数有点杂，最好是把整个script添加，chronux有些代码我改过，但是没拷过来应该，找macbook里的mychronux文件夹添加。

2. 主要用到的是script/compute里的代码，batch_compute是批量计算代码，里面有调用的例子

3. 主要函数讲解
spk_feature   -- 展示每个neuron自己的特征图
channel_lfp   --  不同深度、不同频带的LFP功率
lfp_spectro    --   LFP的spectrogram
di_sp_neuro_tune   --  direction tuning   （返回计数值
oi_sp_neuro_tune   --  orientation tuning （5个方向
oi4_sp_neuro_tune   --  orientation tuning （4个方向
GC_Other_channel  --   granger causality计算
channel_cohegram --   coherence计算 （这和GC都是跨通道的multiscale






