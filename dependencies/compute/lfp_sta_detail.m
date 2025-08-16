function lfp_sta_detail(sitepath)
% %coorperate with lfp_wavelet_spectro.m
% close all;
% combinepath=fullfile(sitepath,'combined');
% cd(combinepath);
% load p_data.mat
% 
%  savePath = fullfile(sitepath,'results\spectro\wavelet');
% rmdir(savePath,'s'); % 选择
% savePath = fullfile(sitepath,'results\spectro\wavelet\nocut');
% mkdir(savePath);

% coorperate with lfp_wavelet_spectro.m
close all;
combinepath=fullfile(sitepath,'combined');
cd(combinepath);
if isfile('Combined_cut_new.mat')
    load Combined_cut_new.mat
else
    load Combined_cut.mat
end


savePath = fullfile(sitepath,'results\sta');
if isfolder(savePath)
    rmdir(savePath,'s');
end
mkdir(savePath);


%Set the parameters of the MTM.
Fs = 1000;                  %Sampling frequency
TW = 3;						%Choose time-bandwidth product of 3.
ntapers = 2*TW-1;			%Choose the # of tapers.
fpass = [0 200];
pad = 2;

params.Fs = Fs;				%... sampling frequency
params.tapers=[TW,ntapers];	%... time-bandwidth product & tapers.
params.fpass=fpass;       %... fpass
params.pad = pad;             %... padding.
params.trialave = 1;		%... trial average.

%speparam=params;
%speparam.fpass=[0 120];

%时间窗参数length steplength
%movingwin=[0.3 0.03];
%%
Fs=1000;

vars = who;
vars_lfp = who('-regexp', '^LFP');
vars_spk = who('-regexp', '^SPKC');

totallfp_sta=struct();
for i=1:numel(vars_lfp)
    lfp=eval(vars_lfp{i});
    lfp=lfp_norm(lfp); %norm lfp
    spk=eval(vars_spk{i});

    smp=0:1/1000:(2.5-1/1000);
    T=[1.2 1.8];
    D=[-0.04 0.04];
    plt='n';
    [s,t] = sta(spk,lfp,smp,plt,0,T,D);
    s=mean(s);
    % figure;
    % plot(t,s-mean(s) ,'k', 'LineWidth', 2);
    % xlabel('Time (s)');ylabel('Voltage (V)');
    % saveas(gcf, fullfile(savePath, ['total_',vars_lfp{i},'.png']));
     totallfp_sta.(vars_lfp{i})=s;
    %[delta,theta,alpha,beta,gamma]=detaild_filter(s);
end
save(fullfile(savePath,'totallfp_sta.mat'),'-struct','totallfp_sta');
% save(fullfile(savePath,'lowlfp_spec.mat'),'-struct','lowlfp_spec');
