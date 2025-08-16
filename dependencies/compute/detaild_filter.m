function [delta,theta,alpha,beta,gamma]=detaild_filter(lfp)
Fs=1000;

 % delta (1–3 Hz), 
 d = designfilt('bandpassiir', ...       % Response type
       'FilterOrder',4, ...    % Frequency constraints
       'HalfPowerFrequency1',1, ...
       'HalfPowerFrequency2',3, ...
       'DesignMethod','butter', ...      % Design method  
       'SampleRate',Fs) ;              % Sample rate
delta=filtfilt(d,lfp);

  % theta (4–8 Hz),
 d = designfilt('bandpassiir', ...       % Response type
       'FilterOrder',4, ...    % Frequency constraints
       'HalfPowerFrequency1',4, ...
       'HalfPowerFrequency2',8, ...
       'DesignMethod','butter', ...      % Design method  
       'SampleRate',Fs) ;              % Sample rate
theta=filtfilt(d,lfp);

  % alpha (9–12 Hz)
 d = designfilt('bandpassiir', ...       % Response type
       'FilterOrder',4, ...    % Frequency constraints
       'HalfPowerFrequency1',9, ...
       'HalfPowerFrequency2',12, ...
       'DesignMethod','butter', ...      % Design method  
       'SampleRate',Fs)  ;             % Sample rate
alpha=filtfilt(d,lfp);

  % beta (12–30 Hz),
 d = designfilt('bandpassiir', ...       % Response type
       'FilterOrder',4, ...    % Frequency constraints
       'HalfPowerFrequency1',12, ...
       'HalfPowerFrequency2',30, ...
       'DesignMethod','butter', ...      % Design method  
       'SampleRate',Fs)    ;           % Sample rate
beta=filtfilt(d,lfp);

  %  gamma (30–80 Hz)
 d = designfilt('bandpassiir', ...       % Response type
       'FilterOrder',4, ...    % Frequency constraints
       'HalfPowerFrequency1',30, ...
       'HalfPowerFrequency2',80, ...
       'DesignMethod','butter', ...      % Design method  
       'SampleRate',Fs)  ;             % Sample rate
gamma=filtfilt(d,lfp);
