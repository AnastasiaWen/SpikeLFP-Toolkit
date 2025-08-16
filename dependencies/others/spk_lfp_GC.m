function [S1,S2,C,Phi,tt,f,GCs2p,GCp2s,ctl]=spk_lfp_GC(lfp,spike,params,movingwin)
%算GC
% lfp: T x TRLS
% spk: structure for spike timing for TRLS
ctl=0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%% spike-lfp spectrum analysis analysis
%Set the parameters of the MTM.
Fs = params.Fs;				%... sampling frequency
fpass = params.fpass;       %... fpass
pad = params.pad;             %... padding.

%time window and fft params
%movingwin=[0.15 0.01];      %Time window length and step
Nwin=round(Fs*movingwin(1)); % number of samples in window
Nstep=round(movingwin(2)*Fs); % number of samples to step through

N=size(lfp,1);Ntrl=size(lfp,2); % size of data   
winstart=1:Nstep:N-Nwin+1;      % time window 
nw=length(winstart);
nfft=max(2^(nextpow2(Nwin)+pad),Nwin);

% frequencies where spectrum is evaluated
f=0:Fs/nfft:Fs;
findx=find(f>=fpass(1) & f<=fpass(end));
f=f(findx);
Nf=length(f);

% initialize 
C=zeros(nw,Nf); %Coherence
Phi=zeros(nw,Nf);
S12=zeros(nw,Nf);S21=zeros(nw,Nf); S1=zeros(nw,Nf); S2=zeros(nw,Nf);

for n=1:nw
   indx=winstart(n):winstart(n)+Nwin-1;
   t=indx/Fs;
   datap=lfp(indx,:);
   datas=extractspk(spike,[t(1) t(end)]);
   [c,phi,s12,s2,s1,f]=coherencycpt(datap,datas,params,0,t); %Compute the MTM coherence and spectrum
   if isempty(find(s1~=0))
       ctl=1;
       tt=[];f=[];GCs2p=[];GCp2s=[];
       return ;
   end
   C(n,:)=c';Phi(n,:)=phi';
   S1(n,:)=s1';S2(n,:)=s2';S12(n,:)=s12';S21(n,:)=s12';
   S(1,2,:)=s12';S(2,1,:)=s12';
   S(1,1,:)=s1';S(2,2,:)=s2'; % spectral matrix for GC calculation
   [gc_s2p,gc_p2s]=GrangerCausality(S,f,Fs);
   GCp2s(n,:)=gc_p2s;GCs2p(n,:)=gc_s2p;
end

%%%%%%%%%%%%%%%%%%%%%%%%% Computing Granger causality using spectrum analysis 
tt=(winstart+round(Nwin/2))/Fs;                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Plot
% figure;
% plot(tt, sum(GCs2p,2))
% hold on
% plot(tt, sum(GCp2s,2),'r')
end