function lfpnorm=lfp_norm(lfp)
%note: LFP format: time amplitude x trial
for i=1:size(lfp,2)
    tmp=lfp(:,i);
    meantmp=mean(tmp);
    lfpnorm(:,i)=tmp-meantmp;
end