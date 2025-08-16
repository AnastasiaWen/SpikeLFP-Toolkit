function n=spk_ts2bin(spk_ts,t)
%弃用，直接用chronux的更好
%n- trial X time
K=size(spk_ts,2);
n=zeros(K,length(t));
for i=1:K
    tmp=spk_ts(i).times;
    for j=1:length(tmp)
        tmpindx=find(abs(t-tmp(j))<eps);
        n(i,tmpindx)=1;
    end
end
end