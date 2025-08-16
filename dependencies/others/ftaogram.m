function [F,phi,tc] = ftaogram(data_spk,data_lfp,t,Tc,Tinc,Tw,wb,fs)

Tc(1) = ceil(Tc(1)/Tinc)*Tinc;
Tc(2) = floor(Tc(2)/Tinc)*Tinc;
tc = Tc(1):Tinc:Tc(2);
Fch=0;
for tt=1:length(tc)
  T = [tc(tt)-Tw/2 tc(tt)+Tw/2];
  if T(1)<Tc(1)||T(2)>Tc(2)
      continue;
  end
  if tt == 1
    [SS,phi] = fta(data_spk, data_lfp, fs,wb,t,T);
    F = zeros(length(tc),length(SS));
    Fch=1;
  else
    [SS,phi] = fta(data_spk, data_lfp, fs,wb,t,T);
    if Fch==0
        F = zeros(length(tc),length(SS));
        Fch=1;
    end
  end
  F(tt,:) = mean(SS,1);
end