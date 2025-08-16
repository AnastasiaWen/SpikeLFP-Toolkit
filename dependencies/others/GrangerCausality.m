function  [GCs2p,GCp2s]=GrangerCausality(S,f,Fs)

%nw=size(S,3);% num of time window
Nf=length(f);
    %%%%%%%%%%%%%%%%%%%%spectral matrix factorization (Willson,1972 && Dhamala,2008)
    f_indx=0;
    for jj=f
        f_indx=f_indx+1;
        Sarr(:,:,f_indx)=S(:,:,f_indx);
        if(f_indx>1)
          Sarr(:,:,2*(Nf-1)+2-f_indx)=S(:,:,f_indx).';%'
        end
    end
    
    % perform ifft to obtain gammas
    for k1=1:2
        for k2=1:2
            gam(k1,k2,:)=real(ifft(squeeze(Sarr(k1,k2,:)))*Fs);
        end
    end
    gam0 = gam(:,:,1); h = chol(gam0); 

    for ind =1:size(Sarr,3),
           psi(:,:,ind) = h; % initialization of  for the 1st iteration
    end
    
    for iter = 1:100
         for ind = 1:size(Sarr,3),
           g(:,:,ind)=inv(psi(:,:,ind))*Sarr(:,:,ind)*inv(psi(:,:,ind))'+eye(2);%'
         end
         gplus = PlusOperator(g,2,Fs,f); % Eq 3.1
       psiold=psi;
       for k = 1:size(Sarr,3),
              psi(:,:,k) = psi(:,:,k)*gplus(:,:,k);
              psierr(k)=norm(psi(:,:,k)-psiold(:,:,k),1);
       end
       iter;
       psierrf=mean(psierr);
       if(psierrf<1E-12),break;end;     
    end 

    for k = 1:Nf
          Snew(:,:,k) = squeeze(psi(:,:,k)*psi(:,:,k)'); %% Snew:improved spectral density matrix 
    end

    for k1=1:2
        for k2=1:2
            gamtmp(k1,k2,:)=real(ifft(squeeze(psi(k1,k2,:))));
        end
    end
    A0=squeeze(gamtmp(:,:,1)); % this is psi_1
    A0inv=inv(A0);
    Znew=A0*A0.'*Fs; % noise covariance(channels x channels x frequency),

    for k = 1:Nf
          Hnew(:,:,k) = squeeze(psi(:,:,k))*A0inv; % transfer function
    end

    %%%%% estimate Granger causality from spectral density matrix, transfer function and noise covariance matrix.
    clear f_indx jj;
    f_indx = 0;
    for mm = f
        f_indx = f_indx + 1;
        Zp2s = Znew(2,2) - Znew(1,2)^2/Znew(1,1); %corrected noise covariance
        Zs2p = Znew(1,1) - Znew(2,1)^2/Znew(2,2);
        GCp2s(f_indx) = log(abs(Snew(1,1,f_indx))/abs(Snew(1,1,f_indx)-(Hnew(1,2,f_indx)*Zp2s*conj(Hnew(1,2,f_indx)))/Fs)); %Geweke's original measure
        GCs2p(f_indx) = log(abs(Snew(2,2,f_indx))/abs(Snew(2,2,f_indx)-(Hnew(2,1,f_indx)*Zs2p*conj(Hnew(2,1,f_indx)))/Fs));
        %GCp2s(kk,f_indx) = log(abs(Snew(1,1,f_indx))/abs(Snew(1,1,f_indx)-(abs(Hnew(1,2,f_indx))*Zp2s)/Fs)); %Geweke's original measure
        %GCs2p(kk,f_indx) = log(abs(Snew(2,2,f_indx))/abs(Snew(2,2,f_indx)-(abs(Hnew(2,1,f_indx))*Zs2p)/Fs));
    end
    


%tt=(winstart+round(Nwin/2))/Fs;