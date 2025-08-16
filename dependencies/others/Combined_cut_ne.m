function Combined_cut_ne(sitepath) 
%%
combinepath=fullfile(sitepath,'combined');
cd(combinepath);
load('Combined_cut.mat');
matData=load('Combined_cut.mat');

savePath = fullfile(sitepath,'combined\Combined_cut_new.mat');
if isfile(savePath)
    delete(savePath);
end
%mkdir(savePath);
%%
uniqSpeed=[5,10,20,25,40,50,80,100,160,200,240,320];
Fs=1000;
uniqAngle=[0,45,90,135,180];% 0,45,90,135,180 per 24times

lowspeed=uniqSpeed(1:6);
highspeed=uniqSpeed(7:12);

deltrial=[];
vars_lfp = who('-regexp', '^LFP');
vars_spk = who('-regexp', '^SPKC');
vars_ang = who('-regexp', '^StimAngle');
vars_spd = who('-regexp', '^StimSpeed');
vars_wf = who('-regexp', '^Waveform');
vars_wft = who('-regexp', '^Wavetime');
vars_or = who('-regexp', '^Order');

for q=1:numel(vars_spk)
varName=vars_spk{q};
totalangle=matData.(['StimAngle',varName(end-2:end)]);
totalspeed=matData.(['StimSpeed',varName(end-2:end)]);
totalnum=zeros(5,12);
for i=1:5
    for j=1:12
        rc=find(totalangle==uniqAngle(i)&totalspeed==uniqSpeed(j));
        totalnum(i,j)=numel(rc);
    end
end
eval(['tiralnum',varName(end-2:end),'=totalnum;'])

sn=sum(totalnum,1);
an=sum(totalnum,2);
if numel(find(sn<5))~=0 | numel(find(an<5))~=0 
    % ldeltrial=[];
    % orderspk=eval(['Order',varName(end-2:end)]);
    % orderlfp=eval(['Order',varName(end-2:end-1)]);
    % tmpspk=eval(varName);
    % for jj=1:numel(tmpspk)
    %     lfpindx=find(all(orderlfp == orderspk(:,jj), 1));
    %     ldeltrial=[ldeltrial,lfpindx];
    % end
    % eval(['LFP',varName(end-2:end-1),'(:,ldeltrial)=[];']);
    % eval(['Order',varName(end-2:end-1),'(:,ldeltrial)=[];'])
    deltrial=[deltrial,q];
    % vars_spk(strcmp(vars_spk, varName))=[];
    vars_ang(strcmp(vars_ang, ['StimAngle',varName(end-2:end)]))=[];
    vars_spd(strcmp(vars_spd, ['StimSpeed',varName(end-2:end)]))=[];
    vars_or(strcmp(vars_or, ['Order',varName(end-2:end)]))=[];
    vars_wft(strcmp(vars_wft, ['Wavetime',varName(end-2:end)]))=[];
    vars_wf(strcmp(vars_wf, ['Waveform',varName(end-2:end)]))=[];
end

end

vars_spk(deltrial)=[];
% % vars_lfp(deltrial)=[];
% vars_ang(deltrial)=[];
% vars_spd(deltrial)=[];
% % vars_or(deltrial)=[];
% vars_wft(deltrial)=[];
% vars_wf(deltrial)=[];

savevar=[vars_lfp;vars_spk;vars_ang;vars_spd;vars_wf;vars_wft;vars_or];
save(savePath, savevar{:});
end

