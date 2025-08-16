function CutTrial_data(sitepath) 
%%
combinepath=fullfile(sitepath,'combined');
cd(combinepath);
load Combined.mat

savePath = fullfile(sitepath,'combined\Combined_cut.mat');
if isfile(savePath)
    delete(savePath);
end
%mkdir(savePath);
%%

vars = who;
vars_lfp = who('-regexp', '^LFP');
vars_spk = who('-regexp', '^SPKC');
vars_ang = who('-regexp', '^StimAngle');
vars_spd = who('-regexp', '^StimSpeed');
vars_wf = who('-regexp', '^Waveform');
vars_wft = who('-regexp', '^Wavetime');
vars_or = who('-regexp', '^Order');

%生成timestamp结构体 ( 1 x trial)
for i = 1:numel(vars_spk)
    tmpname=vars_spk{i};tmpname=tmpname(end-2:end);
    tmpspk=eval(vars_spk{i});
    tmporder=eval(['Order',tmpname]);
    tmpspeed=eval(['StimSpeed',tmpname]);
    tmpangle=eval(['StimAngle',tmpname]);
    [merged_Angle,merged_Speed,merged_SPKC, merged_Order]  = merge_trials(tmpspk, tmporder, tmpangle, tmpspeed);
    eval([vars_spk{i},'=merged_SPKC;']);
    eval(['Order',tmpname,'=merged_Order;']);
    eval(['StimAngle',tmpname,'=merged_Angle;']);
    eval(['StimSpeed',tmpname,'=merged_Speed;']);
end

for i = 1:numel(vars_lfp)
    deltrial=[];
    ldeltrial=[];
    eval([ 'tmplfp=',vars_lfp{i},';']);
    lfporder_name=vars_lfp{i};
    lfporder_name=['Order',lfporder_name(end-1:end)];
    lfporder=eval(lfporder_name);

    tmpname=vars_lfp{i};tmpname=tmpname(end-1:end);tmpname=['^SPKC',tmpname];
    vars_sp = who('-regexp', tmpname);
    sp_indx=cell(1,numel(vars_sp));
    

    for j=1:size(tmplfp,2)
        values = tmplfp(:,j);
        if isnan(values(1))
            deltrial=[deltrial,lfporder(:,j)];
            orderval=lfporder(:,j);
            ldeltrial=[ldeltrial,j];
            for kk=1:numel(vars_sp)
                tmpname=vars_sp{kk};tmpname=tmpname(end-2:end);
                orderspk=eval(['Order',tmpname]);
                spi=find(all(orderspk == orderval, 1));
                if ~isempty(spi)
                    sp_indx(kk)={[sp_indx{kk},spi]};
                end
            end
        
        else
            orderval=lfporder(:,j);
            novalidnum=0;
            for ss=1:numel(vars_sp)
                tmpname=vars_sp{ss};tmpname=tmpname(end-2:end);
                orderspk=eval(['Order',tmpname]);
                eval([ 'tmpspk=',vars_sp{ss},';']);
                spi=find(all(orderspk == orderval, 1));
                if ~isempty(spi)
                    values=tmpspk(spi).times;
                    if isempty(values)
                        sp_indx(ss)={[sp_indx{ss},spi]};
                        novalidnum=+novalidnum+1;
                    end
                else
                    novalidnum=+novalidnum+1;
                end
            end
        
            if novalidnum==numel(vars_sp)
                for kk=1:numel(vars_sp)
                    tmpname=vars_sp{kk};tmpname=tmpname(end-2:end);
                    orderspk=eval(['Order',tmpname]);
                    spi=find(all(orderspk == orderval, 1));
                    if ~isempty(spi)
                        sp_indx(kk)={[sp_indx{kk},spi]};
                    end
                end
                ldeltrial=[ldeltrial,j];
            end
            
            
            
    
        end
    end

    eval([vars_lfp{i},'(:,ldeltrial)=[];']);
    eval([lfporder_name,'(:,ldeltrial)=[];']);
    for ss=1:numel(vars_sp)
        tmpname=vars_sp{ss};tmpname=tmpname(end-2:end);
        orderspk=eval(['Order',tmpname]);
        indx=sp_indx{ss};
        eval([vars_sp{ss},'(indx)=[];']);
        eval(['StimAngle',tmpname,'(indx)=[];']);
        eval(['StimSpeed',tmpname,'(indx)=[];']);
        eval(['Order',tmpname,'(:,indx)=[];']);
    end

end


savevar=[vars_lfp;vars_spk;vars_ang;vars_spd;vars_wf;vars_wft;vars_or];
save(savePath, savevar{:});
end

