vars=who;
vars_LFP=who('-regexp','^before');
lfp_data=[];
for i=1:15
    lfp_data=[lfp_data;eval(vars_LFP{i})];
end