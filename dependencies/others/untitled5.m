%%
close all
vars_lfp = who('-regexp', '^LFP');

for i=1:numel(vars_lfp)
    lfp=eval(vars_lfp{i});
    figure;
    plot(lfp);
    title(vars_lfp{i});
end
