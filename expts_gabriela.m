clear all;

mkdir('figures');
penalties = {'rad', 'radvol', 'con', 'radub'};
fnames = {'regen1.txt', 'regen2.txt', 'regen3.txt'};

vols = [];
for i=1:length(penalties)
    penalty = penalties{i};
    for j=1:length(fnames)
        fname = fnames{j};
        label = sprintf('%s_%s', fname(1:end-4), penalty);
        ddt_demo;
        saveas(h, sprintf('figures/%s.png', label));
        vols = [vols sprintf('\t%s=%f', label, stats.vol)];
        close all;
    end
end
disp(splitstring(vols));
