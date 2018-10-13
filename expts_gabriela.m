%clear all;

mkdir('figures');
penalties = {'rad', 'radvol', 'con', 'radub'};
fnames = {'process1.txt', 'process2.txt', 'process3.txt'};
thresholds = [2.06, 26.27, 82.86];

vols = [];
for i=1:length(penalties)
    penalty = penalties{i};
    for j=1:length(fnames)
        fname = fnames{j};
        threshold = thresholds(j);
        label = sprintf('%s_%s', fname(1:end-4), penalty);
        ddt_demo;
        saveas(h, sprintf('figures/%s.png', label));
        vols = [vols sprintf('\t%s=%f', label, stats.vol)];
        close all;
    end
end
disp(splitstring(vols));
