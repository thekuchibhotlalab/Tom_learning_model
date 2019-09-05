nUnit = 100;
repetition = 20;
load('data/unit100_inhibitory_average.mat');

[~,index] = min(allCost);

bestParam = allParam(:,index);

[target_corr, FA_rate, target_corr_probe, FA_rate_probe, WESTORE, WISTORE, W, act] = ...
    CircuitModel_Stochastic_ZZ(bestParam,nUnit,'average','inhibitory','on', repetition);

%% try training
t0 = tic;
nUnit = 100;
animal = 'average';

[Parameters_fitted, fval] = CircuitModel_FitModel_ZZ (nUnit, 'inhibitory', animal);
t = toc;
disp(t-t0)