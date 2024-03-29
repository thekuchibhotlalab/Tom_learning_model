nUnit = 100;
repetition = 20;
load('data/unit100_inhibitory_average.mat');

[~,index] = min(allCost);

bestParam = allParam(:,index);
%%
% [target_corr, FA_rate, target_corr_probe, FA_rate_probe, WESTORE, WISTORE, W, act] = ...
%     CircuitModel_Stochastic_ZZ(bestParam,nUnit,'average','inhibitory','on', repetition);

[target_corr, FA_rate, target_corr_probe, FA_rate_probe, WESTORE, WISTORE, W, act, EIstore, ctxAct] = ...
    CircuitModel_Stochastic_ZZ(bestParam,nUnit,'average','inhibitoryCtx','on', repetition);

%%
% [target_corr, FA_rate, target_corr_probe, FA_rate_probe, WESTORE, WISTORE, W, act] = ...
%     CircuitModel_Stochastic_ZZ(bestParam,nUnit,'average','inhibitory','on', repetition);

[target_corr, FA_rate, target_corr_probe, FA_rate_probe, WESTORE, WISTORE, W, act, EIstore, ctxAct] = ...
    CircuitModel_Stochastic_ZZ_reverse(bestParam,nUnit,'average','inhibitoryCtx','on', repetition);
%%
cmap = redbluecmap;
newCmap = imresize(cmap, [128,3]);
newCmap = min(max(newCmap,0),1);

modStore = ctxAct{1};
targInhibScaleStore = ctxAct{2};
foilInhibScaleStore = ctxAct{3};
reinfTargStore = ctxAct{4};
probeTargStore = ctxAct{5};
reinfFoilStore = ctxAct{6};
probeFoilStore = ctxAct{7};

figure; 
plot(mean(targInhibScaleStore,2), 'LineWidth',2); hold on; plot(mean(foilInhibScaleStore,2), 'LineWidth',2)
xlim([1 60])
legend('target', 'foil')
title ('scale of modulation on ctx')
xlabel('blocks');

figure; imagesc(reinfFoilStore(:,:,1) - probeFoilStore(:,:,1))

figure; 
probeDiff = probeTargStore(:,:,1) - probeFoilStore(:,:,1);
reinfDiff = reinfTargStore(:,:,1) - reinfFoilStore(:,:,1);
act = mean(reinfDiff(:,1:3),2);
[~, sortIndex] = sort(act);

colorlim = [-0.8 0.8];
figure; subplot(131); imagesc(probeDiff(sortIndex,:))
colormap(newCmap); caxis(colorlim); title('probe (T-F)')
xlabel('block'); ylabel('neurons')
subplot(132); imagesc(reinfDiff(sortIndex,:))
colormap(newCmap); caxis(colorlim); title('reinforce (T-F)')
xlabel('block'); ylabel('neurons')
subplot(133); imagesc(probeDiff(sortIndex,:) - reinfDiff(sortIndex,:))
colormap(newCmap); caxis(colorlim); title('reinf - probe')
xlabel('block'); ylabel('neurons')
colorbar

figure;
%%

targEStore = EIstore{1};
foilEStore = EIstore{2};
targIStore = EIstore{3};
foilIStore = EIstore{4};


figure;
hold on;
plot(targEStore, 'Color', [0.9 0.9 0.9])
plot(foilEStore, 'Color', [0.9 0.9 0.9])
plot(targIStore, 'Color', [0.9 0.9 0.9])
plot(foilIStore, 'Color', [0.9 0.9 0.9])

h(1) = plot(mean(targEStore,2),'LineWidth',2); 
h(2) = plot(mean(foilEStore,2),'LineWidth',2);
h(3) = plot(mean(targIStore,2),'LineWidth',2);
h(4) = plot(mean(foilIStore,2),'LineWidth',2);
legend(h,{'E_t_a_r_g','E_f_o_i_l','I_t_a_r_g','I_f_o_i_l'},'Location','Best')
%% try training
t0 = tic;
nUnit = 100;
animal = 'average';

[Parameters_fitted, fval] = CircuitModel_FitModel_ZZ_reverse (nUnit, 'inhibitory', animal);
t = toc;
disp(t-t0)