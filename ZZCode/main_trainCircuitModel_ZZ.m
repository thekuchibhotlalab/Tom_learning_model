allUnits = [3, 5, 8, 10, 15, 20:10:150];
animal = 'average';
ContextModulation = {'excitatory', 'inhibitory', 'threshold', 'gain'} ;
repetition = 20;

for context = ContextModulation
    for nUnit = allUnits
        allParam = zeros(1,repetition);
        allCost = zeros(1,repetition);
        for i = 1:repetition
            [Parameters_fitted, fval] = CircuitModel_FitModel_ZZ (nUnit, context, animal);
            allParam(i) = Parameters_fitted;
            allCost(i) = fval;
        end
        
    end
    save(['Unit' num2str(nUnit,'%03d') '_' context '_' animal '.mat'],'allParam','allCost')
end

