function [Parameters_fitted, fval] = CircuitModel_FitModel_ZZ (nUnit, ContextModulation, animal)

%set upper and lower bounds for each parameter: [alpha alpha0_NR sigma kappa WI WE WI_S WE_S Context(c)]
%upper bounds can be reduced from infinity for increased accuracy.
UpperBound = [Inf Inf Inf Inf 100 100 100 100 1 ];
LowerBound = [0 0 0 0 0 0 0 0 0 ];

%set the plausible upper and lower bounds for each parameter to increase
%accuracy
PlausibleLowerBound = [1e-3 1e-4 0.05 1.2 1e-3  1e-3  1e-3  1e-3 0];
PlausibleUpperBound = [5e-2 1e-2 0.4 4 2 2 2 2 0.8];

%select function being minimized
CostFun = @(x,animal,ContextModulation) CircuitModel_CostFun_ZZ(x,animal,ContextModulation,nUnit);

%select a random starting point within the plausible bounds (fitting should be run with ~10 different
%starting points).
x0 = PlausibleLowerBound + (PlausibleUpperBound-PlausibleLowerBound).*rand(1,numel(PlausibleLowerBound));
%run optimization
options = [];
options.UncertaintyHandling = 0;
[Parameters_fitted,fval] = bads(CostFun,x0,LowerBound,UpperBound, PlausibleLowerBound, PlausibleUpperBound,[],options,animal,ContextModulation);

end 