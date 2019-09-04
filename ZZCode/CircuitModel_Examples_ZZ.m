%% Examples of context-dependent circuit model implementation
%
% EXAMPLE 1: Figure  2b; Stochastic run, gain modulation, average behavioral data.
% EXAMPLE 2: Figure  2d; Stochastic run, inhibitory scaling, average behavioral data.
% EXAMPLE 3: Figure  3ab; Stochastic run, inhibitory scaling,  mouse202.
% EXAMPLE 4: Fitting procedure; Fitting an inhibitory scaling model to
% mouse 203.
%
% Please ensure that the relevant datafiles ('individual_behavior_v2.mat' |
%   'average_behavior_v2.mat' | 'optimized_parameters.mat') are in the matlab
% path. These files may be retrieved at http://froemkelab.med.nyu.edu/

% EXAMPLE 4 requires installation of BADS (Bayesian Adaptive Direct Search,
% Acerbi & Ma, 2017). BADS may be retrieved at https://github.com/lacerbi/bads

%   possible animals:
%
%      'average'
%      'kkjm202'
%      'kkjm203'
%      'kkjm204'
%      'kkjscam015'
%      'kkpv10'
%      'kkpv11'
%      'kkpv13'
%{
%% EXAMPLE 1: Figure  2b; Stochastic run, gain modulation, average behavioral data.

load optimized_parameters.mat;

animal = 'average';

ContextModulation = 'gain'; %alt: inhibitory, excitatory,

%set to 'off' if only output is needed.
Visualization = 'on';

%retrieve parameters: [alpha alpha_NR sigma kappa WI WE WI_S WE_S Context(c)]
Parameters = optimized_parameters.(ContextModulation).(animal);

%run stochastic model using pre-defined parameters
[ReinforcedHit,ReinforcedFA,ProbeHit,ProbeFA] = CircuitModel_Stochastic(Parameters,animal,ContextModulation,Visualization);
%% EXAMPLE 2: Figure  2d; Stochastic run, inhibitory scaling, average behavioral data.

load optimized_parameters.mat;

animal = 'average';

ContextModulation = 'inhibitory';

%set to 'off' if only output is needed.
Visualization = 'on';

%retrieve parameters: [alpha alpha_NR sigma kappa WI WE WI_S WE_S Context(c)]
Parameters = optimized_parameters.(ContextModulation).(animal);

%run stochastic model using pre-defined parameters
CircuitModel_Stochastic(Parameters,animal,ContextModulation,Visualization);

%% EXAMPLE 3: Figure  3ab; Stochastic run, inhibitory scaling,  mouse202.

load optimized_parameters.mat;

animal = 'kkjm202';

ContextModulation = 'inhibitory';

%set to 'off' if only output is needed.
Visualization = 'on';

%retrieve parameters: [alpha alpha_NR sigma kappa WI WE WI_S WE_S Context(c)]
Parameters = optimized_parameters.(ContextModulation).(animal);

%run stochastic model using pre-defined parameters
CircuitModel_Stochastic(Parameters,animal,ContextModulation,Visualization);
%}
%% EXAMPLE 4: Fitting procedure; Fitting an inhibitory scaling model to mouse 203.

animal = 'average';
ContextModulation = 'inhibitory';

%set upper and lower bounds for each parameter: [alpha alpha0_NR sigma kappa WI WE WI_S WE_S Context(c)]
%upper bounds can be reduced from infinity for increased accuracy.
UpperBound = [Inf Inf Inf Inf 100 100 100 100 1 100 100 ];
LowerBound = [0 0 0 0 -100 -100 -100 -100 0 -100 -100 ];

%set the plausible upper and lower bounds for each parameter to increase
%accuracy
PlausibleLowerBound = [1e-3 1e-4 0.05 1.2 -2 -2 -2 -2 0 -2 -2];
PlausibleUpperBound = [5e-2 1e-2 0.4 4 2 2 2 2 0.8 2 2];

%select function being minimized
CostFun = @(x,animal,ContextModulation) CircuitModel_CostFun_ZZ(x,animal,ContextModulation);

%select a random starting point within the plausible bounds (fitting should be run with ~10 different
%starting points).
x0 = PlausibleLowerBound + (PlausibleUpperBound-PlausibleLowerBound).*rand(1,numel(PlausibleLowerBound));
%run optimization
options = [];
options.UncertaintyHandling = 0;
[Parameters_fitted,fval] = bads(CostFun,x0,LowerBound,UpperBound, PlausibleLowerBound, PlausibleUpperBound,[],options,animal,ContextModulation);
