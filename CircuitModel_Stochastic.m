function varargout = CircuitModel_Stochastic(params,animal,ContextModulation,visualization);
%--------------------------------------------------------------------------
%           Stochastic 3-Neuron Circuit Model (v1.0.2).
%--------------------------------------------------------------------------
%
%   This function runs a stochastic version of the circuit model
%   introduced by Kuchibhotla & Hindmarsh Sten (2017) to account for
%   context-dependent learning trajectories. Stochasticity arises from
%   noise in the decision making function, as well as random tone
%   selection. Synaptic updates are made on a trial-by-trial basis.
%
%   Behavioral datafiles ('individual_behavior.mat' ||
%   'average_behavior.mat') must be added to the Matlab path before function
%   can be run.
%
%   'Visualization' should be set to 'on' if model output is to be
%   visualized in matlab. If set to 'off' the result of repeated model runs
%   is outputed.
%
%   Input vector 'params' contains the 9 parameters for the model, and
%   should be ordered as follows:
%
%    params = [alpha alpha_NR sigma kappa WI WE WI_S WE_S Context(c)]
%
%   
%   possible models are as follows, and should be inserted under
%   'ContextModulation'.
%
%      'excitatory'
%      'inhibitory'
%      'threshold'
%      'gain'
%
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
%
%   available outputs: 
%       
%       Reinforced Hit Rate (1)
%       Reinforced FA Rate (2)
%       Probe Hit Rate (3)
%       Probe FA Rate (4)
%       Excitatory synaptic weights ([S S+ S-]) (5)
%       Inhibitory synaptic weights ([S S+ S-]) (6)
%
%--------------------------------------------------------------------------
%   Paper: Expression of task-knowledge during learning is context-dependent
%   Author: Tom Hindmarsh Sten, 2017
%   e-mail: hindmarsh.sten@nyu.edu
%   Release date: XXX
%   Version: 1.0.3
%   Code repository: http://froemkelab.med.nyu.edu/
%--------------------------------------------------------------------------

%% Retrieve parameters from input vector (params)

%learning rates
alpha = params(1);
alpha_NR = params(2);

%noise
sig = params(3);

%asymptote
kappa = params (4);

%initial weights
%     [S_I    S+_I    S-_I]
WI = [params(7) repmat(params(5),1,2)];
%     [S_I    S+_I    S-_I]
WE = [params(8) repmat(params(6),1,2)];

%inhibitory scaling parameter
c = params(9);

%% Non-Varying parameters

%reward magnitude
reward = 1;

%number of model runs
nloops = 5e1;
%% Load and organize behavioral data;
switch animal
    
    case 'average'
        
        load average_behavior_v2.mat
        
        %assign hit rates and smooth trajectory
        reinforcedhit = smooth(reinforced(:,2),5); probehit = smooth(probe(:,2),3);
        
        %assign false alarm rates and smooth trajectory
        reinforcedfa = smooth(reinforced(:,3),5); probefa = smooth(probe(:,3),3);
        
        
    otherwise
        
        load individual_behavior_v2.mat;
        
        %retreieve data for both contexts for the given animal
        reinforced = individual_behavior.(animal).reinforced;
        probe = individual_behavior.(animal).probe;
        
        %assign hit rates and smooth trajectory
        reinforcedhit = smooth(reinforced(:,2),5); probehit = smooth(probe(:,2),3);
        
        %assign false alarm rates and smooth trajectory
        reinforcedfa = smooth(reinforced(:,3),5); probefa = smooth(probe(:,3),3);
        
end
%% Some preparation
ReinforcedTrialBlocks = reinforced(:,1); ProbeTrialBlocks = probe(:,1);
EndTrial= max(ReinforcedTrialBlocks);

h = waitbar(0,'Running stochastic model with optimal parameters...');
absctr = 0;

%% initiate model
for j = 1:nloops;
    
    %Frequency of probe probing
    ProbeID = [1:EndTrial];
    
    %Number of reinforced trial blocks
    ReinforcedID = [1:EndTrial];
    
    %get initial weights (renew each model run).
    W_I = WI;
    W_E = WE;
    
    
    %% Generate a mix of target and foil tones
    
    trialcount = 100:100:EndTrial*100;
    stim = struct; stim(length(trialcount)).toneID = [];
    
    for i = 1:numel(trialcount);
        if i == 1;
            stim(i).toneID = randi(2,[trialcount(i),1]);
        else
            stim(i).toneID = randi(2,[trialcount(i)-trialcount(i-1),1]);
        end
    end
    
    %% Run through training sessions: block into 100s.
    
    ctr = 0;
    for session = 1:numel(stim);
        
        %preparation // house-keeping
        absctr = absctr +1;
        waitbar(absctr/(numel(stim)*nloops),h);
        ctr = ctr+1;
        hit = []; miss = []; cr = []; fa = [];
        
        %for each trial
        for trial = 1:numel(stim(session).toneID)
            
            %get toneID for given trial
            tone  = stim(session).toneID(trial);
            
            if tone == 1;
                x = [1 1 0]; %[s s+ s-], which neuron is activated
            else
                x = [1 0 1];
            end
            
            %get random number between zero and 1 with uniform probability;
            p = rand(1);
            
            switch ContextModulation
                
                case 'inhibitory'
                    
                    %probability of licking in the reinforced context
                    prob_lick = (1 + exp(-(W_E*x' - c*W_I*x')./sig)).^(-1);
                    
                case 'excitatory'
                    
                    %probability of licking in the reinforced context
                    prob_lick = (1 + exp(-(c*W_E*x' - W_I*x')./sig)).^(-1);
                    
                case 'threshold'
                    
                    %probability of licking in the reinforced context
                    prob_lick = (1 + exp(-(W_E*x' - W_I*x' + c)./sig)).^(-1);
                    
                case 'gain'
                    
                    %probability of licking in the reinforced context
                    prob_lick = (1 + exp(-(c*W_E*x' - c*W_I*x')./sig)).^(-1);
                    
                otherwise
                    error('Error. \nInput must be a valid scaling handle. The following is not a valid handle: %s',upper(type))
            end
            
            
            
            %outcome
            if p < prob_lick;
                y = 1; %Go
            else
                y = 0; %No-Go
            end
            
            %store trial outcome and update synaptic weights;
            if tone == 1 && y == 1; %Hit
                
                hit(trial) = 1;
                
                %update each synapse independently
                for synapse = 1:3;
                    if x(synapse) == 1;
                        
                        alpha_mult = alpha*W_E(synapse);
                        W_E(synapse) = W_E(synapse) + alpha_mult*(reward - kappa^-1*(W_E*x'-W_I*x'));
                        alpha_mult = alpha*W_I(synapse);
                        W_I(synapse) = W_I(synapse) - alpha_mult*(reward - kappa^-1*(W_E*x'-W_I*x'));
                        
                    end
                end
                
            elseif tone == 2 && y == 0; %correct reject; no synaptic update
                
                cr(trial) = 1;
                
            elseif tone == 1 && y == 0; %miss; no synaptic update
                
                miss(trial) = 1;
                
            elseif tone == 2 && y == 1; %False Alarm
                
                fa(trial) = 1;
                
                %update each synapse independently
                for synapse = 1:3;
                    if x(synapse) == 1;
                        
                        alpha_mult = alpha_NR*W_E(synapse);
                        W_E(synapse) = W_E(synapse) + alpha_mult*(-reward - kappa^-1*(W_E*x'-W_I*x'));
                        alpha_mult = alpha_NR*W_I(synapse);
                        W_I(synapse) = W_I(synapse) - alpha_mult*(-reward - kappa^-1*(W_E*x'-W_I*x'));
                        
                    end
                end
                
            end
            
        end;
        
        %store the synaptic weights from each training block
        WESTORE(:,ctr,j) = W_E;
        WISTORE(:,ctr,j) = W_I;
        
        %quantify performance in this reinforced session
        target_corr(ctr,j) = (sum(hit))/(sum(hit) + sum(miss));
        FA_rate(ctr,j) = (sum(fa))/(sum(cr) + sum(fa));
        
        %if we want to probe behavior, initiate 100 probe
        %trials; no synaptic updates
        if ismember(session,ProbeID); %if this is a probe trial
            hit = []; miss = []; cr = []; fa = [];
            
            for trial = 1:length(stim(session).toneID)
                %get toneID for given trial
                tone  = stim(session).toneID(trial);
                
                if tone == 1;
                    x = [1 1 0];
                else
                    x = [1 0 1];
                end
                
                %get random number between zero and 1 with uniform probability;
                p = rand(1);
                prob_lick = (1 + exp(-(W_E*x' - W_I*x')./sig)).^(-1);
                
                if p < prob_lick;
                    y = 1; %Go
                else
                    y = 0; %No-Go
                end
                
                %outcomes
                if tone == 1 && y == 1; %Hit
                    hit(trial) = 1;
                elseif tone == 2 && y == 0; %correct reject; no synaptic update
                    cr(trial) = 1;
                elseif tone == 1 && y == 0; %miss; no synaptic update
                    miss(trial) = 1;
                elseif tone == 2 && y == 1;
                    fa(trial) = 1;
                end
                
            end
            
            %quantify performance in this probe session
            target_corr_probe(ctr,j) = (sum(hit))/(sum(hit) + sum(miss));
            FA_rate_probe(ctr,j) = (sum(fa))/(sum(cr) + sum(fa));
            
        end
    end
end

switch visualization
    case 'on'
        target_corr = mean(target_corr,2);
        FA_rate = mean(FA_rate,2);
        
        target_corr_probe = nanmean(target_corr_probe,2);
        FA_rate_probe = nanmean(FA_rate_probe,2);
        nonpas = setdiff(ReinforcedID,ProbeID);
        target_corr_probe(nonpas) = NaN;
        FA_rate_probe(nonpas) = NaN;
        
        figure;
        set(gcf,'Position',[17,242,1220,463],'color','w')
        
        subplot(221);
        plot(ReinforcedID,target_corr,':g','LineWidth',1.5); %model hit
        hold on
        plot(ReinforcedID,FA_rate,':r','LineWidth',1.5); %model FA
        plot(reinforcedhit,'g','LineWidth',1.5);
        plot(reinforcedfa,'r','LineWidth',1.5);
        ylabel('Action Rate'); xlabel('Trial (100s)');
        legend('Model Hit','Model FA','Real Hit','Real FA');
        ylim([0 1]); xlim ([0 numel(stim)])
        title('Reinforced Performance');
        box off;
        
        subplot(222)
        plot([ProbeID],target_corr_probe(ProbeID),':g','LineWidth',1.5); %model hit
        hold on
        plot([ProbeID],FA_rate_probe(ProbeID),':r','LineWidth',1.5); %model FA
        plot(ProbeTrialBlocks,probehit,'og','LineWidth',1.5)
        plot(ProbeTrialBlocks,probefa,'or','LineWidth',1.5);
        legend('Model Hit','Model FA','Real Hit','Real FA');
        ylabel('Action Rate'); xlabel('Trial (100s)');
        ylim([0 1]); xlim ([0 numel(stim)])
        title('Probe Performance')
        box off;
        
        WESTORE = squeeze(nanmean(WESTORE,3));
        WISTORE = squeeze(nanmean(WISTORE,3));
        
        subplot(234);
        plot(WESTORE(2,:),'b'); hold on
        plot(WISTORE(2,:),'r');
        legend('S+/E','S+/I');
        title('S+ Weights');
        xlim ([0 numel(stim)])
        box off
        
        subplot(235);
        plot(WESTORE(3,:),'b'); hold on
        plot(WISTORE(3,:),'r');
        legend('S-/E','S-/I');
        title('S- Weights')
        xlim ([0 numel(stim)])
        box off
        
        subplot(236);
        plot(WESTORE(1,:),'b'); hold on
        plot(WISTORE(1,:),'r');
        legend('S/E','S/I');
        title('S Weights')
        xlim ([0 numel(stim)])
        box off
        
    otherwise
end
close(h);

varargout{1} = target_corr;
varargout{2} = FA_rate;
varargout{3} = target_corr_probe;
varargout{4} = FA_rate_probe;
varargout{5} = WESTORE;
varargout{6} = WISTORE;
end

