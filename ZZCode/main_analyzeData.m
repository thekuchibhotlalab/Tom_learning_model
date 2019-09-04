%% hello
clear;
nUnits = [3,5,10,15,20:10:100];
animal = 'average';
contextModulation = {'excitatory','inhibitory','threshold','gain'};
repetition = 20;
allModelCost = zeros(length(contextModulation),length(nUnits),repetition);
binaryModelCost = zeros(length(contextModulation),3);
contextColor = [ 0    0.4470    0.7410; ...
                 0.8500    0.3250    0.0980; ...
                 0.9290    0.6940    0.1250;...
                 0.4940    0.1840    0.5560];
figure;

for i = 1:length(contextModulation)
    load(['data/unitBinary_' contextModulation{i} '_' animal '.mat']);
    
    binaryModelCost(i,:) = allCost;
    scatter(zeros(1,3),allCost,10,contextColor(i,:),'filled'); hold on;
    minBinCost = min(binaryModelCost(i,:));
    for j = 1:length(nUnits)
        load(['data/unit' num2str(nUnits(j),'%03d') '_' contextModulation{i} '_' animal '.mat' ]);
        allModelCost(i,j,:) = allCost;
        scatter(repmat(nUnits(j),[1 repetition]),allCost,10,contextColor(i,:),'filled');
        hold on
    end
    minCost = min(allModelCost(i,:,:),[],3);
    %plot(nUnits, minCost,'Color',contextColor(i,:),'LineWidth',1.5);
    
    
    plot([0 nUnits],[minBinCost minCost],'Color',contextColor(i,:),'LineWidth',1.5);
end
set(gca,'yscale','log')

yticks(0:0.1:0.8)
yticklabels({'0','0.1','0.2','0.3','0.4','0.5','0.6','0.7','0.8'})

ylim([0.08 0.8])
xlim([0 100])

%% get single neuron figure
nUnit = 100;
repetition = 20;
load('data/unit100_inhibitory_average.mat');

[~,index] = min(allCost);

bestParam = allParam(:,index);

[target_corr, FA_rate, target_corr_probe, FA_rate_probe, WESTORE, WISTORE, W, act] = ...
    CircuitModel_Stochastic_ZZ(bestParam,nUnit,'average','inhibitory','on', repetition);

%concatWE = reshape(permute(WESTORE,[1 3 2]),size(WESTORE,1)*size(WESTORE,3),size(WESTORE,2));
%[projData,D,V] = doPCA(concatWE);

ETarg = repmat(reshape(act{1},[nUnit, 1, repetition]),[1 size(WESTORE,2) 1]) .* WESTORE;
concatETarg = reshape(permute(ETarg,[1 3 2]),size(ETarg,1)*size(ETarg,3),size(ETarg,2));

ITarg = repmat(reshape(act{1},[nUnit, 1, repetition]),[1 size(WISTORE,2) 1]) .* WISTORE;
concatITarg = reshape(permute(ITarg,[1 3 2]),size(ITarg,1)*size(ITarg,3),size(ITarg,2));

EFoil = repmat(reshape(act{2},[nUnit, 1, repetition]),[1 size(WESTORE,2) 1]) .* WESTORE;
concatEFoil = reshape(permute(EFoil,[1 3 2]),size(EFoil,1)*size(EFoil,3),size(EFoil,2));

IFoil = repmat(reshape(act{2},[nUnit, 1, repetition]),[1 size(WISTORE,2) 1]) .* WISTORE;
concatIFoil = reshape(permute(IFoil,[1 3 2]),size(IFoil,1)*size(IFoil,3),size(IFoil,2));
[projData,D,V] = doPCA({concatETarg,concatITarg,concatEFoil,concatIFoil});

function [projData, D, V] = doPCA(rawData)
    if iscell(rawData)
        processedData = cell2mat(rawData);
        trialTypeName = {'target E', 'target I', 'foil E', 'foil I'};
    else
        processedData = rawData;
    end
    processedData = (processedData - repmat(mean(processedData,2),1,size(processedData,2)));
    %data = (data - repmat(mean(data,2),1,size(data,2)))./repmat(std(data,0,2),1,size(data,2));
    covMatrix = cov(processedData');
    [D,V] = eig(covMatrix);
    [V,sortIndex] = sort(diag(V),'descend');
    V = V / sum(V);
    D = D(:,sortIndex);
    D(:,1) = -D(:,1);

    figure; 
    subplot(1,2,1);
    plot(0:length(V),[0;cumsum(V)],'-o','Linewidth',2);
    xlim([0 10]);ylim([0 1]);
    ylabel('%variance explained'); xlabel('#components');title('First 10 PC')  
    subplot(1,2,2);
    plot(0:length(V),[0;cumsum(V)],'Linewidth',2);
    xlim([0 200]);ylim([0 1]);
    ylabel('%variance explained'); xlabel('#components');title('First 200 PC')  

    for i = 1:size(D,2)
        if sum(D(:,i)) < 0 
            D(:,i) = -D(:,i);
        end
    end
    if iscell(rawData)
        for i = 1:length(rawData)
            projData{i} = D' * rawData{i};
        end           
    else
        projData = D' * processedData;
    end
    visualizeComponent = 6;
    figure;
    for i = 1:visualizeComponent
        subplot(2,3,i)
        if iscell(rawData)
            for j = 1:length(rawData)
                plot(projData{j}(i,:)); hold on;
            end
            legend(trialTypeName,'Location','northwest')
        else
            plot(projData(i,:))
        end
         title(['component' int2str(i) ' EV=' num2str(V(i)*100,'%.1f') '%'])
    end




end

