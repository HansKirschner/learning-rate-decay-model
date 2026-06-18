clear;
path('functions_RL_modeling/',path);
% load data
load('data/exampleData.mat');
for i = 1:length(D)
    Index1 = find(D(i).Stimno ==1 & D(i).reversal ==1);
    D(i).reversal(Index1(end))=0;
end

%%
% Model 6 - LR for confirmation/disconferimation + starting EV scaling + LR decay + stickyness
AF.HyperPriors = [0 2;       0 1;     0 1;      0 1;     0 1;       0 1;       0 1;  -1 1];
AF.ParamNames  = {'beta'   'alpha1'   'alpha2' 'alpha3' 'alpha4'    'omega'    'eta' 'lambda'}; 
AF.Defaults    = [0.1        0.4      0.4       0.4      0.4        0.5        .2    .1];   %default parameters if no parameter is present
AF.Order       = [1          2        3         4        5          6          7     8];   %maps onto model order of parameters
AF.Free        = [1          1        1         0        0          1          1     1];
AF.Cut         = 1e-4; % how estimated
AF.CutP        = 1e-7;  %Determines power of the priors - how estimated - also currently not applied in FitRL?
FitF           = @FitRL_REFIT;
AF.ModelName   = char(FitF);
AF.PlotOn      = 0; % if set, plots Fit information and parameter distribution
AF.DoPrior     = 1;
AF.PH          = 0;
AF.nP          = sum(AF.Free);

AF.prior_distributions = [makedist('gamma', 1.4,0.1) makedist('gamma', 1.4,0.2) makedist('gamma', 1.4,0.2) makedist('normal', .5,.3) makedist('Uniform', 0,1)  makedist('normal', .12,.4)];
AF.prior_functions = @(x1) ([pdf(AF.prior_distributions(1),x1(1)) pdf(AF.prior_distributions(2),x1(2)) pdf(AF.prior_distributions(3),x1(3)) pdf(AF.prior_distributions(4),x1(4)) pdf(AF.prior_distributions(5),x1(5)) pdf(AF.prior_distributions(6),x1(6))]); 

oldOpts = optimset('fmincon');
options=optimset(oldOpts, 'maxFunEvals', 1000000, 'MaxIter', 1000000,'Display', 'off');
reps = 20;

clear FitPara nLL
for s = 1 : length(D)
    tic
    fprintf('Modelling subject %d / %d',s,length(D));

    minNegLL=inf;
    params_fit = [];
    i = 1;
    while i <= reps

        if i == 1
            fminconOutcome.fitParams  = fmincon(@(x)FitF(x,D(s),AF), AF.Defaults, [], [], [], [], AF.HyperPriors(:,1), AF.HyperPriors(:,2),[],options);
            fminconOutcome.negLogLike = FitF(fminconOutcome.fitParams,D(s),AF);
        else
            OK = 0; count = 0;
            while OK == 0
                try
                    AF.newDefaults = unifrnd(AF.HyperPriors(:,1),AF.HyperPriors(:,2))';
                    fminconOutcome.fitParams  = fmincon(@(x)FitF(x,D(s),AF), AF.newDefaults, [], [], [], [], AF.HyperPriors(:,1), AF.HyperPriors(:,2),[],options);
                    fminconOutcome.negLogLike = FitF(fminconOutcome.fitParams,D(s),AF);
                    OK = 1;
                catch
                    OK = 0;

                end
            end
        end

        if fminconOutcome.negLogLike<minNegLL
            bestOutput=fminconOutcome;
            minNegLL=bestOutput.negLogLike;
            FitPara(:,s) = bestOutput.fitParams';
        end

        i = i+1;
    end

    [nLL(s), mo] = FitF(FitPara(:,s),D(s),AF);

    % Add Model parameters to behavioral data
    D(s).RPE   = mo.decisions(:,1)';
    D(s).EV    = mo.decisions(:,2)';
    D(s).PC1   = mo.decisions(:,3)';
    D(s).lr    = mo.LR';

    elapsedTime = toc;
    fprintf('(%.02fsec)\n',elapsedTime);
end

save('data/FitM6_SampleData',"D")
save('data/FitM6_Parameters_SampleData',"FitPara")
ert

%% Plot model predictions
clear

load('data/FitM6_SampleData.mat')
load('data/FitM6_Parameters_SampleData')
load('data/EndPoints_SampleData.mat')

n = size(D,2);
SymStage = [];
for c =  1: n
    US = D(c).SymStage'; If = 1; AllI = [];
    for c2 = 2 : length(US)
        if US(c2)<=US(c2-1) %end of stage
            if US(c2-1)>19
                AllI = [AllI If:c2-1];
            end
            If = c2;
        end
    end
    length(AllI);
    AllIAll(c) = {AllI};
    SymStage = [SymStage D(c).SymStage(AllI)' ];
end
SymStage


% Change all fieldnames into equally size arrays
FN = fieldnames(D); FN(1:2) = [];
for c = [1:length(FN)]
    eval([FN{c} ' = [];'])
    for c2 = 1 : n
        eval([FN{c} ' = [' FN{c} ' D(c2).(FN{c})(AllIAll{c2})''];'])
    end
end

%% Plot choice behavior and chocie predictions
set(0,'defaultAxesFontSize',20);
figure(1);clf;
shadedErrBar([], AGF_running_average(mean(choice(:,:),2),2,2), AGF_running_average(se(choice(:,:)',0.99),2,2),'k'); hold;
plot(mean(Prob,2)/100,'-k',LineWidth=2); 
plot(AGF_running_average(mean(PC1(:,:),2),2,2),LineWidth=2); 
legend({'choice' 'GT'})
ylim([-.05 1.05]);gridxy(find(reversal(:,1)==1));
xlim([1 length(reversal)])
ylabel('Choices'); %xlabel('Tr Nr'); 
title('Example Data');
h= gcf;
print2pdf(fullfile(pwd,'Figures',['Choice_BHV.pdf']),h,300)

%% plot learning rate trjectories

LR_confirmation    = nan(21,length(D));
LR_disconfirmation = nan(21,length(D));


for t = 1:21
    for Vp =1:length(D)
    LR_confirmation(t,Vp)    = (FitPara(7,Vp)*1/t) + ((1-FitPara(7,Vp))*FitPara(2,Vp));
    LR_disconfirmation(t,Vp) = (FitPara(7,Vp)*1/t) + ((1-FitPara(7,Vp))*FitPara(3,Vp));
    end
end


set(0,'defaultAxesFontSize',20);
figure(2);clf;
plot(mean(LR_confirmation(:,:),2),LineWidth=2);hold on
plot(mean(LR_disconfirmation(:,:),2),LineWidth=2)
legend({'Confirmation' 'Disconfirmation'})
xlim([0 22]);ylabel('LR');ylim([0 .6]);

%% Correlation Parameters and Endpoints

[rOmega,pOmega]= corr(EndPOints',FitPara(6,:)');

[rEta,pEta]= corr(EndPOints',FitPara(7,:)');

[rAlphaDis,pAlphaDis]= corr(EndPOints',FitPara(3,:)');


set(0,'defaultAxesFontSize',20);
figure(3);clf;
subplot(131)
scatter(FitPara(6,:)',EndPOints',100,'m','filled')
refline
xlabel('\omega'); ylabel('total points earned')
title(['r = ' num2str(round2(rOmega,0.01)) ', p = ' num2str(round2(pOmega,0.001))])
subplot(132)
scatter(FitPara(7,:),EndPOints',100,'g','filled')
xlabel('\eta');
refline
title(['r = ' num2str(round2(rEta,0.01)) ', p = ' num2str(round2(pEta,0.001))])
subplot(133)
scatter(FitPara(3,:)',EndPOints',100,'y','filled')
refline
xlabel('\alpha disconf'); ylabel('total points earned')
title(['r = ' num2str(round2(rAlphaDis,0.01)) ', p = ' num2str(round2(pAlphaDis,0.001))])

h=gcf;
print2pdf(fullfile(pwd,'Figures',['Bias_and_Performance.pdf']),h,300)
