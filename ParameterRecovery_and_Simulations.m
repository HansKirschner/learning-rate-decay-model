%% This script runs the parameter Recovery for M6

%% 1. Set-up and data loading
clear

path('functions_RL_modeling/',path);
load('data/FitM6_SampleData.mat')
savepath = 'data/Simulated_Data/';
load('data/FitM6_Parameters_SampleData.mat')

%% 2. set-up Full Model and run simulation
%
% Model 6 - full model + starting EV scaling + LR decay + stickyness
AF.HyperPriors = [0 .4;       0 1;     0 1;      0 1;     0 1;       0 1;       0 1;  -1 1];
AF.ParamNames  = {'beta'   'alpha1'   'alpha2' 'alpha3' 'alpha4'    'omega'    'eta'  'lamba'}; 
AF.Defaults    = [0.1        0.4      0.4       0.4      0.4        0.5        .2     .1];   %default parameters if no parameter is present
AF.Order       = [1          2        3         4        5          6          7      8];   %maps onto model order of parameters
AF.Free        = [1          1        1         0        0          1          1      1];
AF.Cut         = 1e-4; % how estimated
AF.CutP        = 1e-7;  %Determines power of the priors - how estimated - also currently not applied in FitRL?
FitF           = @FitRL_REFIT;
AF.ModelName   = char(FitF);
AF.PlotOn      = 0; % if set, plots Fit information and parameter distribution
AF.DoPrior     = 1;
AF.nP          = sum(AF.Free);
AF.PH          = 0;

AF.prior_distributions = [makedist('gamma', 1.4,0.1) makedist('gamma', 1.4,0.2) makedist('gamma', 1.4,0.2) makedist('normal', .5,.3) makedist('Uniform', 0,1)  makedist('normal', .12,.4)];
AF.prior_functions = @(x1) ([pdf(AF.prior_distributions(1),x1(1)) pdf(AF.prior_distributions(2),x1(2)) pdf(AF.prior_distributions(3),x1(3)) pdf(AF.prior_distributions(4),x1(4)) pdf(AF.prior_distributions(5),x1(5)) pdf(AF.prior_distributions(6),x1(6))]); 

oldOpts = optimset('fmincon');
options=optimset(oldOpts, 'maxFunEvals', 1000000, 'MaxIter', 1000000,'Display', 'off');

reps       = 20;
n_search   = 10; % number of iterations for recovery and correlation analyses (just a quick check here, should be a reasonable number (e.g. 500 simulations in total)


GT_parameters(:,1) = FitPara(1,:);
GT_parameters(:,2) = FitPara(2,:);
GT_parameters(:,3) = FitPara(3,:);
GT_parameters(:,4) = FitPara(4,:);
GT_parameters(:,5) = FitPara(5,:);
GT_parameters(:,6) = FitPara(6,:);
GT_parameters(:,7) = FitPara(7,:);
GT_parameters(:,8) = FitPara(8,:);


FitAll              = nan(length(D)*n_search,length(AF.Free));
GT_parameters_out   = nan(length(D)*n_search,length(AF.Free));
BonusSum            = nan(length(D),n_search);
         
AF.DoPrior   = 1;

for n = 1 : length(D)
    
    tic;
    fprintf('Starting simulations for subject %d / %d ...\n',n, length(D));
    
    Choices = [];PC=[];
    
    for s = 1 : n_search

        %simulate data with the parameter set
        [SimD]   =  simFitRL_REFIT(GT_parameters(n,:),D(n),AF); %use model with ground truth parameters & outcomes (same for all subjects)

        Choices(s,:) = SimD.choice;
        PC(s,:)      = SimD.PC;

        minNegLL=inf;
        params_fit = [];
        i = 1;
        while i <= reps

            if i == 1
                fminconOutcome.fitParams  = fmincon(@(x)FitF(x,SimD,AF), AF.Defaults, [], [], [], [], AF.HyperPriors(:,1), AF.HyperPriors(:,2),[],options);
                fminconOutcome.negLogLike = FitF(fminconOutcome.fitParams,SimD,AF);
            else
                OK = 0;
                while OK == 0
                    try
                        AF.newDefaults = unifrnd(AF.HyperPriors(:,1),AF.HyperPriors(:,2))';
                        fminconOutcome.fitParams  =  fmincon(@(x)FitF(x,SimD,AF), AF.Defaults, [], [], [], [], AF.HyperPriors(:,1), AF.HyperPriors(:,2),[],options);
                        fminconOutcome.negLogLike = FitF(fminconOutcome.fitParams,SimD,AF);
                        OK = 1;
                    catch
                        OK = 0;

                    end
                end
            end

            if fminconOutcome.negLogLike<minNegLL
                bestOutput=fminconOutcome;
                minNegLL=bestOutput.negLogLike;
                FitAll((n-1)*n_search+s,:) = bestOutput.fitParams';
            end

            i = i+1;
        end
        
        if s == .5*n_search
            fprintf(' %d / %d of simulations done... \n',s,n_search);
        end


    end

    for GT_out = (n-1)*n_search+1:n*n_search
        GT_parameters_out(GT_out,:) = GT_parameters(n,:);
    end
    
    
    elapsedTime = toc;
    fprintf('total time for subject %d - %.02fsec \n',n,elapsedTime);
    
    D(n).SimChoices = Choices;
    D(n).SimPC      = PC;

end

save([savepath,'ParamRecovery_m6.mat'],'GT_parameters_out','FitAll')


ert

save([savepath,'ParamRecovery_m6.mat'],'GT_parameters_out','FitAll')


%% 3. plot results
load('data/Simulated_Data/ParamRecovery_m6.mat')

n_search   = 10;      % number of iterations for recovery and correlation analyses
n          = 12;      % number of subjects

set(0,'defaultAxesFontSize',20);

for c = 1 : size(GT_parameters_out,2)
    [r1,~] = corr([GT_parameters_out(:,c) FitAll(:,c)]); %correlation between ground truth parameters and parameters estimated from choices
    [GTr(c,:)] = regress(GT_parameters_out(:,c),[ones(n*n_search,1) FitAll(:,c)]);    
    GTc(1,c) = round(r1(2)*100)/100;
end

minB = min(min(FitAll(:,1),GT_parameters_out(:,1)));
maxB = max(max(FitAll(:,1),GT_parameters_out(:,1)));
%swap axes....
h = figure(1);clf;
plot(GT_parameters_out(:,1),FitAll(:,1),'o'); hold on
set(gca, 'ylim',[minB maxB],'xlim',[minB maxB])
text(minB+minB/2 ,maxB-minB*2,{['r = ' num2str(GTc(1))],['{\beta_{0}} = ' num2str(round(GTr(1,1)*100)/100)],['{\beta_{1}} = ' num2str(round(GTr(1,2)*100)/100)]},'FontSize',20)
ylabel('estimated \beta', 'FontSize',20);
xlabel('simulated \beta', 'FontSize',20);
fig = gcf;
fig.PaperUnits = 'centimeters';
fig.PaperPosition = [0 0 9 9];
print2pdf(fullfile(pwd,'Figures',['beta_sim_M6.pdf']),h,300)


minA1 = min(min(FitAll(:,2),GT_parameters_out(:,2)));
maxA1 = max(max(FitAll(:,2),GT_parameters_out(:,2)));
h = figure(2);clf;
plot(GT_parameters_out(:,2),FitAll(:,2),'o'); hold on
%set(gca, 'ylim',[min(FitPara(2,:)) max(FitPara(2,:))],'xlim',[min(FitPara(2,:)) max(FitPara(2,:))])
set(gca, 'ylim',[minA1 maxA1],'xlim',[minA1 maxA1])
text(minA1+minA1/2 ,maxA1-maxA1/2,{['r = ' num2str(GTc(2))],['{\beta_{0}} = ' num2str(round(GTr(2,1)*100)/100)],['{\beta_{1}} = ' num2str(round(GTr(2,2)*100)/100)]},'FontSize',20)
xlabel('estimated \alpha confirmation', 'FontSize',20);
ylabel('simulated \alpha confirmation', 'FontSize',20);
fig = gcf;
fig.PaperUnits = 'centimeters';
fig.PaperPosition = [0 0 9 9];
% print('beta_sim','-dpdf')
print2pdf(fullfile(pwd,'Figures',['alpha1_conf_sim_M6.pdf']),h,300)

minA2 = min(min(FitAll(:,3),GT_parameters_out(:,3)));
maxA2 = max(max(FitAll(:,3),GT_parameters_out(:,3)));
h = figure(3);clf;
plot(GT_parameters_out(:,3),FitAll(:,3),'o'); hold on
%set(gca, 'ylim',[min(FitPara(2,:)) max(FitPara(2,:))],'xlim',[min(FitPara(2,:)) max(FitPara(2,:))])
set(gca, 'ylim',[minA2 maxA2],'xlim',[minA2 maxA2])
text(minA2+minA2/2 ,maxA2-maxA2/2,{['r = ' num2str(GTc(3))],['{\beta_{0}} = ' num2str(round(GTr(3,1)*100)/100)],['{\beta_{1}} = ' num2str(round(GTr(3,2)*100)/100)]},'FontSize',20)
xlabel('estimated \alpha disconfirmation', 'FontSize',20);
ylabel('simulated \alpha disconfirmation', 'FontSize',20);
fig = gcf;
fig.PaperUnits = 'centimeters';
fig.PaperPosition = [0 0 9 9];
% print('beta_sim','-dpdf')
print2pdf(fullfile(pwd,'Figures',['alpha1_disconf_sim_M6.pdf']),h,300)

minEta = min(min(FitAll(:,6),GT_parameters_out(:,6)));
maxEta = max(max(FitAll(:,6),GT_parameters_out(:,6)));
h = figure(4);clf;
plot(GT_parameters_out(:,6),FitAll(:,6),'o'); hold on
%set(gca, 'ylim',[min(FitPara(2,:)) max(FitPara(2,:))],'xlim',[min(FitPara(2,:)) max(FitPara(2,:))])
set(gca, 'ylim',[minEta  maxEta],'xlim',[minEta  maxEta])
text(minEta+minEta/2 ,maxEta-maxEta/2,{['r = ' num2str(GTc(6))],['{\beta_{0}} = ' num2str(round(GTr(6,1)*100)/100)],['{\beta_{1}} = ' num2str(round(GTr(6,2)*100)/100)]},'FontSize',20)
xlabel('estimated \eta', 'FontSize',20);
ylabel('simulated \eta', 'FontSize',20);
fig = gcf;
fig.PaperUnits = 'centimeters';
fig.PaperPosition = [0 0 9 9];
% print('beta_sim','-dpdf')
print2pdf(fullfile(pwd,'Figures',['omega_sim_M6.pdf']),h,300)

minOmega = min(min(FitAll(:,7),GT_parameters_out(:,7)));
maxOmega = max(max(FitAll(:,7),GT_parameters_out(:,7)));
h = figure(5);clf;
plot(GT_parameters_out(:,end),FitAll(:,end),'o'); hold on
%set(gca, 'ylim',[min(FitPara(2,:)) max(FitPara(2,:))],'xlim',[min(FitPara(2,:)) max(FitPara(2,:))])
set(gca, 'ylim',[minOmega maxOmega],'xlim',[minOmega maxOmega])
text(minOmega+minOmega/2 ,maxOmega-maxOmega/2,{['r = ' num2str(GTc(7))],['{\beta_{0}} = ' num2str(round(GTr(7,1)*100)/100)],['{\beta_{1}} = ' num2str(round(GTr(7,2)*100)/100)]},'FontSize',20)
xlabel('estimated \omega', 'FontSize',20);
ylabel('simulated \omega', 'FontSize',20);
fig = gcf;
fig.PaperUnits = 'centimeters';
fig.PaperPosition = [0 0 9 9];
% print('beta_sim','-dpdf')
print2pdf(fullfile(pwd,'Figures',['eta_sim_M6.pdf']),h,300)

minLamda = min(min(FitAll(:,8),GT_parameters_out(:,8)));
maxLamda = max(max(FitAll(:,8),GT_parameters_out(:,8)));
h = figure(6);clf;
plot(GT_parameters_out(:,end),FitAll(:,end),'o'); hold on
%set(gca, 'ylim',[min(FitPara(2,:)) max(FitPara(2,:))],'xlim',[min(FitPara(2,:)) max(FitPara(2,:))])
set(gca, 'ylim',[minLamda maxLamda],'xlim',[minLamda maxLamda])
text(minLamda+minLamda/2 ,maxLamda-maxLamda/2,{['r = ' num2str(GTc(8))],['{\lambda_{0}} = ' num2str(round(GTr(8,1)*100)/100)],['{\lambda_{1}} = ' num2str(round(GTr(8,2)*100)/100)]},'FontSize',20)
xlabel('estimated \lambda', 'FontSize',20);
ylabel('simulated \lambda', 'FontSize',20);
fig = gcf;
fig.PaperUnits = 'centimeters';
fig.PaperPosition = [0 0 9 9];
% print('beta_sim','-dpdf')
print2pdf(fullfile(pwd,'Figures',['lambda_sim_M6.pdf']),h,300)







