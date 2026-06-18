function [SumAIC, SumBIC, SumLL, iBIC, bic, aic] = PlotFit(AF,FitPara,nLL,D)
%keyboard
%%
%Remove fixed parameters
if sum(AF.Free) ~= numel(AF.Order)
    AF.HyperPriors = AF.HyperPriors(logical(AF.Free),:);
    AF.ParamNames = AF.ParamNames(logical(AF.Free));
    FitPara = FitPara(logical(AF.Free),:);
end
[aic, bic] = aicbic(-nLL,AF.nP,[D.Nch]);
SumAIC = sum(aic)
SumBIC = sum(bic)
SumLL  = sum(nLL)

% compute iBIC
Nsj = length(D); sjind=1:Nsj;
Np = AF.nP;
reg=cell(Np,1);
Xreg = repmat(eye(Np),[1 1 Nsj]);
Nreg=0; 
for j=1:Np
	if size(reg{j},2)==Nsj; reg{j} = reg{j}';end
	for k=1:size(reg{j},2)
		Nreg = Nreg+1;
		Xreg(j,Np+Nreg,:) = reg{j}(:,k); 
	end
end
Nch=0; for sj=sjind; Nch = Nch + D(sj).Nch;end
iBIC =  -2*(sum(-nLL) - 1/2*(2*Np+Nreg)*log(Nch));




if AF.PlotOn
    figure(1);clf;
    z = 2; s = AF.nP;
    for c = 1 : AF.nP
        subplot(z,s,c)
        p = c;
        xv = AF.HyperPriors(p,1):0.001:AF.HyperPriors(p,2);
        FV = pdf(AF.prior_distributions(p),xv);
        FV(FV<AF.CutP) = AF.CutP;
        plot(xv,-log(FV))
        ylabel('nLL')
        yyaxis right
        plot(xv,FV)
        ylabel('PDF')
    end
    for c = 1 : AF.nP
        subplot(z,s,AF.nP+c)
        p = c;
        histogram(FitPara(c,:)',12);
        title(['Distribution for: ' AF.ParamNames{c}])
        set(gca, 'XLim', [AF.HyperPriors(p,1) AF.HyperPriors(p,2)]);
        xlabel(AF.ParamNames{c}); ylabel('n sbj');
    end
    
    %%
    figure(2);clf;
    [r,p]=corr(FitPara')
    imagesc(r)
    set(gca,'XTickLabel', AF.ParamNames, 'YTick', 1:AF.nP, 'YTickLabel', AF.ParamNames)
    colorbar
    
end
%%
return