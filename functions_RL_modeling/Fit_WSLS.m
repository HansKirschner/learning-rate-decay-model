function [nLL, mo] = Fit_WSLS(IN,D,AF)
%keyboard

%% inital output
Switch  = nan(length(D.reversal),1);
Feedbackindex = [1 2 2 1];



if ~AF.Free(1) && ~AF.Free(2) %use eta of 0 
    
    Eta1 = 0;
    Eta2 = 0;

elseif AF.Free(1) && AF.Free(2) %use seperate learn rates for confirmation or contradiction
    
    Eta1 = AF.Defaults(1);
    Eta2 = AF.Defaults(2);
     
end

SwitchConf      = AF.Defaults(3);
SwitchDisconf   = AF.Defaults(4);

%%

if AF.Free(1) && AF.Free(2)

    Eta1 = IN(AF.Order(1));
    Eta2 = IN(AF.Order(2));
end

SwitchConf      = IN(AF.Order(3));
SwitchDisconf   = IN(AF.Order(4));


%OParams = AllParameter(AF.Free);
if AF.Free(1) && AF.Free(2)
    OParams = [Eta1 Eta2 SwitchConf SwitchDisconf];
else
   OParams = [SwitchConf SwitchDisconf];
end


%% calculate priors
if AF.DoPrior
    Prior               = AF.prior_functions(OParams);
    Prior(Prior<AF.Cut) = AF.Cut;
    Prior               = -log(Prior); %has been change to only incluce priors for free parameters
else
    Prior = 0;
end


%% Model
bt = 0;

for t = 1 : length(D.reversal)
    
    if D.reversal(t)
          bt = 1;
    end

    if Feedbackindex(D.event(t)) == 1 %confirmational feedback

        Switch(t) =  Eta1/bt + (1-Eta1)*SwitchConf;
   
    elseif Feedbackindex(D.event(t)) == 2 %disconfirmational feedback
        
        Switch(t) =  Eta2/bt + (1-Eta2)*SwitchDisconf;
   
    end
    
     bt = bt+1;
end



%%
% Probability

%bring Switchs in borders
Switch(Switch>1) = 1;
Switch(Switch<0) = 0;

%get choice probabilty and choice --> for the WSLS model this is redundant
%´
% PC                  = 1 ./ (1 + exp((-EV+0.5) ./ Beta)); %softmax of expected value is choice probability - this beta to one ad define the other stuff above

mo.decisions        = Switch;

InputLog = nanrem(1-abs(Switch-make_long([D.RS])));

%
nLL               = -sum(log(InputLog)) + sum(Prior);


mo.OParams        = OParams;
mo.CallParms      = IN;
mo.Priors         = Prior;


return