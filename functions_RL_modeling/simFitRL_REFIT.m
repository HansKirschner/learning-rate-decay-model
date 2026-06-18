function SimD = simFitRL_REFIT(IN,D,AF)
%keyboard

AF.nP          = sum(AF.Free); 

%% inital output
EV      = nan(length(D.reversal),1);
d       = nan(length(D.reversal),1);

%defaults 

Beta = AF.Defaults(1);

if ~AF.Free(3) %use one learn rate for everything 
    Alpha_Start = AF.Defaults(2);

elseif AF.Free(3) && ~AF.Free(4) %use seperate learn rates for confirmation or contradiction
    LRindex = [1 2 2 1];
    Alpha_Start(1) = AF.Defaults(2);
    Alpha_Start(2) = AF.Defaults(3);
   

elseif AF.Free(3) && AF.Free(4) %use seperate learn rates for everything
    
    Alpha_Start(1) = AF.Defaults(2);
    Alpha_Start(2) = AF.Defaults(3);
    Alpha_Start(1) = AF.Defaults(3);
    Alpha_Start(2) = AF.Defaults(4);
   
end

if AF.Free(6)
    Omega = AF.Order(6);
end

if AF.Free(7)
    Eta = AF.Defaults(7);
end


if AF.Free(8)
   Lambda = AF.Defaults(8);
end


%%
Beta        = IN(AF.Order(1));
if ~AF.Free(3) %use one learn rate for everything
    LRindex = [1 1 1 1]; %Index always points to LR Alha(1)   
    Alpha_Start = IN(AF.Order(2));
   

elseif AF.Free(3) && ~AF.Free(4) %use seperate learn rates for confirmation or contradiction
    LRindex = [1 2 2 1];
    Alpha_Start(1) = IN(AF.Order(2));
    Alpha_Start(2) = IN(AF.Order(3));
   

elseif AF.Free(3) && AF.Free(4) %use seperate learn rates for everything
    LRindex = [1 2 3 4];
    Alpha_Start(1) = IN(AF.Order(2));
    Alpha_Start(2) = IN(AF.Order(3));
    Alpha_Start(3) = IN(AF.Order(4));
    Alpha_Start(4) = IN(AF.Order(5));  
end

if AF.Free(6)
    Omega = IN(AF.Order(6));
end

if AF.Free(7)
   Eta = IN(AF.Order(7));
end
if AF.Free(8)
   Lambda = IN(AF.Order(8));
end


%OParams = AllParameter(AF.Free);
if AF.Free(6) && ~(AF.Free(7) == 1)
    OParams = [Beta Alpha_Start Omega];
elseif AF.Free(6) && (AF.Free(7) == 1)
    OParams = [Beta Alpha_Start Omega Eta];
else
   OParams = [Beta Alpha_Start];
end

if AF.Free(7) == 1
    Alpha   = nan(length(D.reversal),1);
end

if AF.Free(8) == 1
    OParams  = [OParams Lambda];
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
    if D.reversal(t)==1
        bt = 1;
        if AF.Free(6)
            EV(t) = Omega;
        else
            EV(t) = 0.5; %Reset expected values to chance when a stimulus is seen for the first time
        end
    end


    EV(EV>1) = 1;
    EV(EV<0) = 0;
    
    PC(t) = 1 ./ (1 + exp((-EV(t)+0.5) ./ Beta)); %softmax of expected value is choice probability

    %simulate choice
    ac(t) = 2 - (PC(t) > rand);

    if ac(t) == 1 && D.outcome(t) == 1
        D.event(t) = 1;
        D.choice(t) = 1;
        if t == 1
            D.BonusSum(t) = 10;
        else
            D.BonusSum(t) = D.BonusSum(t-1) + 10;
        end

    elseif ac(t) == 1 && D.outcome(t) == 0
        D.event(t) = 2;
        D.choice(t) = 1;
        if t == 1
            D.BonusSum(t) = -10;
        else
            D.BonusSum(t) = D.BonusSum(t-1) - 10;
        end
    elseif ac(t) == 2 && D.outcome(t) == 1
        D.event(t) = 3;
        D.choice(t) = 0;
        if t == 1
            D.BonusSum(t) = 0;
        else
            D.BonusSum(t) = D.BonusSum(t-1);
        end
    elseif ac(t) == 2 && D.outcome(t) == 0
        D.event(t) = 4;
        D.choice(t) = 0;
        if t == 1
            D.BonusSum(t) = 0;
        else
            D.BonusSum(t) = D.BonusSum(t-1);
        end
    end


    % Prediction error calculation
    d(t)        = D.outcome(t) - EV(t);


    if t < length(D.reversal)
        if ~AF.Free(7)
            % Update expected value
            if ~AF.Free(8)
                EV(t+1)     = EV(t) + Alpha_Start(LRindex(D.event(t))) * d(t);
            elseif AF.Free(8) && D.choice(t)==D.choice(t+1)

                EV(t+1)     = (EV(t)+Lambda) + Alpha_Start(LRindex(D.event(t))) * d(t);

            elseif AF.Free(8) && D.choice(t)~=D.choice(t+1)
                EV(t+1)     = EV(t) + Alpha_Start(LRindex(D.event(t))) * d(t);
            end

        elseif ~AF.PH && AF.Free(7)
            if ~AF.Free(8)
                % Update expected value
                Alpha(t)  = (Eta*1/bt) + (1-Eta)*Alpha_Start(LRindex(D.event(t)));
                EV(t+1)   = EV(t) + Alpha(t) * d(t);
            elseif AF.Free(8) && D.choice(t)==D.choice(t+1)
                Alpha(t)  = (Eta*1/bt) + (1-Eta)*Alpha_Start(LRindex(D.event(t)));
                EV(t+1)   = (EV(t)+Lambda) + Alpha(t) * d(t);
            elseif AF.Free(8) && D.choice(t)~=D.choice(t+1)
                Alpha(t)  = (Eta*1/bt) + (1-Eta)*Alpha_Start(LRindex(D.event(t)));
                EV(t+1)   = EV(t) + Alpha(t) * d(t);
            end


        elseif AF.PH && AF.Free(7)

            Alpha(t)  = Eta * abs(d(t)) + (1-Eta)*Alpha_Start(LRindex(D.event(t)));
            EV(t+1)   = EV(t) + Alpha(t) * d(t);
        end
    end

    bt = bt+1;
end


%%
SimD.reversal   = D.reversal;
SimD.outcome    = D.outcome;
SimD.choice     = D.choice;
SimD.event      = D.event;
SimD.Nch        = D.Nch;
SimD.Bonus      = D.BonusSum(end);
SimD.PC         = PC;


return