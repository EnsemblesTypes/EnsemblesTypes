function  [Results, Points,Weighting,MaxentStore] = Make_Ensembles(Results,Points,Weighting,Parameters,Set2,~)
%model combine module, based on Combination
warning off
%Reset parameters
Parameters.make_log = 0;
Parameters.ensemble = 1;
MaxentStore = [];
%HIMaxentStore = [];
%% All preparations
% Set data_set
mod_tmp = [];
val_tmp = [];
mod_tmp2 = [];
val_tmp2 = [];
for i=1:1: Parameters.data_set_max
    model_Values = Points.Models.(genvarname(char(Parameters.SetNames(i))));
    model_ValuesS2 = Set2.Models.(genvarname(char(Parameters.SetNames(i))));
    vali_Values = Points.Validation.(genvarname(char(Parameters.SetNames(i))));
    vali_ValuesS2 =  Set2.Validation.(genvarname(char(Parameters.SetNames(i))));
    Parameters.vali_nanS1(i,1) = length(find(isnan(vali_Values)));
    Parameters.vali_nanS2(i,1) = length(find(isnan(vali_ValuesS2)));
    mod_tmp = [mod_tmp,model_Values]; %#ok<*AGROW>
    mod_tmp2 = [mod_tmp2,model_ValuesS2];
    val_tmp = [val_tmp,vali_Values];
    val_tmp2 = [val_tmp2,vali_ValuesS2];
end
model_values_store = mod_tmp;
model_values_Set2 = mod_tmp2;
validation_store =  val_tmp;
validation_S2 =  val_tmp2;
Result_model_txt = {'Datapoints';'RHO';'PVAL';'mean_deviance'};
Results.Ensemble = dataset([1;1;1;1],'ObsNames', Result_model_txt,'Varnames',char(Parameters.Ensemble_Names(1)));
clear i list model_Values clear model_ValuesS2 val_tmp* mod_tmp* Result_model_txt
%% Ensemble stats
for Ensemble = 1:Parameters.Nr_ensembles
    if Ensemble < Parameters.CutOffType
        Sizes = Parameters.Sizes(Parameters.Training);
        Parameters.valiNaN= (find(min(Parameters.vali_nanS1)));
        if isempty(Parameters.valiNaN) == 1
            Parameters.valiNaN = 1;
        end
        [model_values, name,Weighting,MaxentStore,validation_set] = Actual_Ensembles(model_values_store, model_values_Set2,Weighting,...
            Ensemble,Parameters,Results,MaxentStore,validation_store,validation_S2,Sizes);
    else
        Sizes = Parameters.Sizes(Parameters.Training);
        Parameters.valiNaN= find(min(Parameters.vali_nanS2));
        if isempty(Parameters.valiNaN) == 1
            Parameters.valiNaN = 1;
        end
        [model_values, name,Weighting,~,validation_set] = Actual_Ensembles(model_values_store, model_values_Set2, Weighting,...
            Ensemble,Parameters,Results,MaxentStore,validation_store,validation_S2,Sizes);
        Sizes = Parameters.Sizes(Parameters.Validator);
    end
    testArray = [validation_set,model_values];
    [Outputs] = Accuracy_statistics(testArray,Parameters,Sizes);% run statistics
    Results.Ensemble.(genvarname(char(name))) = [Outputs.datapoints;Outputs.RHO;Outputs.PVAL;Outputs.mean_double_deviation]; %output
    Results.Models.Data_set((Parameters.data_set_max+Ensemble),1) = name;
    Results.Models.Datapoints((Parameters.data_set_max+Ensemble),1) = Outputs.datapoints;
    Results.Models.RHO((Parameters.data_set_max+Ensemble),1) = Outputs.RHO;
    Results.Models.PVal((Parameters.data_set_max+Ensemble),1) =  Outputs.PVAL;
    Results.Models.Inversed_deviance((Parameters.data_set_max+Ensemble),1) =   Outputs.mean_double_deviation;
    %if Ensemble < Parameters.CutOffType
    Points.Validation.(genvarname(char(name)))  = Outputs.xes;
    Points.Models.(genvarname(char(name)))  = Outputs.yes;
    Points.Deviation.(genvarname(char(name))) = Outputs.deviation_point;
    %end
    clear testArray Outputs name
end
end
%%
function [model_values, name,Weighting,MaxentStore,validation_set] = Actual_Ensembles(model_values_store, model_values_S2, Weighting,...
    Ensemble,Parameters,Results,MaxentStore,validation_store,validation_S2,Sizes)
name = Parameters.Ensemble_Names(Ensemble);
str = sprintf('Ensemble running = %s ',char(name));
disp(str)
clear InvModels new_weights model_values
validation_set(:,1) = validation_store(:,Parameters.valiNaN);
%% Reference Best models
if Ensemble == 1
    nrM = find(Results.Models.Inversed_deviance==max(Results.Models.Inversed_deviance));
    model_values(:,1) = model_values_store(:,nrM);
    clear validation_set
    validation_set(:,1) = validation_store(:,nrM);
    Weighting.BestDev= nrM;
elseif Ensemble == 2
    nrM = find(Results.Models.RHO(1:(size(Results.Models.RHO)-1))==...
        max( Results.Models.RHO(1:(size(Results.Models.RHO)-1))));
    model_values(:,1) = model_values_store(:,nrM);
    clear validation_set
    validation_set(:,1) = validation_store(:,nrM);
    Weighting.BestRho= nrM;
    %% Uninformed Ensembles
elseif Ensemble == 3
    % Mean
    model_values(:,1) = (nanmean(model_values_store,2)); %The recalculation method
elseif  Ensemble == 4
    % Median stats
    model_values(:,1) = (nanmedian(model_values_store,2)); %The recalculation method
elseif  Ensemble == 5
    %PCA weights
    InvModels = model_values_store;
    [InvModels, ~] = ModelClean(InvModels, []);
    [coefs,~,~,~] = princomp(InvModels);
    PCA_weights = coefs(:,1);
    InvModels = model_values_store; % Restore full dataset
    [model_values,~,new_weights] = Weightin_algo(InvModels,PCA_weights,Parameters);
    Weighting.PCA(:,1) = new_weights;
    clear PCA_weights
elseif Ensemble == 6 % Corrcoef
    disp('      Correlation coefficient among ouputs')
    InvModels = model_values_store;
    [InvModels, ~] = ModelClean(InvModels, []);
    NewWeights = mean(corrcoef(InvModels));
    Weights = ShapeWeights(NewWeights,Parameters);
    InvModels = model_values_store; % Restore full dataset
    [model_values,~,new_weights] = Weightin_algo(InvModels,Weights,Parameters);
    Weighting.CorCoef(:,1) = Weights;
    clear Weights NewWeights
elseif Ensemble == 7 %Regression to the median
    InvModels = model_values_store;
    Predicted(:,1) = (nanmedian(InvModels,2));
    [InvModels, Predicted] = ModelClean(InvModels, Predicted);
    disp('      Regressing fit to median')
    [Median_Models,~]  = RegresAlgo(Parameters,InvModels,Predicted);
    InvModels = model_values_store; % Restore full dataset
    [model_values,~,new_weights] = Weightin_algo(InvModels,Median_Models,Parameters);
    Weighting.RegresstoMedian(:,1) = new_weights; %Store to compare
    clear  Median_Models
elseif Ensemble == 8 %Leave one out regression
    InvModels = model_values_store;
    %Step 1: Estimate all one to all fits, restrict weights and calcuate deviance
    step = 1;
    nrM = size(model_values_store,2);
    for i = 1:1:nrM
        InvModels = model_values_store;
        take_array = RemvOne(nrM,i);
        Predicted(:,1) = InvModels(:,i);
        InvModels =  InvModels(:,take_array);
        [InvModels, Predicted] = ModelClean(InvModels, Predicted);
        str = sprintf('      Leave one out Regression = %s ',char(Parameters.SetNames(i)));
        disp(str)
        [weightings,StdsBeta]  = RegresAlgo(Parameters,InvModels,Predicted);
        %display(weightings)
        OrgWeights(take_array,i)= weightings;
        % Restore full dataset
        InvModels = model_values_store;
        Predicted = InvModels(:,i);
        InvModels =  InvModels(:,take_array);
        [model_values,~,new_weights] = Weightin_algo(InvModels,weightings,Parameters);
        testArray = [Predicted,model_values];
        [Outputs] = Accuracy_statistics(testArray,Parameters,Sizes);% run statistics
        OrgDev(i) = Outputs.mean_double_deviation;
        OrgRho(i) = Outputs.RHO;
        NewWeights(take_array,i) = new_weights;
        NewCV(take_array,i) = StdsBeta./new_weights;
        NewCV(NewWeights == 0) = 0;
        clear Outputs take_array weightings InvModels Predicted testArray model values
    end
    AvgDev = nanmean(OrgDev);
    AvgRho = nanmean(OrgRho);
    for i = 1:1:nrM
        take_array = RemvOne(nrM,i);
        Weights(i,1) = nanmean(NewWeights(i,take_array));
        CVErr(i,1)= nanmean(NewCV(i,take_array));
    end
    Weights = ShapeWeights(Weights,Parameters);
    InvModels = model_values_store; % Restore full dataset
    [model_values,~,new_weights] = Weightin_algo(InvModels,Weights,Parameters);
    Weighting.LeaveOneOut(:,1) = new_weights;
    Weighting.LeaveOneOut(:,2) = CVErr.*new_weights;
    Weighting.LeaveOneOut(:,3) = CVErr;
    Weighting.MaxentStart = [AvgDev; AvgRho];
elseif Ensemble == 9
    disp('      Iterating Cross validation optimal fit')
    % Start algo, make general function by using different names for in and
    % output
    Store.Dev(1) = Weighting.MaxentStart(1);
    Store.Rho(1) = Weighting.MaxentStart(2);
    Store.Least(1) = 0;
    Store.Weights(:,1) = Weighting.LeaveOneOut(:,1);
    Deltas =  Weighting.LeaveOneOut(:,3);
    nrM = size(model_values_store,2);
    for runs = 1:Parameters.RunsMaxent
        Weights = TakeBetaWeights(Parameters,nrM,1);
        %Baseline
        [NewDev, NewRho,NewLeast] = CrossValidate(model_values_store,nrM,Weights,Parameters,Sizes);
        Store.Dev(2) = NewDev;
        Store.Rho(2) = NewRho;
        Store.Least(2) = NewLeast;
        Store.Weights(:,2) = Weights;
        Cur.Dev = NewDev;
        Cur.Rho = NewRho;
        Cur.Least = NewLeast;
        Cur.Weights = Weights;
        % clear Weights NewDev NewRho
        %% Loops
        %Initiate
        its = 1;
        improv = 1;
        precision = 8;
        if Ensemble == 8
            DecCrit = Cur.Dev;
        elseif Ensemble == 9
            DecCrit = Cur.Rho;
        else
            DecCrit = Cur.Least;
        end
        
        while its < Parameters.max_its
            if Parameters.testRun == 1 ||  Parameters.testRun == 2
                clc
                disp(str)
                display(['Runs =' num2str(runs)]);
                display(['Criterion Value =' num2str(DecCrit,precision)])
                display(['Iteration = ',int2str(its)])
                display(['Improvement # = ',int2str(improv)])
            end
            if rem(its,5) == 0
                Cur.Weights = TakeBetaWeights(Parameters,nrM,its);
            end
            for i = 1:1:nrM
                Nwght = Cur.Weights;
                Nwght(i) = DeltaWght(Cur.Weights(i),Deltas(i));
                Nwght = ShapeWeights(Nwght,Parameters);
                % Cross Validate Models
                [NewDev, NewRho,NewLeast] = CrossValidate(model_values_store,nrM,Nwght,Parameters,Sizes);
                if NewRho == -9999
                    break
                end
                DecCrit = Cur.Least;
                Crit = NewLeast;
                [Cur, ~, improv, its] = DecisionTake(Crit,DecCrit,Store,improv,Cur,its,Nwght,NewDev,NewRho,NewLeast,Parameters);
            end
            its = its + 1;
        end
        CurWeights(:,runs) = Cur.Weights;
    end
    NewWeights = nanmedian(CurWeights,2);
    Weights = ShapeWeights(NewWeights,Parameters);
    Cur.Weights = Weights;
    % screen output of better deviance (maxent style, add when it works)
    % calcuate the actual model values
    InvModels = model_values_store; % Restore full dataset
    [model_values,~,new_weights] = Weightin_algo(InvModels,Cur.Weights,Parameters);
    Store.TimesImproved = improv;
    Weighting.MaxentLeastSquares(:,1) = new_weights;
    MaxentStore.LeastSquares = Store;
    clear Store Cur
elseif Ensemble == 10 % GridSize
    disp('      Gridsize differences among maps')
    NewWeights = reshape((1./(log10(Parameters.GridSizes))),length(Parameters.GridSizes),1);
    Weights = ShapeWeights(NewWeights,Parameters);
    InvModels = model_values_store; % Restore full dataset
    [model_values,~,new_weights] = Weightin_algo(InvModels,Weights,Parameters);
    Weighting.GridSize(:,1) = Weights;
    clear Weights NewWeights
elseif Ensemble == 11 % Uniqeness
    disp('      Model Uniqueness (upweighting unique models)')
    NewWeights = zeros(length(Parameters.groups),1);
    tst = unique(Parameters.groups);
    for i = 1: length(tst)
        tl = find(Parameters.groups == tst(i));
        NewWeights(tl) = 1./(length(tl)./length(Parameters.groups));
    end
    Weights = ShapeWeights(NewWeights,Parameters);
    InvModels = model_values_store; % Restore full dataset
    [model_values,~,new_weights] = Weightin_algo(InvModels,Weights,Parameters);
    Weighting.UniquenessUp(:,1) = new_weights;
    clear Weights NewWeights
elseif Ensemble == 12 % Uniqueness Flipped
    disp('      Model Uniqueness (downweighting unique models)')
    NewWeights = Weighting.UniquenessUp(:,1);
    [model_values, new_weights] = FlipAlgo(NewWeights,Parameters,model_values_store);
    Weighting.UniquenessDown(:,1) = new_weights;
    clear Weights NewWeights
    
    %% Half_informed Ensembles
elseif  Ensemble == 13
    % Deviance weighted stats
    %Weights scaling
    Deviances = Weighting.Deviance(:,1);
    Deviances_cor_max = max(Deviances);
    Deviances_cor_range = Deviances_cor_max  - min(Deviances);
    for k = 1:length(Deviances)
        Deviances_cor(k) = Deviances_cor_max - ((Deviances_cor_max - Deviances(k)).*(1./Deviances_cor_range)); %#ok<*SAGROW>
    end
    Deviances = Deviances_cor;
    clear k Deviances_*
    Deviances(Deviances < 0.25) = 0.25;
    Weighting.Deviance(:,2) = Deviances; % Store to compare
    % here we start predicting the other part of the data-set
    InvModels = model_values_S2; % note models_values_used = S2
    [model_values,~,new_weights] = Weightin_algo(InvModels,Deviances,Parameters);
    Weighting.Deviance(:,3) = new_weights;
    validation_set = validation_S2(:,Parameters.valiNaN);
elseif  Ensemble == 14
    % Rho weighted stats
    %Weights scaling
    Rhos_overview =  (Weighting.Rho(:,1)+1)./2;
    Rhos_cor_max = max(Rhos_overview);
    Rhos_cor_range = Rhos_cor_max  -min(Rhos_overview);
    for k = 1:length(Rhos_overview)
        Rhos_cor(k) = Rhos_cor_max - ((Rhos_cor_max - Rhos_overview(k)).*(1./Rhos_cor_range));
    end
    Rhos = Rhos_cor;
    clear k Rhos_*
    
    Rhos(Rhos < 0.25) = 0.25;
    Weighting.Rho(:,2) = Rhos; % Store to compare
    % here we start predicting the other part of the data-set
    InvModels = model_values_S2; % note models_values_used = S2
    [model_values,~,new_weights] = Weightin_algo(InvModels,Rhos,Parameters);
    Weighting.Rho(:,3) = new_weights;
    validation_set = validation_S2(:,Parameters.valiNaN);
elseif Ensemble == 15
    % Linear model Weight training
    InvModels = model_values_store;
    Predicted = validation_set;
    [InvModels, Predicted] = ModelClean(InvModels, Predicted);
    disp('      Regressing fit to Training set')
    [Trained_Weights,StdsBeta]  = RegresAlgo(Parameters,InvModels,Predicted);
    % transfer to Validation set
    InvModels = model_values_S2; % Restore full dataset from clean
    [model_values,~,new_weights] = Weightin_algo(InvModels,Trained_Weights,Parameters);
    Weighting.Trained_Weights(:,1) = new_weights;
    validation_set = validation_S2(:,Parameters.valiNaN);
    clear   Trained_Weights
elseif Ensemble == 16
    % Iterated model Weight training
    InvModels = model_values_store;
    Predicted = validation_set;
    disp('     Iterating maximum fit to Training set')
    % Start algo, make general function by using different names for in and
    % output
    Deltas =  Weighting.LeaveOneOut(:,3);
    nrM = size(model_values_store,2);
    Weights = TakeBetaWeights(Parameters,nrM,1);
    Cur.Weights = Weights;
    Cur.Dev = 0;
    Cur.Rho = 0;
    Cur.Least = -inf;
    Store.Weights(:,1) = Cur.Weights;
    Store.Dev(1) = Cur.Dev;
    Store.Rho(1) = Cur.Rho;
    Store.Least(1) = Cur.Least;
    %Initiate
    its = 1;
    improv = 1;
    precision = 4;
    while its < Parameters.max_itsTRained
        if Parameters.testRun == 1 || Parameters.testRun == 2
            clc
            disp(str)
            if exist('DecCrit','var')~= 0
                display(['Criterion Value =' num2str(DecCrit,precision)])
            end
            display(['Iteration = ',int2str(its)])
            display(['Improvement = ',int2str(improv)])
        end
        if rem(its,5) == 0
            [Cur.Weights] = TakeBetaWeights(Parameters,nrM,its);
        end
        for i = 1:1:nrM
            Nwght = Cur.Weights;
            Nwght(i) = DeltaWght(Cur.Weights(i),Deltas(i));
            Nwght = ShapeWeights(Nwght,Parameters);
            % Calcuate Outpust
            [model_values,~,~] = Weightin_algo(InvModels,Nwght,Parameters);
            testArray = [Predicted,model_values];
            [Outputs] = Accuracy_statistics(testArray,Parameters,Sizes);% run statistics
            if Outputs.RHO ~= -9999;
                NewDev = Outputs.mean_double_deviation;
                NewRho = Outputs.RHO;
                NewLeast = Outputs.LeastSquares;
                DecCrit = Cur.Least;
                Crit = NewLeast;
                [Cur, Store, improv, its] = DecisionTake(Crit,DecCrit,Store,improv,Cur,its,Nwght,NewDev,NewRho,NewLeast,Parameters);
                clear Nwght NewDev NewRh
            end
        end
        its = its + 1;
    end
    % transfer to Validation set
    InvModels = model_values_S2; % Restore full dataset from clean
    [model_values,~,new_weights] = Weightin_algo(InvModels,Cur.Weights,Parameters);
    validation_set = validation_S2(:,Parameters.valiNaN);
    Store.TimesImproved = improv;
    Weighting.Trained_IteratedLeast(:,1) = new_weights;
    MaxentStore.Least = Store;
    clear Store Cure
end
end
%%
function [model_values, new_weights] = FlipAlgo(NewWeights,Parameters,model_values_store)
NewWeights(find(NewWeights == 0)) = NaN; %#ok<FNDSB>
NewWeights = 1./NewWeights;
Weights = ShapeWeights(NewWeights,Parameters);
InvModels = model_values_store; % Restore full dataset
[model_values,~,new_weights] = Weightin_algo(InvModels,Weights,Parameters);
end

function [ValOut,StdOut,WdS] = Weightin_algo(InvM,WdS,Parameters)
for k = 1:1:length(InvM)
    WdS = ShapeWeights(WdS,Parameters);
    nomi = 0;
    denomi = 0;
    for j = 1:1:length(WdS)
        if isnan(InvM(k,j))~= 1
            nomi = nomi + (InvM(k,j).* WdS(j));
            denomi = denomi +  WdS(j);
        end
    end
    % denomi keeps needed since sometimes values drop out as NaN for indiviudal data points
    % so sum of weights < 1 so correction is needed
    ValOut(k,1) = nomi/denomi;
    clear nomi
    nomi = 0; %denomi is same as above
    for j = 1:1:length(WdS)
        if isnan(InvM(k,j))~= 1
            nomi = nomi + (((InvM(k,j)- ValOut(k))^2).*WdS(j));
        end
    end
    StdOut(k,1) =  nomi/denomi;
    clear nomi denomi
end
end

function  InV = ShapeWeights(InV,Parameters)
InV = reshape(InV,(length(InV)),1);
InV= max(InV,0);
InV = InV./sum(InV);
InV= min(max(InV,0),1);
InV = (round(InV.*Parameters.Precision(1)))./Parameters.Precision(1);
end

function  [wght] = TakeBetaWeights(Parameters,nrM,nrIt)
test = 0;
while test == 0
    if nrIt == 1 || rem(nrIt,5) == 0
        name = 'beta';
        BetarRnd= abs(randn);
    elseif  rem(nrIt,10) == 0
        name = 'gamma';
        BetarRnd= abs(2 + randn);
    elseif rem(nrIt,15) == 0
        name = 'gamma';
        BetarRnd= abs(0.5 + randn);
    end
    wght(:,1) = ShapeWeights((random(name,BetarRnd,BetarRnd,nrM,1)),Parameters);
    if isempty(find(isnan(wght),1))== 1 && isempty(find(isinf(wght),1)) == 1
        test = 1;
    end
end
end

function wght = DeltaWght(WghtIn,VarIn)
wghtmp = -1;
while wghtmp<0
    wghtmp = WghtIn +  (WghtIn.*(randn.* VarIn));
end
wght = wghtmp;
end

function [Inv,InvP] = ModelClean(Inv,InvP)
num = [];
for f = size(Inv,1):-1:1
    Inv(f,(isnan(Inv(f,:))==1)) = nanmean(Inv(f,:));
    testArray = find(isnan(Inv(f,:))~=1);
    if isempty(testArray==1)
        num = [num;f];
    end
end
if isempty(InvP)~=1
    num = unique([num;(find(isnan(InvP)==1))]);
    InvP(num) = [];
end
Inv(num,:) =[];
end

function take_array = RemvOne(nrM,i)
take_array = 1:nrM;
take_array(take_array==i) = [];
end

function  [NewDev, NewRho,NewLeast] = CrossValidate(model_in,nrM,wght,Parameters,Sizes)
% Cross Validate Models
for j = 1:nrM
    InvModels = model_in;
    Predicted = InvModels(:,j);
    wghts = wght;
    wghts = ShapeWeights(wghts,Parameters);
    [model_values,~,~] = Weightin_algo(InvModels,wghts,Parameters);
    testArray = [Predicted,model_values];
    [Outputs] = Accuracy_statistics(testArray,Parameters,Sizes);% run statistics
    if Outputs.RHO ~= -9999;
        NewDevs(j) = Outputs.mean_double_deviation;
        NewRhos(j) = Outputs.RHO;
        NewLeasts(j) = Outputs.LeastSquares;
    end
end
if exist('NewRhos','var') == 1
    NewDev = (round(nanmedian(NewDevs).*Parameters.Precision(2)))./Parameters.Precision(2);
    NewRho = (round(nanmedian(NewRhos).*Parameters.Precision(2)))./Parameters.Precision(2);
    NewLeast = (round(nanmedian(NewLeasts).*Parameters.Precision(2)))./Parameters.Precision(2);
else
    NewDev =-9999;
    NewRho = -9999;
    NewLeast = -9999;
end
end

function  [Cur, Store, improv, its] = DecisionTake(Crit,DecCrit,Store,improv,Cur,its,Nwght,NewDev,NewRho,NewLeast,~)
if Crit > DecCrit
    [Cur,Store,improv] = SetPara(Nwght,NewDev,NewRho,NewLeast,Cur,improv);
    % Reset all running parameters
    its = 1;
elseif Crit == DecCrit
    Dice = rand;
    if (Dice>0.5)
        [Cur,Store,improv] = SetPara(Nwght,NewDev,NewRho,NewLeast,Cur,improv);
        its = its; %#ok<ASGSL> %Note no change in iterations to avoid endless loops
    end
end
    function [Cur,Store,improv] = SetPara(Nwght,NewDev,NewRho,NewLeast,Cur,improv)
        Cur.Weights = Nwght;
        Cur.Dev =  NewDev;
        Cur.Rho =  NewRho;
        Cur.Least = NewLeast;
        Store.Dev(improv+2) = NewDev;
        Store.Rho(improv+2) = NewRho;
        Store.Least(improv+2) = NewLeast;
        Store.Weights(:,improv+2) =  Cur.Weights;
        improv = improv + 1;
    end

end

function [WghtOut,StdsBeta] = RegresAlgo(Parameters,InParY,InParX)
[modelFunw,series, prior,Options,Group] = CreateModelFun(InParY);
StdsBeta = ones(size(InParY,2),1).*Parameters.delta(1);
if Parameters.testRun ==1
    Options = statset('FunValCheck','off','Display','off','MaxIter',200,...
        'TolFun',1.0000e-4, 'TolX',1.0e-4, 'DerivStep', 6.0555e-06,'Robust',...
        'off', 'WgtFun', 'bisquare');
    [WghtOut,~,~,~,~]  = nlinfit(InParY,InParX,modelFunw,prior,'Options',Options);
else
    [WghtOut,~,stats] = nlmefit(InParY,InParX,Group,[],modelFunw,prior,'RefineBeta0','on','Options',Options,'FEParamsSelect',series); %#ok<*NASGU>
    StdsBeta = reshape(stats.sebeta,length(stats.sebeta),1);
end
WghtOut = ShapeWeights(WghtOut,Parameters);
end

function [modelFunw,series, prior,Options,Group] = CreateModelFun(InvModels)
series = [1:size(InvModels,2)];
prior = (ones((size(InvModels,2)),1)/(size(InvModels,2)))';
if (size(InvModels,2)) == 11
    modelFun = @(b,InvModels) (b(1).*InvModels(:,1)) + (b(2).*InvModels(:,2)) + (b(3).*InvModels(:,3)) +...
        (b(4).*InvModels(:,4)) + (b(5).*InvModels(:,5))+ (b(6).*InvModels(:,6)) + (b(7).*InvModels(:,7)) + ...
        (b(8).*InvModels(:,8)) + (b(9).*InvModels(:,9))+(b(10).*InvModels(:,10))+(b(11).*InvModels(:,11));
elseif (size(InvModels,2)) == 10
    modelFun = @(b,InvModels) (b(1).*InvModels(:,1)) + (b(2).*InvModels(:,2)) + (b(3).*InvModels(:,3)) +...
        (b(4).*InvModels(:,4)) + (b(5).*InvModels(:,5))+ (b(6).*InvModels(:,6)) + (b(7).*InvModels(:,7)) + ...
        (b(8).*InvModels(:,8)) + (b(9).*InvModels(:,9))+(b(10).*InvModels(:,10));
elseif (size(InvModels,2)) == 9
    modelFun = @(b,InvModels) (b(1).*InvModels(:,1)) + (b(2).*InvModels(:,2)) + (b(3).*InvModels(:,3)) +...
        (b(4).*InvModels(:,4)) + (b(5).*InvModels(:,5))+ (b(6).*InvModels(:,6)) + (b(7).*InvModels(:,7)) + ...
        (b(8).*InvModels(:,8)) + (b(9).*InvModels(:,9));
elseif (size(InvModels,2)) == 8
    modelFun = @(b,InvModels) (b(1).*InvModels(:,1)) + (b(2).*InvModels(:,2)) + (b(3).*InvModels(:,3)) +...
        (b(4).*InvModels(:,4)) + (b(5).*InvModels(:,5))+ (b(6).*InvModels(:,6)) + (b(7).*InvModels(:,7)) + ...
        (b(8).*InvModels(:,8));
elseif (size(InvModels,2)) == 7
    modelFun = @(b,InvModels) (b(1).*InvModels(:,1)) + (b(2).*InvModels(:,2)) + (b(3).*InvModels(:,3)) +...
        (b(4).*InvModels(:,4)) + (b(5).*InvModels(:,5))+ (b(6).*InvModels(:,6)) + (b(7).*InvModels(:,7));
elseif (size(InvModels,2)) == 6
    modelFun = @(b,InvModels) (b(1).*InvModels(:,1)) + (b(2).*InvModels(:,2)) + (b(3).*InvModels(:,3)) +...
        (b(4).*InvModels(:,4)) + (b(5).*InvModels(:,5))+ (b(6).*InvModels(:,6));
elseif (size(InvModels,2)) == 5
    modelFun = @(b,InvModels) (b(1).*InvModels(:,1)) + (b(2).*InvModels(:,2)) + (b(3).*InvModels(:,3)) +...
        (b(4).*InvModels(:,4)) + (b(5).*InvModels(:,5));
elseif (size(InvModels,2)) == 4
    modelFun = @(b,InvModels) (b(1).*InvModels(:,1)) + (b(2).*InvModels(:,2)) + (b(3).*InvModels(:,3)) +...
        (b(4).*InvModels(:,4)) ;
elseif (size(InvModels,2)) == 3
    modelFun = @(b,InvModels) (b(1).*InvModels(:,1)) + (b(2).*InvModels(:,2)) + (b(3).*InvModels(:,3));
end
modelFunw = @(b,InvModels)modelFun(b,InvModels);
Options = statset('FunValCheck','on','Display','off','MaxIter',100,...
    'TolFun',1.0000e-2, 'TolX',1.0e-2, 'Jacobian','off', 'DerivStep', 6.0555e-03, 'OutputFcn',' ');
Group = grp2idx(InvModels(:,1));
end