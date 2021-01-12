function [ResultsCombi,ResultsWeights, VariationWeights] = JoinFunc(Parameters)
%[Parameters, ~,~] = DefintionSet_Carbon(validation_set,Parameters);
%load('necesary_para')
 cd('Output_Dir')
 clear ResultsCombi IB ResultsWeights
 fac = 1000;
 fac2 = 10000;
 ResultsCombi.Models = dataset({'Dummy'},'VarNames',{'Data_set'});
 nrM = length(Parameters.SetNames);
for run = 1:Parameters.runMax
    Output_file = [Parameters.output_file,'_',int2str(run),'.mat'];
    if Parameters.testRun == 2
        Output_file = [Parameters.output_file,'_Full'];
    end
    load(Output_file);
    
    leni = length(Results.Models.Data_set);
  for i = 1:leni
     IB.Dev(run,i) = Results.Models.Inversed_deviance(i); %#ok<*SAGROW>
     IB.Rho(run,i) =  Results.Models.RHO(i);
     IB.PVal(run,i) = Results.Models.PVal(i);
     IB.Name(i,1) = Results.Models.Data_set(i);
     IB.DataPoints(run,i) = Results.Models.Datapoints(i);
     IB.Improvs_Dev(i,1:leni) = (Results.Models.Inversed_deviance./Results.Models.Inversed_deviance(i));
     IB.Improvs_Rho(i,1:leni) = ((Results.Models.RHO+1)./2) ./((Results.Models.RHO(i)+1)./2) ;
     %Against summed Models
     IB.Mod_Improv.Dev(run,i) = 1./(nanmean(Results.Models.Inversed_deviance(i)./(Results.Models.Inversed_deviance(1:nrM))));
     IB.Mod_Improv.Rho(run,i) =  1./(nanmean(((Results.Models.RHO(i)+1)/2)./((Results.Models.RHO(1:nrM)+1)./2)));
  end    
  IB.Dev(run,leni+1) = (sum((IB.Dev(run,1:nrM).*IB.DataPoints(run,(1:nrM)))))./sum(IB.DataPoints(run,(1:nrM)));
  IB.Rho(run,leni+1) =  (sum((IB.Rho(run,1:nrM).*IB.DataPoints(run,(1:nrM)))))./sum(IB.DataPoints(run,(1:nrM)));
  IB.PVal(run,leni+1) = (sum((IB.PVal(run,1:nrM).*IB.DataPoints(run,(1:nrM)))))./sum(IB.DataPoints(run,(1:nrM)));
  IB.Name(leni+1,1) = {'Models_as_group'};
  IB.DataPoints(run,leni+1) = (sum((IB.DataPoints(run,1:nrM).*IB.DataPoints(run,(1:nrM)))))./sum(IB.DataPoints(run,(1:nrM)));
  
  IB.Improvement_Dev(run,1:leni,1:leni) = IB.Improvs_Dev;
  IB.Improvement_Dev(run,1:leni,leni+1) = IB.Mod_Improv.Dev(run,:);
  IB.Improvement_Dev(run,leni+1,1:leni) = (1./(IB.Mod_Improv.Dev(run,:)))';
  IB.Improvement_Dev(run,leni+1,leni+1) = 1;
    
  IB.Improvement_Rho(run,1:leni,1:leni) = IB.Improvs_Rho;
  IB.Improvement_Rho(run,1:leni,leni+1) = IB.Mod_Improv.Rho(run,:);
  IB.Improvement_Rho(run,leni+1,1:leni) = (1./(IB.Mod_Improv.Rho(run,:)))';
  IB.Improvement_Rho(run,leni+1,leni+1) = 1;
  
  IB.SummedWeights.PCA(:,run) = Weighting.PCA(:,1);
  IB.SummedWeights.RegresstoMedian(:,run) = Weighting.RegresstoMedian (:,1);
  IB.SummedWeights.MaxentLeastSquares(:,run) = Weighting.MaxentLeastSquares(:,1);
  IB.SummedWeights.LeaveOneOut(:,run) = Weighting.LeaveOneOut(:,1);
    IB.SummedWeights.CorCoef(:,run) = Weighting.CorCoef(:,1);
  IB.SummedWeights.GridSize(:,run) = Weighting.GridSize(:,1);
  IB.SummedWeights.UniquenessUp(:,run) = Weighting.UniquenessUp(:,1);
  IB.SummedWeights.UniquenessDown(:,run) = Weighting.UniquenessDown(:,1);
  
  IB.VarWeights.PCA(run) = nanstd(Weighting.PCA(:,1));
  IB.VarWeights.RegresstoMedian(run) = nanstd(Weighting.RegresstoMedian (:,1));
   IB.VarWeights.MaxentLeastSquares(run) = nanstd(Weighting.MaxentLeastSquares(:,1));
  IB.VarWeights.LeaveOneOut(run) = nanstd(Weighting.LeaveOneOut(:,1));
   IB.VarWeights.CorCoef(run) = nanstd(Weighting.CorCoef(:,1));
  IB.VarWeights.GridSize(run) = nanstd(Weighting.GridSize(:,1));
  IB.VarWeights.UniquenessUp(run) = nanstd(Weighting.UniquenessUp(:,1));
  IB.VarWeights.UniquenessDown(run) = nanstd(Weighting.UniquenessDown(:,1));
  
if Parameters.testRun ~= 2
  IB.SummedWeights.HI_Deviance(:,run) = Weighting.Deviance(:,3); %Half
  IB.SummedWeights.HI_Rho(:,run) = Weighting.Rho(:,3); % Half
  IB.SummedWeights.BaggedRegres(:,run) =Weighting.Trained_Weights(:,1);% Half
  IB.SummedWeights.Trained_IteratedLeast(:,run) = Weighting.Trained_IteratedLeast(:,1); % half
  IB.VarWeights.HI_Deviance(run) = nanstd(Weighting.Deviance(:,3)); %Half
  IB.VarWeights.HI_Rho(run) = nanstd(Weighting.Rho(:,3)); % Half
  IB.VarWeights.BaggedRegres(run) =nanstd(Weighting.Trained_Weights(:,1));% Half
  IB.VarWeights.Trained_IteratedLeast(run) = nanstd(Weighting.Trained_IteratedLeast(:,1)); % half
  end
end
cd ..

ResultsCombi.Models.Data_set = IB.Name;
ResultsCombi.Models.DataPoints(:,1) = round(nanmean(IB.DataPoints,1));
ResultsCombi.Models.RhoMean(:,1) = (round(nanmean(IB.Rho,1).*fac))./fac;
ResultsCombi.Models.RhoStd(:,1) = (round(nanstd(IB.Rho,0,1).*fac2))./fac2;
ResultsCombi.Models.RhoPval(:,1) = (round(nanmean(IB.PVal,1).*fac2))./fac2;
ResultsCombi.Models.DevianceMean(:,1) = (round(nanmean(IB.Dev,1).*fac))./fac;
ResultsCombi.Models.DevianceStd(:,1) = (round(nanstd(IB.Dev,0,1).*fac2))./fac2;
%weighted means for Best models
ResultsCombi.Models.DevianceMean((nrM+1),1) = (round((sum((IB.Dev(:,nrM+1).*IB.DataPoints(:,(nrM+1)))))./sum(IB.DataPoints(:,(nrM+1))).*fac))./fac;
ResultsCombi.Models.RhoMean((nrM+2),1) = (round((sum((IB.Rho(:,nrM+2).*IB.DataPoints(:,(nrM+2)))))./sum(IB.DataPoints(:,(nrM+2))).*fac))./fac;

% Improvements
for iter = 1:1:(Parameters.ImprovementNrTake*2)
    Leni = size(IB.Improvement_Rho,1);
    PerNr = (randperm(Leni));
    TkPerNr = PerNr(1:Parameters.ImprovementNrTake);
    for j = 1:1:size(IB.Improvement_Rho,2)
        for k = 1:1:size(IB.Improvement_Rho,3)
            [~,Signi_Rho(iter,j,k),~,~]=ttest(IB.Improvement_Rho(TkPerNr,j,k),1); %#ok<AGROW>
             [~,Signi_Dev(iter,j,k),~,~]=ttest(IB.Improvement_Dev(TkPerNr,j,k),1); %#ok<AGROW>
        end
    end
end

ResultsCombis.Improvement.Rho.Significance =(round(reshape(nanmean(Signi_Rho,1),leni+1,leni+1).*100000))./100000;
ResultsCombis.Improvement.Deviance.Significance =(round(reshape(nanmean(Signi_Dev,1),leni+1,leni+1).*100000))./100000;
clear Signi* Leni PerNr TkPerNr iter

ResultsCombis.Improvement.Rho.Mean = (round(reshape(nanmean(IB.Improvement_Rho,1),leni+1,leni+1).*fac))./fac;
ResultsCombis.Improvement.Rho.Std = (round(reshape(nanstd(IB.Improvement_Rho,0,1),leni+1,leni+1).*fac2))./fac2;
ResultsCombis.Improvement.Deviance.Mean = (round(reshape(nanmean(IB.Improvement_Dev,1),leni+1,leni+1).*fac))./fac;
ResultsCombis.Improvement.Deviance.Std = (round(reshape(nanstd(IB.Improvement_Dev,0,1),leni+1,leni+1).*fac2))./fac2;

Names = ResultsCombi.Models.Data_set';
Names(leni+1) = {'Models_as_group'};
ResultsCombi.Improvement.Rho.Significance = dataset(ResultsCombis.Improvement.Rho.Significance(:,1),'ObsNames', Names','VarNames',{char(Names(1))});
ResultsCombi.Improvement.Deviance.Significance = dataset(ResultsCombis.Improvement.Deviance.Significance(:,1),'ObsNames', Names','VarNames',{char(Names(1))});
ResultsCombi.Improvement.Rho.Mean = dataset(ResultsCombis.Improvement.Rho.Mean(:,1),'ObsNames', Names','VarNames',{char(Names(1))});
ResultsCombi.Improvement.Rho.Std = dataset(ResultsCombis.Improvement.Rho.Std(:,1),'ObsNames', Names','VarNames',{char(Names(1))});
ResultsCombi.Improvement.Deviance.Mean = dataset(ResultsCombis.Improvement.Deviance.Mean(:,1),'ObsNames', Names','VarNames',{char(Names(1))});
ResultsCombi.Improvement.Deviance.Std = dataset(ResultsCombis.Improvement.Deviance.Std(:,1),'ObsNames', Names','VarNames',{char(Names(1))});
for t = 1:1:(leni+1)
    ResultsCombi.Improvement.Rho.Significance.(genvarname(char(Names(t))))= ResultsCombis.Improvement.Rho.Significance(:,t);
    ResultsCombi.Improvement.Deviance.Significance.(genvarname(char(Names(t))))= ResultsCombis.Improvement.Deviance.Significance(:,t); 
    ResultsCombi.Improvement.Rho.Mean.(genvarname(char(Names(t)))) =  ResultsCombis.Improvement.Rho.Mean(:,t);
    ResultsCombi.Improvement.Rho.Std.(genvarname(char(Names(t)))) = ResultsCombis.Improvement.Rho.Std(:,t);
    ResultsCombi.Improvement.Deviance.Mean.(genvarname(char(Names(t)))) = ResultsCombis.Improvement.Deviance.Mean(:,t);
    ResultsCombi.Improvement.Deviance.Std.(genvarname(char(Names(t)))) = ResultsCombis.Improvement.Deviance.Std(:,t);
end

%% weights
ResultsWeights = dataset((nanmean(IB.SummedWeights.PCA,2)),'Varnames',{'PCA'});
ResultsWeights.RegresstoMedian(:,1) = nanmean(IB.SummedWeights.RegresstoMedian,2);
ResultsWeights.LeaveOneOut(:,1) = nanmean(IB.SummedWeights.LeaveOneOut,2);
ResultsWeights.MaxentLeastSquares(:,1) =  nanmean(IB.SummedWeights.MaxentLeastSquares,2);
ResultsWeights.CorCoef(:,1) = nanmean(IB.SummedWeights.CorCoef,2);
ResultsWeights.GridSize(:,1) = nanmean(IB.SummedWeights.GridSize,2);
ResultsWeights.Upweighted_Uniqueness(:,1) = nanmean(IB.SummedWeights.UniquenessUp,2);
ResultsWeights.Downweighted_Uniqueness(:,1) = nanmean(IB.SummedWeights.UniquenessDown,2);

 if Parameters.testRun ~= 2
ResultsWeights.HI_Deviance(:,1) = nanmean(IB.SummedWeights.HI_Deviance,2);
ResultsWeights.HI_Rho(:,1) = nanmean(IB.SummedWeights.HI_Rho,2);
ResultsWeights.BaggedRegres(:,1) = nanmean(IB.SummedWeights.BaggedRegres,2);
ResultsWeights.Trained_IteratedLeast(:,1) =  nanmean(IB.SummedWeights.Trained_IteratedLeast,2);
ResultsWeights.HI_Deviance_Std(:,1) = nanstd(IB.SummedWeights.HI_Deviance,0,2);
ResultsWeights.HI_Rho_Std(:,1) = nanstd(IB.SummedWeights.HI_Rho,0,2);
ResultsWeights.BaggedRegres_Std(:,1) = nanstd(IB.SummedWeights.BaggedRegres,0,2);
ResultsWeights.Trained_IteratedLeast_Std(:,1) =  nanstd(IB.SummedWeights.Trained_IteratedLeast,0,2);
 end
 
ResultsWeights.PCA_Std(:,1) = nanstd(IB.SummedWeights.PCA,0,2);
ResultsWeights.RegresstoMedian_Std(:,1) = nanstd(IB.SummedWeights.RegresstoMedian,0,2);
ResultsWeights.LeaveOneOut_Std(:,1) = nanstd(IB.SummedWeights.LeaveOneOut,0,2);
ResultsWeights.MaxentLeastSquares_Std(:,1) =  nanstd(IB.SummedWeights.MaxentLeastSquares,0,2);
ResultsWeights.CorCoef_Std(:,1) = nanstd(IB.SummedWeights.CorCoef,0,2);
ResultsWeights.GridSize_Std(:,1) = nanstd(IB.SummedWeights.GridSize,0,2);
ResultsWeights.Upweighted_Uniqueness_Std(:,1) = nanstd(IB.SummedWeights.UniquenessUp,0,2);
ResultsWeights.Downweighted_Uniqueness_Std(:,1) = nanstd(IB.SummedWeights.UniquenessDown,0,2);

 if Parameters.testRun ~= 2
VariationWeights.HI_Deviance(:,1) = nanmean(IB.VarWeights.HI_Deviance);
VariationWeights.HI_Rho(:,1) = nanmean(IB.VarWeights.HI_Rho);
VariationWeights.BaggedRegres(:,1) = nanmean(IB.VarWeights.BaggedRegres);
VariationWeights.Trained_IteratedLeast(:,1) =  nanmean(IB.VarWeights.Trained_IteratedLeast);
 end
 
VariationWeights = dataset((nanmean(IB.VarWeights.PCA)),'Varnames',{'PCA'});
VariationWeights.RegresstoMedian(:,1) = nanmean(IB.VarWeights.RegresstoMedian);
VariationWeights.LeaveOneOut(:,1) = nanmean(IB.VarWeights.LeaveOneOut);
VariationWeights.MaxentLeastSquares(:,1) =  nanmean(IB.VarWeights.MaxentLeastSquares);
VariationWeights.CorCoef(:,1) = nanmean(IB.VarWeights.CorCoef);
VariationWeights.GridSize(:,1) = nanmean(IB.VarWeights.GridSize);
VariationWeights.Upweighted_Uniqueness(:,1) = nanmean(IB.VarWeights.UniquenessUp);
VariationWeights.Downweighted_Uniqueness(:,1) = nanmean(IB.VarWeights.UniquenessDown); 



end