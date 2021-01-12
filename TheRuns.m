function TheRuns(Parameters,Models,Comparator,run)
display(Parameters.output_file) % will only show when in test mode
display(['Running = ',int2str(run),'th run']) % will only show when in test mode
Leni = length(Comparator.Service);
PerNr = randperm(Leni);
% jack knife sample is equal for all models!
Parameters.Training = sort(PerNr(1:ceil(Parameters.TProp*Leni)));
if  Parameters.testRun < 2
Parameters.Validator = sort(PerNr((ceil(Parameters.TProp*Leni)+1):Leni));
else
Parameters.Validator = Parameters.Training;    
end
clear PerNr Leni
for data_set = 1:1:Parameters.data_set_max
    %% Make all values per hectare
    str = sprintf('Running Model= %s ',char(Parameters.SetNames(data_set)));
    disp(str)
    TmpX = (Comparator.Service./Parameters.Sizes);
    TmpY = (Models.Service(:,data_set)./Models.Areas(:,data_set));  
    if data_set == 1 % Initiate outputs
        Points = [];
        Weighting= [];
        Results.Models = dataset({'Dummy'},'Varnames',char({'Data_set'}));
    end
    % Jack knifed SET 1
    Training.Vali(:,1) =  TmpX(Parameters.Training);
    Training.Model(:,1) = TmpY(Parameters.Training);
    Parameters.Names(:,1) = Parameters.names(Parameters.Training);
    Sizes = Parameters.Sizes(Parameters.Training);
    clear testArray Outputs % make sure they are not lingering around
    testArray = [Training.Vali,Training.Model]; % X, Y
    [Outputs] = Accuracy_statistics(testArray,Parameters,Sizes);% run statistics
    [Results,Points,Weighting] = Make_ModelResults(Outputs,Parameters,Results,Points,Weighting,data_set);
    clear Training testArray Outputs Sizes
    if Parameters.testRun == 3 && data_set == Parameters.data_set_max
        disp(Results.Models)
       % error('myApp:argChk', 'Input arguments error');
    end
    
    % Jack knifed SET 2; what is needed is normalised model points and acc. validators
    Set2.Vali(:,1) =  TmpX(Parameters.Validator);
    Set2.Model(:,1) = TmpY(Parameters.Validator);
    Set2.Names(:,1)= Parameters.names(Parameters.Validator); 
    Sizes = Parameters.Sizes(Parameters.Validator);
    clear Tmp*
    if data_set == 1
        S2Points.Models = dataset(Set2.Names,'Varnames',char('Datapoint_name'));
    end
    clear testArray Outputs % make sure they are not lingering around
    testArray = [Set2.Vali,Set2.Model]; % X, Y
    [Outputs] = Accuracy_statistics(testArray,Parameters,Sizes);% run statistics
    S2Points.Models.(genvarname(char(Parameters.SetNames(data_set)))) = Outputs.yes;
    S2Points.Validation.(genvarname(char(Parameters.SetNames(data_set)))) = Outputs.xes;
    clear Set2 testArray Outputs
    clear Sizes
    %%  Ensembles run on results from individual model runs
    Sizes =Parameters.Sizes; %Reset Sizes
  
    if data_set ==  Parameters.data_set_max
        display ('  ') % will only show when in test mode
        display ('Ensemble calculations') % will only show when in test mode
        [Results, Points,Weighting,~] = Make_Ensembles(Results,Points,Weighting,Parameters,S2Points,Sizes); 
        % Replace Sensitivity by meta-Ensemble TO MAKE
        %Results = Sensitivity(Results,Points,Models,Parameters,Parameters.data_set_max);
        Points.Comparator = Parameters.Names;
        Points.Sizes =  Models.Areas;
        cd('Output_Dir')
        Output_file = [Parameters.output_file,'_',int2str(run)];
        if Parameters.testRun == 2
            Output_file = [Parameters.output_file,'_Full'];
        end
        save(Output_file,'Results','Points','Weighting','Parameters');
        cd ..
    end
end
display('   ')
display('   ')
end