function [Results,Points,Weighting] = Make_ModelResults(Outputs,Parameters,Results,Points,Weighting,data_set)
%% All calculations
    Weighting.Deviance(data_set,1) = Outputs.mean_double_deviation;% for use in weighting
    Weighting.Rho(data_set,1) = Outputs.RHO; % for use in weighting
    Results.Models.Data_set(data_set,1) = (Parameters.SetNames(data_set));
    Results.Models.Datapoints(data_set,1) = Outputs.datapoints;
    Results.Models.RHO(data_set,1) = Outputs.RHO;
    Results.Models.PVal(data_set,1) = Outputs.PVAL;
    Results.Models.Inversed_deviance(data_set,1) = Outputs.mean_double_deviation;
    % Collect Individual data points
    if data_set == 1
        clear Points
        Points.Validation = dataset(Parameters.Names,'Varnames',char('Datapoint_name'));
        Points.Models = dataset(Parameters.Names,'Varnames',char('Datapoint_name'));
        Points.Deviation = dataset(Parameters.Names,'Varnames',char('Datapoint_name'));
    end
    Points.Validation.(genvarname(char(Parameters.SetNames(data_set))))  = Outputs.xes;
    Points.Models.(genvarname(char(Parameters.SetNames(data_set)))) = Outputs.yes;
    Points.Deviation.(genvarname(char(Parameters.SetNames(data_set))))= Outputs.deviation_point;
end