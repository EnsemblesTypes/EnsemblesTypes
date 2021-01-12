function   [Parameters, Models,Comparator] = DefintionSet(validation_set,Parameters)

%% Parameters various needed for testrun == 0 = bootstrap run
Parameters.CutOffType = min(find((strncmp(Parameters.Ensemble_Names,'Half',4)))); %#ok<*MXFND>
Parameters.Nr_ensembles = length(Parameters.Ensemble_Names);
Parameters.max_its = 250;
Parameters.max_itsTRained = 500;
Parameters.delta = [0.1;0.25];
%Parameters.testRun = 0; % To be added with function call
%Parameters.runMax = 250; % To be added with function call
Parameters.TProp = 0.5; % = 50% Training set
Parameters.Precision = 0.000001; % rounding factors
Parameters.NrCat = 4;
Parameters.ImprovementNrTake = 50;
Parameters.TimeOut = 120;
Parameters.RunsMaxent =1;

Parameters.Precision(2) = Parameters.Precision;%.*10;
Parameters.Precision = (1./Parameters.Precision);
Parameters.make_log = 0;
Parameters.ensemble = 0;
 
if Parameters.testRun == 1 % Just testing in screen mode
    Parameters.runMax = 5;
    Parameters.TimeOut = 5;
    Parameters.max_its =10;
    Parameters.ImprovementNrTake = 5;
    Parameters.RunsMaxent =2;
elseif Parameters.testRun == 2 %Full set (SI-3)
     Parameters.runMax = 1;
     Parameters.TProp = 1;
     Parameters.ImprovementNrTake = 1;
     Parameters.Nr_ensembles = length(Parameters.Ensemble_Names);
     Parameters.max_its = 100;
     Parameters.max_itsTRained = 500;
     Parameters.RunsMaxent =2;
else 
    % Parameters.testRun == 0; Standard run with parameters above
end
%% load data
if validation_set ==1
     Parameters.ServiceName = 'Carbon_tons';
    %GB set
    cd('TheData')
    load(Parameters.ServiceName)
    cd ..
    Parameters.make_log = 0;
    Parameters.output_file = 'Results_Carbon_GB';
    Parameters.SetNames = {'Aries';'Invest_L2015';'transfer_L2015';'Henrys';'Baredo';'Kindermann';'Density';'NEInventory';'LPJveg';'Luci'};
        Models.Service = [Carbon.GB.(genvarname(char(Parameters.SetNames(1)))), Carbon.GB.(genvarname(char(Parameters.SetNames(2)))),...
        Carbon.GB.(genvarname(char(Parameters.SetNames(3)))),Carbon.GB.(genvarname(char(Parameters.SetNames(4)))),...
        Carbon.GB.(genvarname(char(Parameters.SetNames(5)))),Carbon.GB.(genvarname(char(Parameters.SetNames(6)))),...
        Carbon.GB.(genvarname(char(Parameters.SetNames(7)))),Carbon.GB.(genvarname(char(Parameters.SetNames(8)))),...
        Carbon.GB.(genvarname(char(Parameters.SetNames(9)))),Carbon.GB.(genvarname(char(Parameters.SetNames(10))))];
    Models.Areas = [Sizes.GB.(genvarname(char(Parameters.SetNames(1)))), Sizes.GB.(genvarname(char(Parameters.SetNames(2)))),...
        Sizes.GB.(genvarname(char(Parameters.SetNames(3)))),Sizes.GB.(genvarname(char(Parameters.SetNames(4)))),...
        Sizes.GB.(genvarname(char(Parameters.SetNames(5)))),Sizes.GB.(genvarname(char(Parameters.SetNames(6)))),...
        Sizes.GB.(genvarname(char(Parameters.SetNames(7)))),Sizes.GB.(genvarname(char(Parameters.SetNames(8)))),...
        Sizes.GB.(genvarname(char(Parameters.SetNames(9)))),Sizes.GB.(genvarname(char(Parameters.SetNames(10))))];
    Parameters.GridSizes = [100,25,25,1000,1000,1000,20,20,10000,25];
    Parameters.data_set_max =  size(Models.Service,2);
    %Models.Areas = ones(length(Carbon.GB.Invest_L2007),Parameters.data_set_max);
    Parameters.names = Carbon.GB.Forest;
    Comparator.Service = Carbon.GB.ForestResearch;
    Parameters.Sizes = Sizes.GB.ForestResearchArea;
    Parameters.UniqueMatrix =  Carbon.UniquenessMatrix;
elseif validation_set ==2
    Parameters.ServiceName = 'Water_flows';
    %All set
     cd('TheData')
    load(Parameters.ServiceName)
    cd ..
    Parameters.make_log = 1;
    Parameters.output_file  = 'Results_water_All';
    Parameters.SetNames= {'Invest_L2015_GEAR';'WaterWorldBasin';'Growth_days';'transfer_L2015';'G2G';....
        'Aquaduct';'Decipher';'LPJ';'Luci'};
    Models.Service = [Flows.All.(genvarname(char(Parameters.SetNames(1)))), Flows.All.(genvarname(char(Parameters.SetNames(2)))),...
        Flows.All.(genvarname(char(Parameters.SetNames(3)))),((Flows.All.(genvarname(char(Parameters.SetNames(4)))))./Sizes.All.Population),...
        Flows.All.(genvarname(char(Parameters.SetNames(5)))),Flows.All.(genvarname(char(Parameters.SetNames(6)))),...
        Flows.All.(genvarname(char(Parameters.SetNames(7)))),Flows.All.(genvarname(char(Parameters.SetNames(8)))),...
        Flows.All.(genvarname(char(Parameters.SetNames(9))))];
    Models.Areas = [Sizes.All.(genvarname(char(Parameters.SetNames(1)))), Sizes.All.(genvarname(char(Parameters.SetNames(2)))),...
        Sizes.All.(genvarname(char(Parameters.SetNames(3)))),Sizes.All.(genvarname(char(Parameters.SetNames(4)))),...
        Sizes.All.(genvarname(char(Parameters.SetNames(5)))),Sizes.All.(genvarname(char(Parameters.SetNames(6)))),...
        Sizes.All.(genvarname(char(Parameters.SetNames(7)))),Sizes.All.(genvarname(char(Parameters.SetNames(8)))),...
        Sizes.All.(genvarname(char(Parameters.SetNames(9))))];
    Parameters.GridSizes = [25,1000,1000,25,1000,10000,25,10000,25];
    Parameters.data_set_max =  size(Models.Service,2);
    %Models.Areas = ones(length(Flows.All.Invest_L2015),data_set_max);
    Parameters.names = Flows.All.Catchment;
    Comparator.Service = Flows.All.TwentyYearFlow;
    Parameters.Sizes = Sizes.All.NFRASize;
    Parameters.UniqueMatrix =  Flows.UniquenessMatrix;
elseif validation_set == 3
     Parameters.ServiceName = 'Carbon_tons';
    %All set
     cd('TheData')
    load(Parameters.ServiceName)
    cd ..
    Parameters.make_log = 0;
    Parameters.output_file  = 'Results_Carbon_AfricaPlots';
    Parameters.SetNames= {'Invest';'CostingNature';'LPJ';'Benefit'};
    Models.Service = [Carbon.Africa.(genvarname(char(Parameters.SetNames(1)))), Carbon.Africa.(genvarname(char(Parameters.SetNames(2)))),...
        Carbon.Africa.(genvarname(char(Parameters.SetNames(3)))),Carbon.Africa.(genvarname(char(Parameters.SetNames(4))))];
    Models.Areas = [Carbon.Africa.Size, Carbon.Africa.Size,Carbon.Africa.Size,Carbon.Africa.Size];
     Parameters.GridSizes = [1000,1000,1000,1000];
    Parameters.data_set_max =  size(Models.Service,2);
    %Models.Areas = ones(length(Carbon.All.Invest_L2015),data_set_max);
    Parameters.names = Carbon.Africa.PlotID;
    Comparator.Service = Carbon.Africa.Validation;
    Parameters.Sizes = Carbon.Africa.Size;
    Parameters.NrCat = 2;
    Parameters.UniqueMatrix =  Carbon.AfricaUniquenessMatrix;
elseif validation_set == 4
     Parameters.ServiceName = 'Water_flows';
    %All set
     cd('TheData')
    load(Parameters.ServiceName)
    cd ..
    Parameters.make_log = 1;
    Parameters.output_file  = 'Results_water_GRDC';
    Parameters.SetNames= {'Invest';'WaterWorld';'CostingNature';'LPJ';'Benefit';'GrowthDays'};
    Models.Service = [Flows.Africa.(genvarname(char(Parameters.SetNames(1)))), Flows.Africa.(genvarname(char(Parameters.SetNames(2)))),...
        Flows.Africa.(genvarname(char(Parameters.SetNames(3)))),Flows.Africa.(genvarname(char(Parameters.SetNames(4)))),...
        Flows.Africa.(genvarname(char(Parameters.SetNames(5)))),Flows.Africa.(genvarname(char(Parameters.SetNames(6))))];
    Models.Areas = [Flows.Africa.Sizes, Flows.Africa.Sizes,Flows.Africa.Sizes,...
                       Flows.Africa.Sizes,Flows.Africa.Sizes,Flows.Africa.Sizes];
    Parameters.GridSizes = [1000,1000,1000,1000,1000,1000];
    Parameters.data_set_max =  size(Models.Service,2);
    %Models.Areas = ones(length(Flows.All.Invest_L2015),data_set_max);
    Parameters.names = Flows.Africa.CathmentID;
    Comparator.Service = Flows.Africa.Validator;
    Parameters.Sizes = Flows.Africa.Sizes;
    Parameters.NrCat = 2;
    Parameters.UniqueMatrix =  Flows.AfricaUniquenessMatrix;
end
end