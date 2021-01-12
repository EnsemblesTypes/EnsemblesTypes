function Parameters = ClusterTest(Parameters)
reply = 'N';
while reply ~= 'Y'
      ShowTextReply(Parameters)
    test = Parameters.UniqueMatrix(Parameters.SetNames,:);
    Y =squareform(pdist(double(test),'spearman'));
    Y(isnan(Y)==1) = 0;
    Z = linkage(Y);
    Parameters.groups = clusterdata(double(test),'distance','spearman','maxclust',Parameters.NrCat);
    if Parameters.testRun ~= 1
        reply = showFigure(Parameters,Z);
        str = sprintf('      The reply was: %s ',reply);
        disp(str)
        if strcmp('man',reply) == 1 || strcmp('Man',reply) ==1
            [Parameters,reply] = ManGroups(Parameters,Z);
        elseif strcmp('Y',reply) ~= 1 
            reply = 'N';
            Parameters.NrCat = str2double(input('Enter new number of categories: ', 's'));
        end
    else
          reply = 'Y';
    end 
end
close
end

function reply = showFigure(Parameters,Z)
close
for i = 1:length(Parameters.SetNames)
    Text =[int2str(Parameters.groups(i)),'_',char(Parameters.SetNames(i))];
    LabelsText(i) = {Text};     %#ok<*AGROW>
    clear Text
end
dendrogram(Z,'labels',LabelsText)
reply = input('      Are those the right category assignments Y/N/Man: ', 's');
if strcmp('y',reply) == 1
    reply = 'Y';
end
end

function [Parameters,reply] = ManGroups(Parameters,Z)
 Parameters.groups = (input('Enter array of categories [.;.;.](check order): '));
 Parameters.NrCat = length(unique(Parameters.groups));
 reply = showFigure(Parameters,Z);
 end

function ShowTextReply(Parameters)
str = sprintf('      Uniqueness will be done with %d categories ',Parameters.NrCat);
disp(str)
end
