function [pat_table] = getBestRuns(pat_table, compCriteria,minOrMax)
%getBestRuns
%Returns a table, containing only the best runs of a model


%params, IN:
%   par_table
%           -> parameter table for the referred patient
%   compCriteria, minOrMax
%           -> criteria that shall be used for comparison, and if the best
%           model is the min or max regarding that criteria. The patients
%           table has to have a column named as the compCriteria

%OUT:
%   BestRuns
%           -> Table structure as input table, but from duplicate models
%           with different runs, only the best run is returned
    names=pat_table.Properties.VariableNames;
    if ~isempty(find(names=="model_combined", 1))
        model_names=names(~cellfun(@isempty,regexp(names,"model_combined")));
    else
        model_names=names(~cellfun(@isempty,regexp(names,"model")));
    end
    
    for models=unique(pat_table(:,model_names)).Variables'
        index=find(sum(string(pat_table(:,model_names).Variables)==string(models'),2)==length(models));
%    for model=unique(pat_table.model)'
%        index= find(pat_table.model==model);
        bestRun=getBestMeasure(pat_table(index,:), compCriteria,minOrMax);
        indexNotToDelete=find(abs(pat_table(index,compCriteria).Variables-bestRun)==0,1);
        index(indexNotToDelete)=[];
        pat_table(index,:)=[];    
    end
end