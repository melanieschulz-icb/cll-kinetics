function [bestMeasure] = getBestMeasure(par_table, compCriteria,minOrMax)
%getBestMeasure
%Returns best measure


%params, IN:
%   par_table
%           -> parameter table for the referred patient
%   compCriteria, minOrMax
%           -> criteria that shall be used for comparison, and if the best
%           model is the min or max regarding that criteria. The patients
%           table has to have a column named as the compCriteria

%OUT:
%   BestMeasure
%           -> The compCriteria value of the best model
    
    if minOrMax=="min"
        bestMeasure=par_table(par_table.(compCriteria)==min(par_table.(compCriteria)),compCriteria).Variables;

    elseif minOrMax=="max"
        bestMeasure=par_table(par_table.(compCriteria)==max(par_table.(compCriteria)),compCriteria);

    end

    
end