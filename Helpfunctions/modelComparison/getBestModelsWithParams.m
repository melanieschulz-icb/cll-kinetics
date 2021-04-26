function [BestModels] = getBestModelsWithParams(patients, inputFolder, compCriteria,minOrMax,threshold)
%getBestModelsWithParams
%Returns Best Models with parameters, where best Models is defined by a
%given compCriteria and a threshold, defining the range where two models
%are considered to be equally good


%params, IN:
%   patients
%           -> patients that shall be taken intwo account
%   inputFolder
%           -> basefolder -> for each patient that shall be taken into
%           accunt, a patient folder must exist within the basefolder,
%           within that patient folder, the patient table must be placed
%   compCriteria, minOrMax
%           -> criteria that shall be used for comparison, and if the best
%           model is the min or max regarding that criteria. The patients
%           table has to have a column named as the compCriteria
%   threshold
%           -> models are returned that differ at max thr. from the best
%           model
%OUT:
%   BestModels
%           -> One table containing the best models of all patients

    counter=1;
    
    for patient=patients
        par_table=readtable(strcat(pwd,inputFolder,'patient_',int2str(patient),'/Params_Pat',int2str(patient)),"Sheet",2);
        bestMeasure=getBestMeasure(par_table, compCriteria,minOrMax);
        if minOrMax=="min"
            bestModels=par_table(par_table{:,compCriteria}-bestMeasure<=threshold,:);
        elseif minOrMax=="max"
            bestModels=par_table(bestMeasure-par_table{:,compCriteria}<=threshold,:);
        end
        if counter==1
            BestModels=bestModels;
        else
            BestModels=union(BestModels,bestModels);
        end
        counter=counter+1;
    end

    
end