function [modelParams] = getModelParams(patients, inputFolder, models,sheet)
%getModelParams
%Returns parameters of the specified models



%params, IN:
%   patients
%           -> patients that shall be taken intwo account
%   inputFolder
%           -> basefolder -> for each patient that shall be taken into
%           accunt, a patient folder must exist within the basefolder,
%           within that patient folder, the patient table must be placed
%   models
%           -> models that shall be returned
%OUT:
%   modelParams
%           -> One table containing the params of all patients/models

    modelParams=table;
    
    for patient=patients
        par_table=readtable(strcat(pwd,inputFolder,'patient_',int2str(patient),'/Params_Pat',int2str(patient)),'sheet',sheet);
        for model=models
            modelPatParams=par_table(par_table.model==model,:);
            modelPatParams=modelPatParams(modelPatParams.bic==min(modelPatParams.bic),:);
            modelParams=[modelParams;modelPatParams(1,:)];
        end
    end
        

    
end