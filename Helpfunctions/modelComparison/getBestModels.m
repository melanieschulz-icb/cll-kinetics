function [goodModels] = getBestModels(patient, folder, compCriteria, threshold,sheet)
    [allNums, txt_x_y_z,~] = xlsread(strcat(pwd,folder,['patient_',int2str(patient),'/Params_Pat',int2str(patient)]),'sheet'=sheet); 
    

    txt_x_y_z=txt_x_y_z(1,:);
    [~, model]=find(~cellfun(@isempty, regexp(txt_x_y_z, 'model')));
    model=model(1);

    [~, measure]=find(~cellfun(@isempty, regexp(txt_x_y_z, compCriteria)));
    measure=measure(1);
    
    %find the best model
    [bestMeasure,i]=min(allNums(:,measure));
%     bestModel=allNums(i,model);
    
    %models with measure difference <goodThreshold
    goodModels=unique(allNums(allNums(:,measure)-bestMeasure<=threshold,model))';
%     [~,i]=unique(goodModels(:,model));
%     goodModels=goodModels(i,:);
    
end