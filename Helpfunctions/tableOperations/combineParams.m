function combineParams(sourceFolder,saveFolder,patients)
    patient=patients(1);
    %Par_table = readtable(strcat(pwd,sourceFolder,['patient_',int2str(patient),'/Params_Pat',int2str(patient)]));
    Par_table=load(strcat(pwd,"/",sourceFolder,['patient_',int2str(patient),'/Params_Pat',int2str(patient),'.mat']));
    Par_table=Par_table.tab;
    for patient=patients(2:end)
        Par_table = [Par_table;load(strcat(pwd,"/",sourceFolder,['patient_',int2str(patient),'/Params_Pat',int2str(patient),'.mat'])).tab];
    end
    writetable(Par_table, strcat(pwd,"/",saveFolder,"allPatients.xlsx"));
    save(strcat(pwd,"/",saveFolder,"allPatients"),'Par_table');
end

