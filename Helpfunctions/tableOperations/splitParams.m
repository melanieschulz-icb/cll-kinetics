function splitParams(sourceTable,saveFolder)

    allPat_table = readtable(strcat(pwd,sourceTable));
    for patient=unique(allPat_table.ACC)'
        pat_table=allPat_table(allPat_table.ACC==patient,:);
        mkdir(strcat(saveFolder,'patient_',int2str(patient),'/'));
        writetable(pat_table,strcat(saveFolder,"patient_",int2str(patient),"/Params_Pat",int2str(patient),".xlsx"));
    end

end