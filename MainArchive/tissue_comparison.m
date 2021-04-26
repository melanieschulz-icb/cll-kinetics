XdataComp=[];
cll_volume=166;
for patient=[1:11,13:30]
    [~,XdataComplex]=getTissueDataComplex(patient,cll_volume);
    [~, XdataSimple]=getTissueData(patient, cll_volume);
    
    dif=XdataComplex(2)/XdataSimple(2);
    XdataComp=[XdataComp;patient, XdataComplex(2), XdataSimple(2), dif, 1/dif];
end

writetable(array2table(XdataComp,'VariableNames',{'patient', 'sec_ts_estimated', 'sec_ts_raw','est/raw','raw/est'}),'Tissue_comparison.xlsx');
    
   