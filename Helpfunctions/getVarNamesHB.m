function [varNames] = getVarNamesHB(hbmodel)
   
    varNames=["HB_0", "HB_norm", "groth", "death"];
    varNames=varNames(getRelevantHBParams(hbmodel));
    
end