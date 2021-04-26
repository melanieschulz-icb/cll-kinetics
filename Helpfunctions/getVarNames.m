function [varNames] = getVarNames(model)
   
    varNames=["d_PB", "d_LT","m_LT_PB","LT_0","PB_0","c_LT","d_BM","m_BM_PB","BM_0","c_BM","m_LT_BM","c_PB"]; 
    %           1       2       3       4       5       6       7       
           %    9       10      11      12
   
    if model>181 && model <=309 || model>1021&&model<=1037 || model>2021&&model<=2037
        varNames(11)="m_BM_LT";
    elseif model >309 && model <=437 || model>1037&&model<=1053 || model>2037&&model<=2053
        varNames(8)="m_PB_BM";
    elseif model > 437 && model <=565 || model>1053&&model<=1069 || model>2053&&model<=2069
        varNames(3)="m_PB_LT";
        varNames(11)="m_BM_LT";
    elseif model >565 && model <=693 || model>1069&&model<=1085 || model>2069&&model<=2085
        varNames(3)="m_PB_LT";
        varNames(8)="m_PB_BM";
    elseif model >693 && model <=821 || model>1085&&model<=1101 || model>2085&&model<=2101
        varNames(3)="m_PB_LT";
        varNames(8)="m_PB_BM";
        varNames(11)="m_BM_LT";
    end
    %if model >1000
        varNames(6)="LT_healthy";
        varNames(10)="BM_healthy";
        varNames(12)="Blood_healthy";
    %end
end