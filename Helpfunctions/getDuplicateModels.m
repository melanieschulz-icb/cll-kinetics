function [duplModels]=getDuplicateModels(model)
    if model ==3
        duplModels=readtable(strcat(pwd,'/duplModels3CompsSimple.xlsx'));
    elseif model == 1003
        duplModels=readtable(strcat(pwd,'/duplModels3CompsSeparation.xlsx'));
    elseif model == 2003
        duplModels=readtable(strcat(pwd,'/duplModels3CompsSeparation.xlsx'));
        duplModels.duplModels=duplModels.duplModels+1000;
    else
        error('choose 3 or 1003 or 2003 to find duplmodels')     
    end
end
        