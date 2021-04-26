function par_positions = getRelevantHBParams(model)
    if isempty(model)
        par_positions=[];
    elseif model == 1 || model == 4 || model ==5
        par_positions=1:3;
    elseif model == 2 || model==3 || model == 6 || model==9 || model==10
        par_positions=1:4;
    elseif model==7 || model == 8
        par_positions=[1,3,4];
    else
        error("HB model not specified");
    end
end