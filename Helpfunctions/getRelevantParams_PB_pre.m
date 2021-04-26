function par_positions = getRelevantParams_PB_pre(model)
    if isempty(model) || model==-1
        par_positions=[];
    elseif ~isempty(find(1:6==model,1))
        par_positions=1:2;
    elseif model==7
        par_positions=1;
    else
        error("PB_pre model not specified");
    end
end