function par_positions = getRelevantParams_comb(model,combine_pre_post)
    par_with_pb0=getRelevantParams(model);
    if combine_pre_post
        names=getVarNames(model);

        ind_pb0=find(names(getRelevantParams(model))=="PB_0");

        par_positions=par_with_pb0([1:ind_pb0-1,ind_pb0+1:end]);
    else
        par_positions=par_with_pb0;
    end
end