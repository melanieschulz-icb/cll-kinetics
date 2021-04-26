function [y_out] = solution_PB_pre(par,tspanY_pre, model)
    
    tspanY_pre_shifted=tspanY_pre-min(tspanY_pre);
    if isempty(model)
        error("please specify a model");
    end
    
    %put params into modelParams
    if model==-1
        y_out=[];
    elseif model == 1
        pb_0=par(1);
        k=par(2);

        y_out=pb_0*exp(k*tspanY_pre_shifted);
    elseif model == 2
        pb_0=par(1);
        k=par(2);

        y_out=pb_0*exp(-k*tspanY_pre_shifted);
    elseif model == 3
        pb_0=par(1);
        k=par(2);

        y_out=pb_0*2.^(k*tspanY_pre_shifted);
    elseif model == 4
        pb_0=par(1);
        k=par(2);

        y_out=pb_0*exp(k*tspanY_pre);
    elseif model == 5
        pb_0=par(1);
        k=par(2);

        y_out=pb_0*exp(-k*tspanY_pre);
    elseif model == 6
        pb_0=par(1);
        k=par(2);

        y_out=pb_0*2.^(k*tspanY_pre);   
    elseif model == 7
        pb_0=par(1);
        k=0;

        y_out=pb_0*exp(k*tspanY_pre_shifted);
    else
        error("Selected PB_pre model not implemented");
    end
    
end