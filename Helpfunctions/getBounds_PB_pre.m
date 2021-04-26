function [lb, ub,lb_sd, ub_sd, par_init, sdXYZ] = getBounds_PB_pre(PB_One, model)
%Set ParameterBounds and Initialization_Params            
            
        %all available Params so far, from Model 1-5
        %HB_0, HB_max, k (growing), d (decrease)
        lb_all = [0.5*PB_One, 1e-5]; 
        ub_all = [1.5*PB_One, 20]; 


        lb_sd = 1e-10;
        ub_sd = 100;

        %construct start (later add multi start)
        %all available Params so far, from Model 1-5
        %start parameters
        par_init_all = [PB_One 0.1]; 
        sdXYZ=[2];

        %select relevant Bounds/Params, based on Model
%             lb=lb_all(getRelevantHBParams(model));
%             ub=ub_all(getRelevantHBParams(model));
%             par_init=par_init_all(getRelevantHBParams(model));

        lb=lb_all(getRelevantParams_PB_pre(model));
        ub=ub_all(getRelevantParams_PB_pre(model));
        par_init=par_init_all(getRelevantParams_PB_pre(model));


end