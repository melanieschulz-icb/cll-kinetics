function [lb, ub,lb_sd, ub_sd, par_init, sdXYZ] = getBoundsHB(HBOne, model)
%Set ParameterBounds and Initialization_Params            
            
            %all available Params so far, from Model 1-5
            %HB_0, HB_max, k (growing), d (decrease)
            lb_all = [0.5*HBOne, 0.5*HBOne, 1e-5, 1e-5]; 
            ub_all = [1.5*HBOne, 20, 10^3, 5]; 
            

            lb_sd = 1e-10;
            ub_sd = 100;
            
            %construct start (later add multi start)
            %all available Params so far, from Model 1-5
            %start parameters
            par_init_all = [HBOne HBOne 5 0.01]; 
            sdXYZ=[2];

            %select relevant Bounds/Params, based on Model
            lb=lb_all(getRelevantHBParams(model));
            ub=ub_all(getRelevantHBParams(model));
            par_init=par_init_all(getRelevantHBParams(model));


            
            %Adapt x0, y0 -> take first measurements as base, allow to vary
            %by x percent
%             varNames=getVarNames(model);
%             ind_x0=find(varNames(getRelevantParams(model))=="LT_0");
%             ind_y0=find(varNames(getRelevantParams(model))=="PB_0");
%             
%             lb(ind_x0) = lb(ind_x0) * xOne ; ub(ind_x0) = ub(ind_x0) * xOne; par_init(ind_x0) = xOne;          
%             lb(ind_y0) = lb(ind_y0) * yOne ; ub(ind_y0) = ub(ind_y0) * yOne; par_init(ind_y0) = yOne;
end