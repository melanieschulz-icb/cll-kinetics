function [lb, ub,lb_sd, ub_sd, par_init, sdXYZ] = getBounds(xOne, yOne, model)
%Set ParameterBounds and Initialization_Params            
            
            %all available Params so far, from Model 1-5
            %d2, d1,m1,x0,y0,c1, d3, m2,z0,c2,m3
            lb_all = [1e-10 2e-15 3e-15 0.5 0.5 1e2 4e-15 9e-16 1e3 1e3 8e-16]; 
            ub_all = [5 20 2 1.5 1.5 1e18 3 4 1e20 1e20 1]; 
            

            lb_sd = [1e-10 1e-10 1e-10];
            ub_sd = [100 100 100];
            
            %construct start (later add multi start)
            %all available Params so far, from Model 1-5
            %start parameters
            par_init_all = [0.1 0.2 0.25 1e10 1e10 1e10 0.4 0.6 1e10 1e10 0.25]; 
            sdXYZ=[20 20 20];

            %select relevant Bounds/Params, based on Model
            lb=lb_all(getRelevantParams(model));
            ub=ub_all(getRelevantParams(model));
            par_init=par_init_all(getRelevantParams(model));
            
            %Adapt x0, y0 -> take first measurements as base, allow to vary
            %by x percent
            varNames=getVarNames();
            ind_x0=find(varNames(getRelevantParams(model))=="x0");
            ind_y0=find(varNames(getRelevantParams(model))=="y0");
            
            lb(ind_x0) = lb(ind_x0) * xOne ; ub(ind_x0) = ub(ind_x0) * xOne; par_init(ind_x0) = xOne;          
            lb(ind_y0) = lb(ind_y0) * yOne ; ub(ind_y0) = ub(ind_y0) * yOne; par_init(ind_y0) = yOne;
end