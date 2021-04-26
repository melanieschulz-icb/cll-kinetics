function [lb, ub,lb_sd, ub_sd, par_init, sdXYZ] = getBounds(xOne,xEnd,yOne,zOne,model,patient,combine_pre_post,bounds)
%Set ParameterBounds and Initialization_Params            
            
            if bounds == 1 %broad bounds
                %all available Params so far, from Model 1-5
                %d2, d1,m1,x0,y0,c1, d3, m2,z0,c2,m3,c3
                lb_all = [1e-8 1e-8 1e-8 0.5 0.5 1e2 1e-8 1e-8 0.2 1e3 1e-8 1e3]; 
                ub_all = [20 20 20 1.5 1.5 1e15 20 20 8 1e15 20 1e12];
            elseif bounds == 2 %rates Burger, c/X0/Y0 Burger
                %all available Params so far, from Model 1-5
                %d2, d1,m1,x0,y0,c1, d3, m2,z0,c2,m3,c3
                lb_all = [0.001 0.001 1e-4 0.9 0.9 1e-3 1e-8 1e-8 0.2 1e3 1e-8 1e3]; 
                ub_all = [0.5 20 0.1 1.1 1.1 xEnd 0.8 0.8 8 1e15 0.8 1e12];
            elseif bounds == 3 %rates Burger, c/x0,y0 wodarz
                volume=estimateBloodVolume(patient);
                lb_all = [0.001 0.001 1e-4 0.85 0.85 1e-3 1e-8 1e-8 0.2 1e3 1e-8 1e3]; 
                ub_all = [0.5 20 0.1 1.15 1.15 80*volume 0.8 0.8 8 1e15 0.8 1e12];
            elseif bounds ==4 %rates wodarz, c/x0,y0 burger
                lb_all = [0.001 0.001 1e-4 0.9 0.9 1e-3 1e-8 1e-8 0.2 1e3 1e-8 1e3]; 
                ub_all = [0.1 0.1 0.5 1.1 1.1 xEnd 0.8 0.8 8 1e15 0.8 1e12];
            elseif bounds == 5 %rates wodarz, c/x0,y0 wodarz
                volume=estimateBloodVolume(patient);
                lb_all = [0.001 0.001 1e-4 0.85 0.85 1e-3 1e-8 1e-8 0.2 1e3 1e-8 1e3]; 
                ub_all = [0.1 0.1 0.5 1.15 1.15 80*volume 0.8 0.8 8 1e15 0.8 1e12];
            else
                error("Boundaries not specified");
            end
%             lb_all = [1e-2 1e-2 1e-2 0.5 0.5 1e10 1e-2 1e-2 1e8 1e8 1e-2 1e5]; 
%             ub_all = [1 1 1 1.5 1.5 1e12 1 1 1e14 1e14 1 1e10];
            
            

            lb_sd = [1e-10 1e-10 1e-10];
            ub_sd = [100 100 100];
            
            %construct start (later add multi start)
            %all available Params so far, from Model 1-5
            %start parameters

           
%             par_init_all = [0.1 0.15 0.20 1e10 1e10 1e10 0.35 0.30 1e10 1e10 0.25 1e8]; 
            par_init_all = [0.1 0.1 0.1 1e10 1e10 1e10 0.1 0.1 1e10 1e10 0.1 2];             
%               rates=par_init_all([1,2,3,7,8,11]);
%             par_init_all([1,2,3,7,8,11])=rates(randperm(6));
            
            sdXYZ=[20 20 20];
            
            %turn bounds for normal relative PB counts into absolute values
%             lb_all(end)=turnToAbsoluteCounts(lb_all(end),patient);
%             ub_all(end)=turnToAbsoluteCounts(ub_all(end),patient);
             par_init_all(end)=turnToAbsoluteCounts(par_init_all(end),patient);
            
            
            %select relevant Bounds/Params, based on Model
            lb=lb_all(getRelevantParams_comb(model,combine_pre_post));
            ub=ub_all(getRelevantParams_comb(model,combine_pre_post));
            par_init=par_init_all(getRelevantParams_comb(model,combine_pre_post));
            
            %Adapt x0, y0 -> take first measurements as base, allow to vary
            %by x percent
            varNames=getVarNames(model);
            ind_x0=find(varNames(getRelevantParams_comb(model,combine_pre_post))=="LT_0");
            ind_y0=find(varNames(getRelevantParams_comb(model,combine_pre_post))=="PB_0");
            ind_z0=find(varNames(getRelevantParams_comb(model,combine_pre_post))=="BM_0");
            
            lb(ind_x0) = lb(ind_x0) * xOne ; ub(ind_x0) = ub(ind_x0) * xOne; par_init(ind_x0) = xOne;          
            lb(ind_y0) = lb(ind_y0) * yOne ; ub(ind_y0) = ub(ind_y0) * yOne; par_init(ind_y0) = yOne;
            lb(ind_z0) = lb(ind_z0) * zOne * estimateBMCells(patient); ub(ind_z0)=ub(ind_z0) * zOne * estimateBMCells(patient);...
                par_init(ind_z0)=zOne*estimateBMCells(patient);

end