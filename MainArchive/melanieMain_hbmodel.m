%% melanie_main_Burger_sls
% can be used to fit the Burger data
%-> use lsq with normalization
close all; clearvars;

%% Choose parameter

%models = [3:51,52,55,58,59,61,63,65]; %:66;
%models = [1001:1006];
%models=182:309;
models=2;

%patients you want to fit; 12 will be excluded as ibrutinib dose was paused
%(available patients: 1-30)
patients=[1:2];
runns=1:2;



%if one should normaize when using the least square optimizer
normalize_lls=0;

%assume multiplicative noise
mult_noise=0;

%if optimizing for sd
opt_sd=1;

%if computing/optimizing joint or seperate sd for each compartment
seperate_sd=0;

%if one assumes that vanishing CLL cells in BM are instantly replaced by
%healthy cells
refill=1;

%where to save output
folderbase='/Figures/12_02/';




%% Read raw data

[num, txt, ~] = xlsread('Heavywater-all counts-decoded.xlsx');
[numPlat, txtPlat, ~] = xlsread('Heavywater-all counts-decoded.xlsx','Platelets');
[chars, charDesc,~]=xlsread('Heavy Water patients characterists.xlsx');
%% Prepare Variables for storing the opt. results
Xmulti=[];

%% Perform fitting and plot the results
for patient=patients
    mkdir([folderbase(2:end),'patient_',int2str(patient),'/']);
    folder=[folderbase,'patient_',int2str(patient),'/'];
    mkdir([folder(2:end),'/plots/']);
    plotFolder=[folder,'/plots/'];
    
    for model=models
        for runn=runns
            if patient==12
                'Skipping patient_nr 12 as therapy was paused in between';
            else
                'Fitting Patient_nr '+ string(patient);
               clearvars -except num txt chars charDesc patient use_constraint patients cll_volume Xmulti Xllh Errormulti Errorllh runn runns tss_measurements x normalize_lls tss_measurements use_complex_tss weight mult_noise opt_sd folder folderbase plotFolder model models seperate_sd refill params numPlat txtPlat Xmulti
               close all; 
                %% optimization using integral

                %% Set data

                %use available mess points as time span
                %use available mess values as true data

%                 [startDate, ~] = aggBloodData(num, txt, chars, charDesc, patient,false);
                [~,tspanYPlats,Platelets,PlateletsScaled, HB,HBScaled, Neutroph,NeutrophScaled, RBC, RBCScaled] = getBloodCounts(numPlat, txtPlat, chars, charDesc,patient, false);
                [startDate, tspanY, ydata, ty_outlier, y_outlier] = aggBloodData(num, txt, chars, charDesc, patient);
                tspanYPlats=tspanY';
                HB=ydata';

                %% Set Bounds

                [lb, ub,lb_sd, ub_sd, par_init, sdXYZ] = getBoundsHB(HB(1), model);  

                n_starts=60;


                %% test log pars



                lb=log10(lb);
                ub=log10(ub);
                par_init=log10(par_init);

                %% Get latin hypercube sampled starts

                [par_starts,~] = lhsdesign_modified(n_starts,lb,ub);

                startpoints = CustomStartPointSet(par_starts(:,:));

                %% define objective function.
                %if normalize_lls=true, use the appropriate objective function
                if normalize_lls==true
                    if mult_noise==true
                            lsqfun=@(par)diff_log_norm_summed(par,xdata,ydata, tspanX, tspanY);
                    else
                        lsqfun = @(par)diff_integral_norm_summed(10.^par,xdata,ydata,zdatanorm,tspanX, tspanY,tspanZ);
                    end
                else
                    if mult_noise==true
                        if opt_sd==false
                            lsqfun=@(par)diff_log_summed(par,xdata,ydata,tspanX, tspanY);
                        else
                            lsqfun=@(par)diff_log_summed_sd(10.^par,xdata,ydata,zdata,tspanX, tspanY,tspanZ, model, refill, seperate_sd);
                        end
                    else
                        if opt_sd==false
                            lsqfun = @(par)diff_integral_summed(par,xdata,ydata,tspanX, tspanY);
                        else
                            lsqfun = @(par)diff_integral_summed_sd(10.^par,[],[],[],[],[],[],[],refill,seperate_sd,HB,tspanYPlats,model);
                        end
                    end 
                end

                %% create OptimProblem

                oldOptions = optimoptions('fmincon');
                options = optimoptions(oldOptions,'MaxFunctionEvaluations',6000,'MaxIterations',3000,'OptimalityTolerance',1e-8);

                problem = createOptimProblem('fmincon','x0',par_init,'objective',lsqfun,...
                'lb',lb,'ub',ub);


                %% run the optimization problem 
                ms = MultiStart('PlotFcns',{@gsplotbestf, @errorlandscape},'FunctionTolerance',1e-8);
%                 ms = MultiStart('FunctionTolerance',1e-8);
                [xmulti,errormulti,exitflag,output,solutions] = run(ms,problem,startpoints);
                xmulti=10.^xmulti;
%                 if mult_noise==1
%                     [sd_x, sd_y, sd_z]=compute_mult_sd(xmulti,xdata,ydata,zdata,tspanX, tspanY, tspanZ, model, refill,seperate_sd);
%                     likelihood=likelihood_mult_noise(xmulti,xdata,ydata,zdata,tspanX, tspanY,tspanZ, model, refill,seperate_sd);
% 
%                 else
%                     [sd_x, sd_y, sd_z]=compute_add_sd(xmulti,xdata,ydata,zdata,tspanX, tspanY, tspanZ, model, refill,seperate_sd);
%                     likelihood=likelihood_add_noise(xmulti,xdata,ydata,zdata,tspanX, tspanY,tspanZ, model,refill ,seperate_sd);
%                 end
%                 bic=(length(xmulti)+1+seperate_sd*2)*log(length([xdata, ydata, zdata]))-2*log(likelihood);
%                 %bic=(length(xmulti))*log(length([xdata, ydata]))-2*log(likelihood);
%                 aicr=2*(length(xmulti)+1+seperate_sd*2)-2*log(likelihood);



                %% save params and values
                 storeParams=zeros(1,4);
                 storeParams(getRelevantHBParams(model))=xmulti;
%                 storeParams(13:15)=[sd_x, sd_y, sd_z];
%                 Xmulti=[Xmulti;[model,patient,storeParams,bic, aicr,likelihood]];
                 Xmulti=[Xmulti;[model,patient,storeParams,errormulti]];


%                 %% plot data
% 
                tspanPlats = 0:1:600;
                
                [hb_opt] = solutionHB(xmulti,tspanPlats,[],model);

                 figure;
                subplot(2,1,1)
                plot(tspanYPlats, HB, '*k'); hold on;
                %outliner not taken into account


                %plot(tspan, y_true, 'color',[0 0.5 0]); hold on;
                plot(tspanPlats, hb_opt, '--m',xlabel('days'),ylabel('HB value')); hold on;

                legend('measurements','solution lsq')

                plotTitle=title('ACC_'+string(patient)+'_model_'+string(model)+'_HB_run_'+string(runn),'Interpreter','none');
                %plotTitle.Color='red';
                plotTitle.FontSize=16;


                %save figure as PDF
                f=gcf;
                exportgraphics(f,[pwd,strcat(plotFolder,plotTitle.String,'.png')],'Resolution',300)

          

                %% Profile Likelihood Analysis

    %             if mult_noise==1 && opt_sd==1
    %                 original_lsqfun=@(par)diff_log_summed_sd(par,xdata,ydata,zdatanorm,tspanX, tspanY, tspanZ);
    %                 [parval, parres] = Profiles_lsq(original_lsqfun, xmulti, 10*n_starts, [10.^lb,1e-8,1e-8, 1e-8], [10.^ub,0.2,1,1]);
    %             elseif mult_noise==0 && opt_sd==1
    %                 original_lsqfun=@(par)diff_integral_summed_sd(par,xdata,ydata,zdatanorm,tspanX, tspanY,tspanZ);
    %                 [parval, parres] = Profiles_lsq(original_lsqfun, xmulti, 10*n_starts, [10.^lb,10e4,10e4, 10e4], [10.^ub,10e14,10e14, 10e4]);
    %             else
    %                 [parval, parres] = Profiles_lsq(lsqfun, log10(xmulti), 10*n_starts, lb, ub);
    %             end
    %             
    % 
    %             figure;
    %             plotTitle=sgtitle('parval_ls_ACC_'+string(patient)+'_'+string(runn),'Interpreter','none');
    %             %plotTitle.Color='red';
    %             plotTitle.FontSize=16;
    %             subplot(3,5,1)
    %             plot(parval(1,:),parres(1,:));xlabel('parameter value');ylabel('sum of squared differences');title('d2');
    %             subplot(3,5,2)
    %             plot(parval(2,:),parres(2,:));xlabel('parameter value');ylabel('sum of squared differences');title('d1');
    %             subplot(3,5,3)
    %             plot(parval(3,:),parres(3,:));xlabel('parameter value');ylabel('sum of squared differences');title('m');
    %             subplot(3,5,4)
    %             plot(parval(4,:),parres(4,:));xlabel('parameter value');ylabel('sum of squared differences');title('x0');
    %             subplot(3,5,5)
    %             plot(parval(5,:),parres(5,:));xlabel('parameter value');ylabel('sum of squared differences');title('y0');
    %             subplot(3,5,6)
    %             plot(parval(6,:),parres(6,:));xlabel('parameter value');ylabel('sum of squared differences');title('c');
    % 
    %             subplot(3,5,7)
    %             plot(parval(7,:),parres(7,:));xlabel('parameter value');ylabel('sum of squared differences');title('d3');
    %             subplot(3,5,8)
    %             plot(parval(8,:),parres(8,:));xlabel('parameter value');ylabel('sum of squared differences');title('m2');
    %             subplot(3,5,9)
    %             plot(parval(9,:),parres(9,:));xlabel('parameter value');ylabel('sum of squared differences');title('z0');
    %             subplot(3,5,10)
    %             plot(parval(10,:),parres(10,:));xlabel('parameter value');ylabel('sum of squared differences');title('c2');
    % 
    %             if opt_sd==true
    %                subplot(3,5,11)
    %                plot(parval(11,:),parres(11,:));xlabel('parameter value');ylabel('sum of squared differences');title('sdx');
    %                subplot(3,5,12)
    %                plot(parval(12,:),parres(12,:));xlabel('parameter value');ylabel('sum of squared differences');title('sdy');
    %                subplot(3,5,13)
    %                plot(parval(13,:),parres(13,:));xlabel('parameter value');ylabel('sum of squared differences');title('sdz');
    % 
    %             end
    % 
    % 
    %             %save figure as PDF
    %             f=gcf;
    %             exportgraphics(f,[pwd,strcat(folder,plotTitle.String,'.png')],'Resolution',300)
    %             savefig(f,[pwd,strcat(folder,plotTitle.String)],'compact')
            end 
        end
    end
    
 %% Store Parameters and errors
%     varNames=getVarNames(model);
%     if opt_sd==false
%         writetable(array2table(Xmulti, 'VariableNames',["model","ACC",varNames]),[pwd,strcat(folder,'Params.xlsx')]);
%     else
%         writetable(array2table(Xmulti, 'VariableNames',["model","ACC",varNames,"sdx","sdy","sdz", "bic", "aicr","likelihood"]),[pwd,strcat(folder,['Params_Pat',int2str(patient),'.xlsx'])]);
%     end
%writetable(array2table(Xmulti, 'VariableNames',["model","ACC","HB_0", "HB_max","k","negLogLik"]),[pwd,strcat(folder,['Params_Pat',int2str(patient),'.xlsx'])]);

end

writetable(array2table(Xmulti, 'VariableNames',["model","ACC","HB_0", "HB_max","growth_rate","death_rate","negLogLik"]),[pwd,strcat(folder,['Params_Pat',int2str(patient),'.xlsx'])]);





function [parval, parres] = Profiles_lsq(lsqfun, par, num, lb, ub)
%calculate the likelihood for a certain parameter at num values between 
%the lb and the ub.

%Input:
%   llhfun      function handle of the likelihood function
%   par         optimized parameters
                                %  , parnum parnum      parameter to optimize
%   num         number of points to calculate in the profile
%   lb          lb of the problem
%   ub          ub of the problem

%Output:
%   parval      values of the parameter
%   parres      llh value for the parameter
    
% can also look into MATLAB profile function
    originalpar = par; %save original parameters

    parnum = length(par);
    
    parres = zeros(parnum,num) ;
    parval = zeros(parnum,num) ;

    for j = 1:parnum
        
        parval(j,:) = linspace(lb(j), ub(j),num);
        
        for i = 1:num
            par(j) = parval(j,i);
            
            parres(j,i) = lsqfun(par);
        end
        
        par = originalpar; %reset to original input parameters
        
    end 
    

end