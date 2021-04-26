%% melanie_main_Burger_sls
% can be used to fit the Burger data
%-> use lsq with normalization
close all; clearvars;

%% Choose parameter

%patients you want to fit; 12 will be excluded as ibrutinib dose was paused
%(available patients: 1-30)
patients=4;
runns=1;
%if to use the constraint d2<d1
use_constraint=false;

%assumed volume of a CLL cell (in fl)
cll_volume=166;

%tissue measurements to be taken into account (max 3 available, Burger
%takes only one)
tss_measurements=1;

%if one should use the complex computation of the tissue cll burden
use_complex_tss=0;

%if one should normaize when using the least square optimizer
normalize_lls=1;

%weight for the second tissue measurement (if second measurement is taken
%into account)
weight=1;

%assume multiplicative noise
mult_noise=0;

%if optimizing sd, when assuming multiplicative noise
opt_sd=0;

params=[tss_measurements, use_complex_tss,normalize_lls, weight, mult_noise, opt_sd];

%folder to save figures
folder='/Figures/08_23/';

%save settings
writetable(array2table(params, 'VariableNames',{'tss_measurements', 'use_complex_tissue','normalize_lls', 'weight', 'mult_noise', 'opt_sd'}),[pwd,strcat(folder,'Settings.txt')]);

%% Read the PB data
[num, txt, ~] = xlsread('Heavywater-all counts-decoded.xlsx');
[chars, charDesc,~]=xlsread('Heavy Water patients characterists.xlsx');
%% Prepare Variables for storing the opt. results
Xmulti=[];
Errormulti=[];


%% Perform fitting and plot the results
for patient=patients
    for runn=runns
        if patient==12
            'Skipping patient_nr 12 as therapy was paused in between';
        else
            'Fitting Patient_nr '+ string(patient);
           clearvars -except num txt chars charDesc patient use_constraint patients cll_volume Xmulti Xllh Errormulti Errorllh runn runns tss_measurements x normalize_lls tss_measurements use_complex_tss weight mult_noise opt_sd folder
           close all;
            %% optimization using integral

            %settings

            rng(0);
            parameters.name = {'d2','d1','m','x0','y0','c'};

            parameters.number = length(parameters.name);
  
            lb = [1e-5 0 0 0.5 0.5 1e8]; 
            ub = [0.5 0.5 0.5 1.5 1.5 1e12]; 
 
%             end
            lb_sd = [1e-7 1e-7];
            ub_sd = [1e14 1e13];
            

            n_starts = 800;


            %construct start (later add multi start)
            par_init = [0.1 0.1 0.000 0 0 1e10]; %start parameters
            sdXY=[1e10 1e10];

            %% Set data

            %use available mess points as time span
            %use available mess values as true data
            

            [startDate,tspanY, ydata]=getPatientData(num, txt, chars, charDesc, patient);
            tspanY=tspanY';
            ydata=ydata';
            ty_outlier=[];
            y_outlier=[];
            
            %fix data
            if patient == 3
                outlier=4;
                ty_outlier=tspanY(outlier);
                tspanY=[tspanY(1:outlier-1),tspanY(outlier+1: end)];
                y_outlier=ydata(outlier);
                ydata=[ydata(1:outlier-1),ydata(outlier+1: end)];
            elseif patient == 7
                %datum falsch eingetragen (+1Jahr)
                tspanY(5)=tspanY(5)-365;
            elseif patient == 8
                outlier=4;
                ty_outlier=tspanY(outlier);
                tspanY=[tspanY(1:outlier-1),tspanY(outlier+1: end)];
                y_outlier=ydata(outlier);
                ydata=[ydata(1:outlier-1),ydata(outlier+1: end)];
            end
            
            if use_complex_tss==1    
                [tspanXavailable,Xdata]=getTissueDataComplex(patient,cll_volume);
            else
                [tspanXavailable, Xdata]=getTissueData(patient, cll_volume);
            end
            
            tspanXavailable=(tspanXavailable-startDate)';
            tspanXavailable(tspanXavailable<0)=0;
            tspanX=tspanXavailable(1:tss_measurements);
            xdata=Xdata(1:tss_measurements)';
            if patient==24
                %second measurement missing!!!
                tspanX=tspanX(1:1);
                xdata=xdata(1:1);
            end

            %use first measured values as initial parameter guess for x0, y0
            %allow them to vary +-10%

            lb(4) = lb(4) * Xdata(1) ; ub(4) = ub(4) * Xdata(1); par_init(4) = Xdata(1);          
            lb(5) = lb(5) * ydata(1) ; ub(5) = ub(5) * ydata(1); par_init(5) = ydata(1);
            
            %use the last tissue datapoint as upper bound for c
            [~,simple]=getTissueData(patient, cll_volume);
            %ub(6)=simple(end);
            %ub(6)=1.1*Xdata(end);
            ub(6)=1e14;
            %ub(6)=Xdata(end);
            
            parameters.min = lb;
            parameters.max = ub;
            %% account for sd
            if opt_sd==true && mult_noise==0
                lb=[lb, lb_sd];
                ub=[ub, ub_sd];
                par_init=[par_init,sdXY];
            end
                

            %% optimization for integral

            %% define objective function.

            
            objectiveFunction = @(par)diff_integral_norm_summed(par,xdata,ydata,tspanX, tspanY, weight);
            
            %% create OptimProblem
            %if use_constraint==true,
            %add the linear inequality condition d_1>d_2, that is d2-d1<0
            %that is par(1)-par(2)<0
 
            optionsPesto = PestoOptions();
            optionsPesto.obj_type = 'log-posterior';
            optionsPesto.n_starts = 200;
            optionsPesto.comp_type = 'sequential';
            optionsPesto.mode = 'visual';

%            optionsPesto.plot_options.add_points.logPost = objectiveFunction(theta_true);
            %optionsPesto.plot_options.add_points.prop = nan(properties.number,1);
            %for j = 1 : properties.number
            %    optionsPesto.plot_options.add_points.prop(j) = properties.function{j}(optionsPesto.plot_options.add_points.par);
            %end
            
            if strcmp(optionsPesto.comp_type, 'parallel') && (n_workers >= 2)
                parpool(n_workers); 
            else
                optionsPesto.comp_type = 'sequential';
            end
            
            parameters = getMultiStarts(parameters, objectiveFunction, optionsPesto);
            parameterss = getParameterProfiles(parameters, objectiveFunction, optionsPesto);
        end
    end
end

%% Functions



% function [sensitivity] = diff_integral_summed_sd(par,xdata,ydata,tspanX, tspanY, weight)
% %Calculates the differences between the data and simulation output,
% %normalized via the l2_norm of the measurements
% %can be used as input for lsqnonlin/least squares optimizer
% %outputs:                          
% %   sensitivity = output feature, differences of both outputs
% 
% %inputs:
% %   tspan       timespan for integral
% %   par =       adapted input parameters
% %   data =      data vectors for x and y respectively
% 
% % determine output for adapted parameters
% [x_out, y_out] = intLea(par(1:6), tspanX, tspanY);
%  
% %featX = sqrt((1/weight))*(xdata - x_out)/par(7); %calculate differences
% featX = (xdata - x_out)/par(7);
% featY = (ydata - y_out)/par(8); %calculate differences
% %featY(1)=sqrt((1/weight))*featY(1);
% sensitivity = length(xdata)*log(2*pi*par(7)^2)+length(ydata)*log(2*pi*par(8)^2)+norm([featX, featY])^2;
% 
% featX = (xdata - x_out); %calculate differences
% featY = (ydata - y_out); %calculate differences
% 
% llhX = - 0.5 * log( 2 * pi * par(7)^2) - 0.5 * (featX.^2/par(7)^2);
% llhY = - 0.5 * log( 2 * pi * par(8)^2) - 0.5 * (featY.^2/par(8)^2);
% test=-2*(sum([llhX,llhY]));
% 
% end



% function [sensitivity] = diff_log_summed(par,xdata,ydata,tspanX, tspanY, weight)
% %Calculates the differences between the data and simulation output
% %can be used as input for lsqnonlin/least squares optimizer
% %outputs:                          
% %   sensitivity = output feature, differences of both outputs
% 
% %inputs:
% %   tspan       timespan for integral
% %   par =       adapted input parameters
% %   data =      data vectors for x and y respectively
% 
% % determine output for adapted parameters
% [x_out, y_out] = intLea(par, tspanX, tspanY);
% x_out_logged=log(x_out);
% y_out_logged=log(y_out);
% xdata_logged=log(xdata);
% ydata_logged=log(ydata);
% %determine value of output of interest
% %featX = zeros(length(tspanX),1);
% %featY = zeros(length(tspanY),1);  
% 
% featX = (xdata_logged - x_out_logged);
% %calculate differences and normalize
% featY = sqrt((1/weight))*(ydata_logged - y_out_logged); %calculate differences
% featX(1)=sqrt((1/weight))*featX(1);
% 
% sensitivity = norm([featX, featY])^2;
% 
% end

% function [sensitivity] = diff_log_summed_sd(par,xdata,ydata,tspanX, tspanY, weight)
% %Calculates the differences between the data and simulation output
% %can be used as input for lsqnonlin/least squares optimizer
% %outputs:                          
% %   sensitivity = output feature, differences of both outputs
% 
% %inputs:
% %   tspan       timespan for integral
% %   par =       adapted input parameters
% %   data =      data vectors for x and y respectively
% 
% % determine output for adapted parameters
% [x_out, y_out] = intLea(par(1:6), tspanX, tspanY);
% x_out_logged=log(x_out);
% y_out_logged=log(y_out);
% xdata_logged=log(xdata);
% ydata_logged=log(ydata);
% %determine value of output of interest
% %featX = zeros(length(tspanX),1);
% %featY = zeros(length(tspanY),1);  
% 
% featX = (xdata_logged - x_out_logged)./(par(7));
% %calculate differences and normalize
% featY = sqrt((1/weight))*(ydata_logged - y_out_logged)./(par(8)); %calculate differences
% featX(1)=sqrt((1/weight))*featX(1);
% 
% %sensitivity = norm([featX, featY])^2+sum(log((2*pi*par(7)^2).*(xdata.^2)))+sum(log((2*pi*par(8)^2).*(ydata.^2)));
% sensitivity = norm([featX, featY])^2+length(xdata)*log(2*pi*par(7)^2)+length(ydata)*log(2*pi*par(8)^2);
% 
% end

% function [sensitivity] = diff_log_summed_sd_hier(par,xdata,ydata,tspanX, tspanY, weight)
% %Calculates the differences between the data and simulation output
% %can be used as input for lsqnonlin/least squares optimizer
% %outputs:                          
% %   sensitivity = output feature, differences of both outputs
% 
% %inputs:
% %   tspan       timespan for integral
% %   par =       adapted input parameters
% %   data =      data vectors for x and y respectively
% 
% % determine output for adapted parameters
% [x_out, y_out] = intLea(par(1:6), tspanX, tspanY);
% x_out_logged=log(x_out);
% y_out_logged=log(y_out);
% xdata_logged=log(xdata);
% ydata_logged=log(ydata);
% %determine value of output of interest
% %featX = zeros(length(tspanX),1);
% %featY = zeros(length(tspanY),1);  
% 
% 
% 
% featX = (xdata_logged - x_out_logged);
% featY = sqrt((1/weight))*(ydata_logged - y_out_logged); %calculate differences
% 
% sd_x=(1/length(xdata))*(norm(featX)^2);
% sd_y=(1/length(ydata))*(norm(featY)^2);
% %calculate differences and normalize
% featX(1)=sqrt((1/weight))*featX(1);
% sd_x=max(sd_x, 1e-7);
% featX=featX./sqrt(sd_x);
% featY=featY./sqrt(sd_y);
% 
% %sensitivity = norm([featX, featY])^2+sum(log((2*pi*par(7)^2).*(xdata.^2)))+sum(log((2*pi*par(8)^2).*(ydata.^2)));
% sensitivity = norm([featX, featY])^2+length(xdata)*log(2*pi*sd_x)+length(ydata)*log(2*pi*sd_y);
% 
% end

% function [sd_x, sd_y] = compute_mult_sd(par,xdata,ydata,tspanX, tspanY)
% [x_out, y_out] = intLea(par(1:6), tspanX, tspanY);
% %calculates the opt sd from the set of params and available data
% x_out_logged=log(x_out);
% y_out_logged=log(y_out);
% xdata_logged=log(xdata);
% ydata_logged=log(ydata);
% %determine value of output of interest
% %featX = zeros(length(tspanX),1);
% %featY = zeros(length(tspanY),1);  
% 
% 
% 
% featX = (xdata_logged - x_out_logged);
% featY = (ydata_logged - y_out_logged); %calculate differences
% 
% sd_x=sqrt((1/length(xdata))*(norm(featX)^2));
% sd_y=sqrt((1/length(ydata))*(norm(featY)^2));
% end

% function [sensitivity] = diff_log_norm_summed(par,xdata,ydata,tspanX, tspanY, weight)
% %Calculates the differences between the data and simulation output
% %can be used as input for lsqnonlin/least squares optimizer
% %outputs:                          
% %   sensitivity = output feature, differences of both outputs
% 
% %inputs:
% %   tspan       timespan for integral
% %   par =       adapted input parameters
% %   data =      data vectors for x and y respectively
% 
% % determine output for adapted parameters
% [x_out, y_out] = intLea(par, tspanX, tspanY);
% x_out_logged=log(x_out);
% y_out_logged=log(y_out);
% xdata_logged=log(xdata);
% ydata_logged=log(ydata);
% %determine value of output of interest
% %featX = zeros(length(tspanX),1);
% %featY = zeros(length(tspanY),1);  
% 
% featX = (xdata_logged - x_out_logged)./norm(xdata_logged);
% %calculate differences and normalize
% featY = sqrt((1/weight))*(ydata_logged - y_out_logged)./norm(ydata_logged); %calculate differences
% featX(1)=sqrt((1/weight))*featX(1);
% 
% sensitivity = norm([featX, featY])^2;
% 
% end



% function [parval, parres] = Profiles_lsq(lsqfun, par, num, lb, ub)
% %calculate the likelihood for a certain parameter at num values between 
% %the lb and the ub.
% 
% %Input:
% %   llhfun      function handle of the likelihood function
% %   par         optimized parameters
%                                 %  , parnum parnum      parameter to optimize
% %   num         number of points to calculate in the profile
% %   lb          lb of the problem
% %   ub          ub of the problem
% 
% %Output:
% %   parval      values of the parameter
% %   parres      llh value for the parameter
%     
% % can also look into MATLAB profile function
%     originalpar = par; %save original parameters
% 
%     parnum = length(par);
%     
%     parres = zeros(parnum,num) ;
%     parval = zeros(parnum,num) ;
% 
%     for j = 1:parnum
%         
%         parval(j,:) = linspace(lb(j), ub(j),num);
%         
%         for i = 1:num
%             par(j) = parval(j,i);
%             
%             parres(j,i) = lsqfun(par);
%         end
%         
%         par = originalpar; %reset to original input parameters
%         
%     end 
% 
% end