%% melanie_main_Burger_v4
% can be used to fit the Burger data
%-> use lsq with normalization
close all; clearvars;

%% Choose parameter

%patients you want to fit; 12 will be excluded as ibrutinib dose was paused
%(available patients: 1-30)
patients=9;
runns=1:3;
%if to use the constraint d2<d1
use_constraint=false;

%assumed volume of a CLL cell (in fl)
cll_volume=166;

%tissue measurements to be taken into account (max 3 available, Burger
%takes only one)
tss_measurements=1;

%if one should normaize when using the least square optimizer
normalize_lls=1;

%% Read the PB data
[num, txt, ~] = xlsread('Heavywater-all counts-decoded.xlsx');
[chars, charDesc,~]=xlsread('Heavy Water patients characterists.xlsx');
%% Prepare Variables for storing the opt. results
Xmulti=[];
Errormulti=[];


%% Perform fitting and plot the results
for patient=patients
    for runn=runns
%for patient=first_pat:last_pat
        if patient==12
            'Skipping patient_nr 12 as therapy was paused in between';
        else
            'Fitting Patient_nr '+ string(patient);
            close all; clearvars -except num txt chars charDesc patient use_constraint patients cll_volume Xmulti Xllh Errormulti Errorllh runn runns tss_measurements x normalize_lls tss_measurements
            %% optimization using integral

            %settings

            lb = [1e-3 1e-3 1e-6 0.7 0.9 0]; 
            ub = [5 20 0.1 1.3 1.1 1e12]; 



            n_starts = 1500;


            %construct start (later add multi start)
            par_init = [0.1 0.1 0.000 0 0 1e10]; %start parameters


            %% Set data

            %use available mess points as time span
            %use available mess values as true data

            %Burger ACC1
            [startDate,tspanY, ydata]=getPatientData(num, txt, chars, charDesc, patient);
            tspanY=tspanY';
            ydata=ydata';
            ty_outlier=[];
            y_outlier=[];
            if patient == 3
                tspanY=[tspanY(1:3),tspanY(5:end)];
                ydata=[ydata(1:3),ydata(5:end)];
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
            

            [tspanXavailable,Xdata]=getTissueDataComplex(patient,cll_volume);
            tspanXavailable=(tspanXavailable-startDate)';
            tspanXavailable(tspanXavailable<0)=0;
            tspanX=tspanXavailable(1:tss_measurements);
            xdata=Xdata(1:tss_measurements)';
            
            %xdata=Xdata(2);

            %use first measured values as initial parameter guess for x0, y0
            %allow them to vary +-10%

            lb(4) = lb(4) * Xdata(1) ; ub(4) = ub(4) * Xdata(1); par_init(4) = Xdata(1);
            %lb(4) = lb(4) * X0 ; ub(4) = ub(4) * X0; par_init(4) = X0;
          
            lb(5) = lb(5) * ydata(1) ; ub(5) = ub(5) * ydata(1); par_init(5) = ydata(1);
            %use the last tissue datapoint as upper bound for c
            ub(6)=Xdata(end);
            %% Get latin hypercube sampled starts

            [par_starts,~] = lhsdesign_modified(n_starts,lb,ub);
            startpoints = CustomStartPointSet(par_starts(:,1:6));

            %% optimization for integral

            %% define objective function.
            %if normalize_lls=true, use the appropriate objective function
            if normalize_lls==true
                lsqfun = @(par)diff_integral_norm_summed(par,xdata,ydata,tspanX, tspanY);
            else
                lsqfun = @(par)diff_integral_summed(par,xdata,ydata,tspanX, tspanY);
            end 
            
            %% create OptimProblem
            %if use_constraint==true,
            %add the linear inequality condition d_1>d_2, that is d2-d1<0
            %that is par(1)-par(2)<0
 
            options = optimoptions('fmincon');
            if use_constraint==true
                problem = createOptimProblem('fmincon','x0',par_init,'objective',lsqfun,...
                'lb',lb,'ub',ub,'xdata',xdata,'ydata',ydata,...
                'Aineq',[1 -1 0 0 0 0], 'bineq', 0, 'options', options);
            else
                problem = createOptimProblem('fmincon','x0',par_init,'objective',lsqfun,...
                'lb',lb,'ub',ub,'xdata',xdata,'ydata',ydata);
            end  

            %% run the optimization problem 
            ms = MultiStart('PlotFcns',{@gsplotbestf, @errorlandscape},'FunctionTolerance',1e-8);
            [xmulti,errormulti,exitflag,output,solutions] = run(ms,problem,startpoints)
            
            %% save params and values
            Xmulti=[Xmulti;[patient,xmulti]];
            Errormulti=[Errormulti;[patient,errormulti]];
           

            %% plot data

            tspan = 0:1:600;
            tspanx=0:1:500;
            [x_opt, y_opt] = intLea(xmulti, tspanx, tspan);

            figure;

            subplot(2,1,1)

            %plot lsq fitting
            %tissue measurements taken into account
            plot(tspanXavailable(1:tss_measurements), Xdata(1:tss_measurements), '*k'); hold on;
            %tissue measurements not taken into account
            plot(tspanXavailable(tss_measurements+1:end), Xdata(tss_measurements+1:end), 'og'); hold on;
            %plot(tspan, x_true, 'g'); hold on;
            plot(tspanx, x_opt, '--m',xlabel('days'),ylabel('absolute lymphocyte count')); 
            legend('measurements used for fitting','additional measurements','solution lsq');

            plotTitle=title('tissue_ACC_'+string(patient)+'_'+string(runn),'Interpreter','none');
            %plotTitle.Color='red';
            plotTitle.FontSize=16;

            %save figure as PDF
            f=gcf;
            exportgraphics(f,[pwd,strcat('/Figures/',plotTitle.String,'.png')],'Resolution',300)

            figure;
            subplot(2,1,1)
            plot(tspanY, ydata, '*k'); hold on;
            %outliner not taken into account
           

            %plot(tspan, y_true, 'color',[0 0.5 0]); hold on;
            plot(tspan, y_opt, '--m',xlabel('days'),ylabel('absolute lymphocyte count')); hold on;
            plot(ty_outlier, y_outlier, 'og'); 
            legend('measurements','solution lsq')

            plotTitle=title('PB_ACC_'+string(patient)+'_'+string(runn),'Interpreter','none');
            %plotTitle.Color='red';
            plotTitle.FontSize=16;


            %save figure as PDF
            f=gcf;
            exportgraphics(f,[pwd,strcat('/Figures/',plotTitle.String,'.png')],'Resolution',300)
            %% Profile Likelihood Analysis


            [parval, parres] = Profiles_lsq(lsqfun, xmulti, 100, lb, ub);
            figure;
            plotTitle=sgtitle('parval_ls_ACC_'+string(patient)+'_'+string(runn),'Interpreter','none');
            %plotTitle.Color='red';
            plotTitle.FontSize=16;
            subplot(2,3,1)
            plot(parval(1,:),parres(1,:));xlabel('parameter value');ylabel('sum of squared differences');title('d2');
            subplot(2,3,2)
            plot(parval(2,:),parres(2,:));xlabel('parameter value');ylabel('sum of squared differences');title('d1');
            subplot(2,3,3)
            plot(parval(3,:),parres(3,:));xlabel('parameter value');ylabel('sum of squared differences');title('m');
            subplot(2,3,4)
            plot(parval(4,:),parres(4,:));xlabel('parameter value');ylabel('sum of squared differences');title('x0');
            subplot(2,3,5)
            plot(parval(5,:),parres(5,:));xlabel('parameter value');ylabel('sum of squared differences');title('y0');
            subplot(2,3,6)
            plot(parval(6,:),parres(6,:));xlabel('parameter value');ylabel('sum of squared differences');title('c');

            %save figure as PDF
            f=gcf;
            exportgraphics(f,[pwd,strcat('/Figures/',plotTitle.String,'.png')],'Resolution',300)
        end 
    end
end

%% Store Parameters and errors
writetable(array2table(Xmulti, 'VariableNames',{'ACC','d2', 'd1','m','x0','y0','c'}),'Params_lss.xlsx');

writetable(array2table(Errormulti, 'VariableNames',{'ACC','ss'}),'Errors_ss.xlsx');

%% %%% Function to determine features %%%%%

function [spac_days] = setSpacing(days, n_samples, distribution)
%Input:
%   days        amount of days in total
%   n_samples   amount of samples
%   distribution distritbution of the days, linear, log or log10

%Output:
%   spac_days   spacing of days

    if distribution == 1 %'linear'
        spacing = [0, linspace(0, days, n_samples-1)];
    elseif distribution == 2 % 'log10'
        spacing = [0, logspace(0, log10(days), n_samples-1)];
    elseif distribution == 3 %'log'
        spacing = [0, exp(linspace(1,log(days),n_samples-1))];
        
    end
    
    spac_days = round(spacing);
    
    for i = 1:(length(spac_days)-1)
        if spac_days(i) == spac_days(i+1)
            spac_days(i+1) = spac_days(i+1) + 1;
        end
    end
      
end 


function [sensitivity] = diff_integral(par,xdata,ydata,tspanX, tspanY)
%Calculates the differences between the data and simulation output,
%normalized via the l2_norm of the measurements
%can be used as input for lsqnonlin/least squares optimizer
%outputs:                          
%   sensitivity = output feature, differences of both outputs

%inputs:
%   tspan       timespan for integral
%   par =       adapted input parameters
%   data =      data vectors for x and y respectively

% determine output for adapted parameters
[x_out, y_out] = intLea(par, tspanX, tspanY);

%determine value of output of interest
featX = zeros(length(tspanX),1);
featY = zeros(length(tspanY),1);  

featX = (xdata - x_out); %calculate differences
featY = (ydata - y_out); %calculate differences

sensitivity = [featX, featY];

end 

function [sensitivity] = diff_integral_norm(par,xdata,ydata,tspanX, tspanY)
%Calculates the differences between the data and simulation output
%can be used as input for lsqnonlin/least squares optimizer
%outputs:                          
%   sensitivity = output feature, differences of both outputs

%inputs:
%   tspan       timespan for integral
%   par =       adapted input parameters
%   data =      data vectors for x and y respectively

% determine output for adapted parameters
[x_out, y_out] = intLea(par, tspanX, tspanY);

%determine value of output of interest
%featX = zeros(length(tspanX),1);
%featY = zeros(length(tspanY),1);  

featX = (xdata - x_out)./norm(xdata); %calculate differences and normalize
featY = (ydata - y_out)./norm(ydata); %calculate differences

sensitivity = [featX, featY];

end

function [sensitivity] = diff_integral_summed(par,xdata,ydata,tspanX, tspanY)
%Calculates the differences between the data and simulation output,
%normalized via the l2_norm of the measurements
%can be used as input for lsqnonlin/least squares optimizer
%outputs:                          
%   sensitivity = output feature, differences of both outputs

%inputs:
%   tspan       timespan for integral
%   par =       adapted input parameters
%   data =      data vectors for x and y respectively

% determine output for adapted parameters
[x_out, y_out] = intLea(par, tspanX, tspanY);

%determine value of output of interest
featX = zeros(length(tspanX),1);
featY = zeros(length(tspanY),1);  

featX = (xdata - x_out); %calculate differences
featY = (ydata - y_out); %calculate differences

sensitivity = norm([featX, featY])^2;

end 

function [sensitivity] = diff_integral_norm_summed(par,xdata,ydata,tspanX, tspanY)
%Calculates the differences between the data and simulation output
%can be used as input for lsqnonlin/least squares optimizer
%outputs:                          
%   sensitivity = output feature, differences of both outputs

%inputs:
%   tspan       timespan for integral
%   par =       adapted input parameters
%   data =      data vectors for x and y respectively

% determine output for adapted parameters
[x_out, y_out] = intLea(par, tspanX, tspanY);

%determine value of output of interest
%featX = zeros(length(tspanX),1);
%featY = zeros(length(tspanY),1);  

featX = (xdata - x_out)./norm(xdata); %calculate differences and normalize
featY = (ydata - y_out)./norm(ydata); %calculate differences

sensitivity = norm([featX, featY])^2;

end

function [negllh] = negLLH(par, xdata, ydata, tspanX, tspanY, optimize_sd)
%calculates the negative log likelihood for minimization
%outputs:                          
%   sensitivity = output feature, squared difference of both outputs

%inputs:
%   tspan       timespan for integral
%   par =       adapted input parameters, including standard deviation
%   data =      data vectors for x and y respectively

% determine output for adapted parameters
[x_out, y_out] = intLea(par(1:6), tspanX, tspanY);

if optimize_sd
    sigmX = par(7); 
    sigmY = par(8);
else
    sigmX=1;
    sigmY=1;
end

%determine value of output of interest
  
featX = (xdata - x_out); %calculate differences
featY = (ydata - y_out); %calculate differences

llhX = - 0.5 * log( 2 * pi * sigmX^2) - 0.5 * (featX.^2/sigmX^2);
llhY = - 0.5 * log( 2 * pi * sigmY^2) - 0.5 * (featY.^2/sigmY^2);
%http://jrmeyer.github.io/machinelearning/2017/08/18/mle.html

negllh = -1*(sum(llhX) + sum(llhY)) ; % possibility to apply weight factors

end 

function [negllh] = negLLH_sd(par, xdata, ydata, tspanX, tspanY, xmultillh)
%calculates the negative log likelihood for minimization
%outputs:                          
%   sensitivity = output feature, squared difference of both outputs

%inputs:
%   tspan       timespan for integral
%   par =       adapted input parameters, including standard deviation
%   data =      data vectors for x and y respectively

% determine output for adapted parameters
[x_out, y_out] = intLea(xmultillh(1:6), tspanX, tspanY);


sigmX = par; 
sigmY = par;


%determine value of output of interest
  
featX = (xdata - x_out); %calculate differences
featY = (ydata - y_out); %calculate differences

llhX = - 0.5 * log( 2 * pi * sigmX^2) - 0.5 * (featX.^2/sigmX^2);
llhY = - 0.5 * log( 2 * pi * sigmY^2) - 0.5 * (featY.^2/sigmY^2);
%http://jrmeyer.github.io/machinelearning/2017/08/18/mle.html

negllh = -1*(sum(llhX) + sum(llhY)) ; % possibility to apply weight factors

end 

function [negllh] = negLLH_norm(par, xdata, ydata, tspanX, tspanY, optimize_sd)
%calculates the negative log likelihood for minimization, taking into
%account normalization
%outputs:                          
%   sensitivity = output feature, squared difference of both outputs

%inputs:
%   tspan       timespan for integral
%   par =       adapted input parameters, including standard deviation
%   data =      data vectors for x and y respectively

% determine output for adapted parameters
[x_out, y_out] = intLea(par(1:6), tspanX, tspanY);

if optimize_sd==true
    sigmX = par(7); 
    sigmY = par(8);
else
    sigmX=1;
    sigmY=1;
end

%determine value of output of interest


featX = (xdata - x_out)./norm(ydata); %calculate differences
featY = (ydata - y_out)./norm(xdata); %calculate differences

llhX = - 0.5 * log( 2 * pi * sigmX^2) - 0.5 * (featX.^2/sigmX^2);
llhY = - 0.5 * log( 2 * pi * sigmY^2) - 0.5 * (featY.^2/sigmY^2);
%http://jrmeyer.github.io/machinelearning/2017/08/18/mle.html

negllh = -1*(sum(llhX) + sum(llhY)) ; % possibility to apply weight factors

end

function [x_out, y_out] = intLea(par, tspanX, tspanY)
% Integraded function, integration established by Lea
%Outputs:
%   time            time span used in the integral
%   x/y _out        output from integral for respectively x and y

%Input:
%   par             parameters for integral
%   tspan X/Y       time span for respectively x and y

d2 = par(1);
d1 = par(2);
m = par(3);
x0 = par(4);
y0 = par(5);
c = par(6);

A = d1*c;
B = (m+d1);

x_out = A/B + (x0 - A/B)*exp(-B*tspanX);

y_out = (m*A)/(B*d2) + ...
   ( m/(d2 - B) * (x0 - A/B)*exp(-B*tspanY) ) + ...
    (y0 - (m*A)/(B*d2) - (m/(d2 - B) * (x0 - A/B)))*exp(-d2 * tspanY);

end


function [parval, parres] = Profiles_llh(llhfun, par, num, lb, ub)
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
    originalpar = par %save original parameters

    parnum = length(par);
    
    parres = zeros(parnum,num) ;
    parval = zeros(parnum,num) ;

    for j = 1:parnum
        
        parval(j,:) = linspace(lb(j), ub(j),num);
        
        for i = 1:num
            par(j) = parval(j,i);
            
            parres(j,i) = llhfun(par);
        end
        
        par = originalpar %reset to original input parameters
        
    end 
    
end


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
            
            parres(j,i) = sum(lsqfun(par).^2);
        end
        
        par = originalpar; %reset to original input parameters
        
    end 

end