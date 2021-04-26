%% melanie_main_Burger_sls
% can be used to fit the Burger data
%-> use lsq with normalization

%only intermediate, to be ablbe to run twice a day
close all; clearvars;

%% Choose parameter

%patients you want to fit; 12 will be excluded as ibrutinib dose was paused
%(available patients: 1-30)
patients=[1:11,13:30];
runns=1:3;
%if to use the constraint d2<d1
use_constraint=false;

%assumed volume of a CLL cell (in fl)
cll_volume=166;

%tissue measurements to be taken into account (max 3 available, Burger
%takes only one)
tss_measurements=1;

%if one should use the complex computation of the tissue cll burden
use_complex_tss=1;

%if one should normaize when using the least square optimizer
normalize_lls=0;

%weight for the second tissue measurement (if second measurement is taken
%into account)
weight=1;

%assume multiplicative noise
mult_noise=1;

%if optimizing sd, when assuming multiplicative noise
opt_sd=1;


%folder to save figures
folder='/Figures/08_20_2Comp_1tissue/';
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

           if patient == 13
             lb = [1e-5 1e-5 1e-6 0.9 0.9 0]; 
             ub = [2 0.5 1 1.1 1.1 1e12];
           else  
             lb = [1e-5 1e-7 1e-7 0.7 0.7 1e9]; 
             ub = [5 20 1 1.5 1.5 1e12]; 
           end 
%             end
            lb_sd = [1e-7 1e-7];
            ub_sd = [1 1];
            
            n_starts = 2500;


            %construct start (later add multi start)
            par_init = [0.1 0.1 0.000 0 0 1e10]; %start parameters
            sdXY=[0.1 0.1];

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
            %% account for sd
            if opt_sd==true && mult_noise==0
                lb=[lb, lb_sd];
                ub=[ub, ub_sd];
                par_init=[par_init,sdXY];
            end
                
            %% Get latin hypercube sampled starts

            [par_starts,~] = lhsdesign_modified(n_starts,lb,ub);

            startpoints = CustomStartPointSet(par_starts(:,:));

            %% optimization for integral

            %% define objective function.
            %if normalize_lls=true, use the appropriate objective function
            if normalize_lls==true
                if mult_noise==true
                        lsqfun=@(par)diff_log_norm_summed(par,xdata,ydata,tspanX, tspanY, weight);
                else
                    lsqfun = @(par)diff_integral_norm_summed(par,xdata,ydata,tspanX, tspanY, weight);
                end
            else
                if mult_noise==true
                    if opt_sd==false
                        lsqfun=@(par)diff_log_summed(par,xdata,ydata,tspanX, tspanY, weight);
                    else
                        lsqfun=@(par)diff_log_summed_sd_hier(par,xdata,ydata,tspanX, tspanY, weight);
                    end
                else
                    if opt_sd==false
                        lsqfun = @(par)diff_integral_summed(par,xdata,ydata,tspanX, tspanY, weight);
                    else
                        lsqfun = @(par)diff_integral_summed_sd(par,xdata,ydata,tspanX, tspanY, weight);
                    end
                end 
            end
            
            %% create OptimProblem
            %if use_constraint==true,
            %add the linear inequality condition d_1>d_2, that is d2-d1<0
            %that is par(1)-par(2)<0
 
            options = optimoptions('fmincon');
            if use_constraint==true
                problem = createOptimProblem('fmincon','x0',par_init,'objective',lsqfun,...
                'lb',lb,'ub',ub,'xdata',xdatas,'ydata',ydata,...
                'Aineq',[1 -1 0 0 0 0], 'bineq', 0, 'options', options);
            else
                problem = createOptimProblem('fmincon','x0',par_init,'objective',lsqfun,...
                'lb',lb,'ub',ub,'xdata',xdata);
            end  

            %% run the optimization problem 
            ms = MultiStart('PlotFcns',{@gsplotbestf, @errorlandscape},'FunctionTolerance',1e-8);
            [xmulti,errormulti,exitflag,output,solutions] = run(ms,problem,startpoints);
            
            if mult_noise==1 && opt_sd==1
                [sd_x, sd_y]=compute_mult_sd(xmulti,xdata,ydata,tspanX, tspanY);
                xmulti=[xmulti, sd_x, sd_y];
            end
            
            %% save params and values
            Xmulti=[Xmulti;[patient,xmulti]];
            Errormulti=[Errormulti;[patient,errormulti]];
           

            %% plot data

            tspan = 0:1:600;
            tspanx=0:1:500;
            [x_opt, y_opt] = intLea(xmulti(1:6), tspanx, tspan);

            figure;

            subplot(2,1,1)

            %plot lsq fitting
            %tissue measurements taken into account
            plot(tspanXavailable(1:tss_measurements), Xdata(1:tss_measurements), '*k'); hold on;
            %tissue measurements not taken into account
            if length(tspanXavailable)>tss_measurements
                plot(tspanXavailable(tss_measurements+1:end), Xdata(tss_measurements+1:end), 'og'); hold on;
            end
                %plot(tspan, x_true, 'g'); hold on;
            plot(tspanx, x_opt, '--m',xlabel('days'),ylabel('absolute lymphocyte count')); 
            legend('measurements used for fitting','additional measurements','solution lsq');

            plotTitle=title('tissue_ACC_'+string(patient)+'_'+string(runn),'Interpreter','none');
            %plotTitle.Color='red';
            plotTitle.FontSize=16;

            %save figure as PDF
            f=gcf;
            exportgraphics(f,[pwd,strcat(folder,plotTitle.String,'.png')],'Resolution',300)

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
            exportgraphics(f,[pwd,strcat(folder,plotTitle.String,'.png')],'Resolution',300)
            %% Profile Likelihood Analysis

            if mult_noise==1 && opt_sd==1
                original_lsqfun=@(par)diff_log_summed_sd(par,xdata,ydata,tspanX, tspanY, weight);
                [parval, parres] = Profiles_lsq(original_lsqfun, xmulti, 10*n_starts, [lb,1e-8,1e-8], [ub,0.2,1]);
            else
                [parval, parres] = Profiles_lsq(lsqfun, xmulti, 10*n_starts, lb, ub);
            end
            figure;
            plotTitle=sgtitle('parval_ls_ACC_'+string(patient)+'_'+string(runn),'Interpreter','none');
            %plotTitle.Color='red';
            plotTitle.FontSize=16;
            subplot(2,4,1)
            plot(parval(1,:),parres(1,:));xlabel('parameter value');ylabel('sum of squared differences');title('d2');
            subplot(2,4,2)
            plot(parval(2,:),parres(2,:));xlabel('parameter value');ylabel('sum of squared differences');title('d1');
            subplot(2,4,3)
            plot(parval(3,:),parres(3,:));xlabel('parameter value');ylabel('sum of squared differences');title('m');
            subplot(2,4,4)
            plot(parval(4,:),parres(4,:));xlabel('parameter value');ylabel('sum of squared differences');title('x0');
            subplot(2,4,5)
            plot(parval(5,:),parres(5,:));xlabel('parameter value');ylabel('sum of squared differences');title('y0');
            subplot(2,4,6)
            plot(parval(6,:),parres(6,:));xlabel('parameter value');ylabel('sum of squared differences');title('c');
            if opt_sd==true
               subplot(2,4,7)
               plot(parval(7,:),parres(7,:));xlabel('parameter value');ylabel('sum of squared differences');title('sdx');
               subplot(2,4,8)
               plot(parval(8,:),parres(8,:));xlabel('parameter value');ylabel('sum of squared differences');title('sdy');
            end
            %save figure as PDF
            f=gcf;
            exportgraphics(f,[pwd,strcat(folder,plotTitle.String,'.png')],'Resolution',300);
            savefig(f,[pwd,strcat(folder,plotTitle.String)],'compact')
        end 
    end
end

%% Store Parameters and errors
if opt_sd==false
    writetable(array2table(Xmulti, 'VariableNames',{'ACC','d2', 'd1','m','x0','y0','c'}),'Params_lss.xlsx');
else
    writetable(array2table(Xmulti, 'VariableNames',{'ACC','d2', 'd1','m','x0','y0','c','sdx','sdy'}),'Params_lss.xlsx');
end
writetable(array2table(Errormulti, 'VariableNames',{'ACC','ss'}),'Errors_ss.xlsx');

%% Functions

function [sensitivity] = diff_integral_summed(par,xdata,ydata,tspanX, tspanY, weight)
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
 
featX = sqrt((1/weight))*(xdata - x_out); %calculate differences
featY = (ydata - y_out); %calculate differences
featY(1)=sqrt((1/weight))*featY(1);
sensitivity = norm([featX, featY])^2;

end 

function [sensitivity] = diff_integral_summed_sd(par,xdata,ydata,tspanX, tspanY, weight)
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
[x_out, y_out] = intLea(par(1:6), tspanX, tspanY);
 
%featX = sqrt((1/weight))*(xdata - x_out)/par(7); %calculate differences
featX = (xdata - x_out)/par(7);
featY = (ydata - y_out)/par(8); %calculate differences
%featY(1)=sqrt((1/weight))*featY(1);
sensitivity = length(xdata)*log(2*pi*par(7)^2)+length(ydata)*log(2*pi*par(8)^2)+norm([featX, featY])^2;

featX = (xdata - x_out); %calculate differences
featY = (ydata - y_out); %calculate differences

llhX = - 0.5 * log( 2 * pi * par(7)^2) - 0.5 * (featX.^2/par(7)^2);
llhY = - 0.5 * log( 2 * pi * par(8)^2) - 0.5 * (featY.^2/par(8)^2);
test=-2*(sum([llhX,llhY]));

end

function [sensitivity] = diff_integral_norm_summed(par,xdata,ydata,tspanX, tspanY, weight)
%Calculates the differences between the data and simulation output
%can be used as input for lsqnonlin/least squares optimizer
%outputs:                          
%   sensitivity = output feature, differences of both outputs

%inputs:
%   tspan       timespan for integral
%   par =       adapted input parameters
%   data =      data vectors for x and y respectively

% determine output for adapted parameters
[x_out, y_out] = intLea(par(1:6), tspanX, tspanY);

%determine value of output of interest
%featX = zeros(length(tspanX),1);
%featY = zeros(length(tspanY),1);  

featX = (xdata - x_out)./norm(xdata);
%calculate differences and normalize
featY = sqrt((1/weight))*(ydata - y_out)./norm(ydata); %calculate differences
%featX(1)=sqrt((1/weight))*featX(1);

sensitivity = norm([featX, featY])^2;

end

function [sensitivity] = diff_log_summed(par,xdata,ydata,tspanX, tspanY, weight)
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
x_out_logged=log(x_out);
y_out_logged=log(y_out);
xdata_logged=log(xdata);
ydata_logged=log(ydata);
%determine value of output of interest
%featX = zeros(length(tspanX),1);
%featY = zeros(length(tspanY),1);  

featX = (xdata_logged - x_out_logged);
%calculate differences and normalize
featY = sqrt((1/weight))*(ydata_logged - y_out_logged); %calculate differences
featX(1)=sqrt((1/weight))*featX(1);

sensitivity = norm([featX, featY])^2;

end

function [sensitivity] = diff_log_summed_sd(par,xdata,ydata,tspanX, tspanY, weight)
%Calculates the differences between the data and simulation output
%can be used as input for lsqnonlin/least squares optimizer
%outputs:                          
%   sensitivity = output feature, differences of both outputs

%inputs:
%   tspan       timespan for integral
%   par =       adapted input parameters
%   data =      data vectors for x and y respectively

% determine output for adapted parameters
[x_out, y_out] = intLea(par(1:6), tspanX, tspanY);
x_out_logged=log(x_out);
y_out_logged=log(y_out);
xdata_logged=log(xdata);
ydata_logged=log(ydata);
%determine value of output of interest
%featX = zeros(length(tspanX),1);
%featY = zeros(length(tspanY),1);  

featX = (xdata_logged - x_out_logged)./(par(7));
%calculate differences and normalize
featY = sqrt((1/weight))*(ydata_logged - y_out_logged)./(par(8)); %calculate differences
featX(1)=sqrt((1/weight))*featX(1);

%sensitivity = norm([featX, featY])^2+sum(log((2*pi*par(7)^2).*(xdata.^2)))+sum(log((2*pi*par(8)^2).*(ydata.^2)));
sensitivity = norm([featX, featY])^2+length(xdata)*log(2*pi*par(7)^2)+length(ydata)*log(2*pi*par(8)^2);

end

function [sensitivity] = diff_log_summed_sd_hier(par,xdata,ydata,tspanX, tspanY, weight)
%Calculates the differences between the data and simulation output
%can be used as input for lsqnonlin/least squares optimizer
%outputs:                          
%   sensitivity = output feature, differences of both outputs

%inputs:
%   tspan       timespan for integral
%   par =       adapted input parameters
%   data =      data vectors for x and y respectively

% determine output for adapted parameters
[x_out, y_out] = intLea(par(1:6), tspanX, tspanY);
x_out_logged=log(x_out);
y_out_logged=log(y_out);
xdata_logged=log(xdata);
ydata_logged=log(ydata);
%determine value of output of interest
%featX = zeros(length(tspanX),1);
%featY = zeros(length(tspanY),1);  



featX = (xdata_logged - x_out_logged);
featY = sqrt((1/weight))*(ydata_logged - y_out_logged); %calculate differences

sd_x=(1/length(xdata))*(norm(featX)^2);
sd_y=(1/length(ydata))*(norm(featY)^2);
%calculate differences and normalize
featX(1)=sqrt((1/weight))*featX(1);
sd_x=max(sd_x, 1e-7);
featX=featX./sqrt(sd_x);
featY=featY./sqrt(sd_y);

%sensitivity = norm([featX, featY])^2+sum(log((2*pi*par(7)^2).*(xdata.^2)))+sum(log((2*pi*par(8)^2).*(ydata.^2)));
sensitivity = norm([featX, featY])^2+length(xdata)*log(2*pi*sd_x)+length(ydata)*log(2*pi*sd_y);

end

function [sd_x, sd_y] = compute_mult_sd(par,xdata,ydata,tspanX, tspanY)
[x_out, y_out] = intLea(par(1:6), tspanX, tspanY);
%calculates the opt sd from the set of params and available data
x_out_logged=log(x_out);
y_out_logged=log(y_out);
xdata_logged=log(xdata);
ydata_logged=log(ydata);
%determine value of output of interest
%featX = zeros(length(tspanX),1);
%featY = zeros(length(tspanY),1);  



featX = (xdata_logged - x_out_logged);
featY = (ydata_logged - y_out_logged); %calculate differences

sd_x=sqrt((1/length(xdata))*(norm(featX)^2));
sd_y=sqrt((1/length(ydata))*(norm(featY)^2));
end

function [sensitivity] = diff_log_norm_summed(par,xdata,ydata,tspanX, tspanY, weight)
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
x_out_logged=log(x_out);
y_out_logged=log(y_out);
xdata_logged=log(xdata);
ydata_logged=log(ydata);
%determine value of output of interest
%featX = zeros(length(tspanX),1);
%featY = zeros(length(tspanY),1);  

featX = (xdata_logged - x_out_logged)./norm(xdata_logged);
%calculate differences and normalize
featY = sqrt((1/weight))*(ydata_logged - y_out_logged)./norm(ydata_logged); %calculate differences
featX(1)=sqrt((1/weight))*featX(1);

sensitivity = norm([featX, featY])^2;

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