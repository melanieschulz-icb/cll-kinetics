close all; clear all

%% optimization using integrals
%settings
par_slow = [0.002 0.027 0.0096 3034e9 153e9 0.4e9]; %wodarz 1
par_fast = [0.022 0.015 0.0177 3064e9 58e9  47.5e9]; %wodarz 2

lb = [1e-3 1e-3 1e-4 0.75 0.75 1e5 1e-5 1e-5]; 
ub = [0.1  0.1  0.5  1.25 1.25 1   1e10 1e10]; 


days =     [100, 100,  100, 70, 70, 360, 360,    100, 100, 100, 70, 70, 360, 360, 100] ;
samplesX = [1,   1,    1,   1,  1,  1,   1,      1,   1,   1,   1,  1,  1  , 1,  1] ;
samplesY = [8,   15,   8,   11, 11, 26,  17,     8,   15,  8,   11, 11, 26,  17,  8];

spacing =  [1,   1,   1,   1,  1,  1,  1,  1, 1,  1,  1, 1, 1, 1, 1]; % 1 = linear; 2 = log10; 3 = log
taily = {[130,160,190,220,250,280,310,340,370];[130,160,190,220,250,280,310,340,370]; [130,190,250,310,370];[100,130,160,190,220,250,280,310,340,370];[130,190,250,310,370];[];[]; ...
    [130,160,190,220,250,280,310,340,370,600];[130,160,190,220,250,280,310,340,370,600];[130,190,250,310,370,600];[100,130,160,190,220,250,280,310,340,370,600];[130,190,250,310,370,600];[600];[600];[130,160,190,220,250,280,310,340,370,400,430,460,490,520,550,580]};
tailx = {[0];[0];[0];[0];[0];[0];[0];[0];[0];[0];[0];[0];[0];[0];[0];[0];[0]};

n_starts = 100;

tests = 9; %times to do all runs, make sure it is uneven for medians

%construct start (later add multi start)
par_init = [0.001 0.001 0.0005 0 0 1e10 1e2 1e2]; %start parameters

plotlength = 1:700;

SSEframe = 1:500;

runs = length(samplesY);

%% Compare different amount of samples in Y

par_listslow = zeros(runs, tests, 8);
par_listfast = zeros(runs, tests, 8);

resnorm_list = [];

parval_l_slow = [];
parres_l_slow = [];
parval_l_fast = [];
parres_l_fast = [];

tspanX_all = [];
tspanY_all = [];

for i = 1: runs
    for q = 1:tests
        
        [tspanX] = [0,160,250];%setSpacing(days(i), samplesX(i), tailx{i}, spacing(i));
        [tspanY] = setSpacing(days(i), samplesY(i), taily{i}, spacing(i));
        
        [par_opt_slow, resnorm, exitflag,output,solutions, optfunsl,lbsl,ubsl] = ...
            llhFits(par_slow, lb, ub, tspanX, tspanY, n_starts, par_init);
        [par_opt_fast, resnorm, exitflag,output,solutions, optfunfs,lbfs,ubfs] = ...
            llhFits(par_fast, lb, ub, tspanX, tspanY, n_starts, par_init);
        
        [parval_sl, parres_sl] = Profiles_llh(log(par_slow), par_opt_slow, SSEframe, 100, log(lbsl), log(ubsl),optfunsl);
        [parval_fs, parres_fs] = Profiles_llh(log(par_fast), par_opt_fast, SSEframe, 100, log(lbfs), log(ubfs),optfunfs);
        
        parval_l_slow{i,q} = parval_sl; parval_l_fast{i,q} = parval_fs;
        parres_l_slow{i,q} = parres_sl; parres_l_fast{i,q} = parres_fs;
        
        par_listslow(i,q,:) = par_opt_slow;
        par_listfast(i,q,:) = par_opt_fast;
        
        tspanX_all{i,q} = tspanX;
        tspanY_all{i,q} = tspanY;
        
        clear parval_sl parval_fs parres_sl parres_fs lbsl lbfs ubsl ubfs
    
    end
end



%% Plot data distribution
[x_slow, y_slow] = intLea(par_slow, plotlength, plotlength);
[x_fast, y_fast] = intLea(par_fast, plotlength, plotlength);

figure; 
subplot(1,2,1); title('Y fast responder'); xlabel('days');title('Spacing');hold on

set(gca, 'YTick',[])
grid minor
subplot(1,2,2); 
set(gca, 'XTick', 0.75:0.5:runs+0.5,'YAxisLocation','right','XGrid','on','YMinorGrid','on'); hold on; 
set(gca, 'Xticklabel', {'slow','fast'},'fontsize',6);ylabel('normalized sum of squared errors','fontsize',10);
posmedianfs = zeros(1,runs); %store the positions of the medians for plots later
posmediansl = zeros(1,runs);

for p = 1:runs
    
    SSEfast = zeros(1,tests);
    SSEslow = zeros(1,tests);
    
    for s = 1: tests
    
        SSEslow(s) = sumofsquares(par_listslow(p,s,:),log(par_slow),SSEframe);
        
        SSEfast(s) = sumofsquares(par_listfast(p,s,:),log(par_fast),SSEframe);
        
    end
    
    meanfast = mean(SSEfast); meanslow = mean(SSEslow);
    stdfast = std(SSEfast);   stdslow = std(SSEslow);
    
    posmedianfs(p) = find(SSEfast == median(SSEfast));
    posmediansl(p) = find(SSEslow == median(SSEslow));
    
    subplot(1,2,1); %plot spacing
    plot([tspanY_all{p,s}, tspanX_all{p,s}], [-p*0.3*ones(1,length(tspanY_all{p,s}))+0.1, -p*0.3*ones(1,length(tspanX_all{p,s}))], '*','DisplayName',[string(p)]); hold on;
        
    subplot(1,2,2); %plot SSE values
    xline(p+0.5);
    errorbar([p-0.25 p+0.25], [meanslow meanfast], [stdslow stdfast],'*'); hold on

    clear SSEfast SSEslow meanfast meanslow stdfast stdslow
end
xline(7.5, '.-','LineWidth', 1.2); hold on

subplot(1,2,1);
plot(plotlength, y_slow/2e12,'color',[0.7 0.7 0.7],'LineWidth',1.2);hold on
plot(plotlength, x_slow/2e12,'color',[0.7 0.7 0.7],'LineWidth',1.2);hold on
plot(plotlength, y_fast/2e12,'color',[0.3 0.3 0.3],'LineWidth',1.2);hold on
plot(plotlength, x_fast/2e12,'color',[0.3 0.3 0.3],'LineWidth',1.2);hold on


%% plot resulting models

tspan = plotlength;

figure;subplot(2,2,1);
[x_true, y_true] = intLea(par_fast, tspan, tspan);
for j = 1:runs
    
    [x_opt, y_opt] = intLea(exp(par_listfast(j,posmedianfs(j),:)),tspan,tspan);
    plot(tspan, y_opt, '--','DisplayName',['y = '+string(samplesY(j))],'LineWidth',1.2); hold on;
    
end

plot(tspan, y_true, 'g','LineWidth',1.2);xlabel('days');ylabel('absolute lymphocyte counts');title('blood compartment fast');hold on
legend

subplot(2,2,2);
for k = 1:runs
    
    [x_opt, y_opt] = intLea(exp(par_listfast(k,posmedianfs(k),:)),tspan,tspan);
    plot(tspan, x_opt, '--','DisplayName',['y = '+string(samplesY(k))],'LineWidth',1.2); hold on;
       
end

plot(tspan, x_true, 'g','LineWidth',1.2);xlabel('days');ylabel('absolute lymphocyte counts');title('tissue compartment fast');hold on

subplot(2,2,3);
[x_true2, y_true2] = intLea(par_slow, tspan, tspan);
for j = 1:runs
    
    [x_opt, y_opt] = intLea(exp(par_listslow(j,posmediansl(j),:)),tspan,tspan);
    plot(tspan, y_opt, '--','DisplayName',['y = '+string(samplesY(j))],'LineWidth',1.2); hold on;
    
end

plot(tspan, y_true2, 'g','LineWidth',1.2);xlabel('days');ylabel('absolute lymphocyte counts');title('blood compartment slow');hold on

subplot(2,2,4);
for k = 1:runs
    
    [x_opt, y_opt] = intLea(exp(par_listslow(k,posmediansl(k),:)),tspan,tspan);
    plot(tspan, x_opt, '--','DisplayName',['y = '+string(samplesY(k))],'LineWidth',1.2); hold on;
       
end

plot(tspan, x_true2, 'g','LineWidth',1.2);xlabel('days');ylabel('absolute lymphocyte counts');title('tissue compartment slow');hold on


%% Plot resulting profiles
figure;
subplot(2,8,1);xline(log(par_slow(1)));xlabel('parameter values slow');ylabel('normalized sum of squares');title('d2');hold on
subplot(2,8,2);xline(log(par_slow(2)));xlabel('parameter value');title('d1');hold on
subplot(2,8,3);xline(log(par_slow(3)));xlabel('parameter value');title('m');hold on
subplot(2,8,4);xline(log(par_slow(4)));xlabel('parameter value');title('x0');hold on
subplot(2,8,5);xline(log(par_slow(5)));xlabel('parameter value');title('y0');hold on
subplot(2,8,6);xline(log(par_slow(6)));xlabel('parameter value');title('c');hold on
subplot(2,8,7);xlabel('parameter value');ylabel('negative log likelihood');title('sd x');hold on
subplot(2,8,8);xlabel('parameter value');title('sd y');hold on
subplot(2,8,9);xline(log(par_fast(1)));xlabel('parameter values fast');ylabel('normalized sum of squares');title('d2');hold on
subplot(2,8,10);xline(log(par_fast(2)));xlabel('parameter value');title('d1');hold on
subplot(2,8,11);xline(log(par_fast(3)));xlabel('parameter value');title('m');hold on
subplot(2,8,12);xline(log(par_fast(4)));xlabel('parameter value');title('x0');hold on
subplot(2,8,13);xline(log(par_fast(5)));xlabel('parameter value');title('y0');hold on
subplot(2,8,14);xline(log(par_fast(6)));xlabel('parameter value');title('c');hold on
subplot(2,8,15);xlabel('parameter value');ylabel('negative log likelihood');title('sd x');hold on
subplot(2,8,16);xlabel('parameter value');title('sd y');hold on

for l = 8:runs
    
    subplot(2,8,1)
    plot(parval_l_slow{l,posmediansl(l)}(1,:),(parres_l_slow{l,posmediansl(l)}(1,:)));hold on
    subplot(2,8,2)
    plot(parval_l_slow{l,posmediansl(l)}(2,:),(parres_l_slow{l,posmediansl(l)}(2,:)));hold on
    subplot(2,8,3)
    plot(parval_l_slow{l,posmediansl(l)}(3,:),(parres_l_slow{l,posmediansl(l)}(3,:)));hold on
    subplot(2,8,4)
    plot(parval_l_slow{1,posmediansl(l)}(4,:),(parres_l_slow{l,posmediansl(l)}(4,:)));hold on
    subplot(2,8,5)
    plot(parval_l_slow{l,posmediansl(l)}(5,:),(parres_l_slow{l,posmediansl(l)}(5,:)));hold on
    subplot(2,8,6)
    plot(parval_l_slow{l,posmediansl(l)}(6,:),(parres_l_slow{l,posmediansl(l)}(6,:)));hold on
    subplot(2,8,7)
    plot(parval_l_slow{l,posmediansl(l)}(7,:),(parres_l_slow{l,posmediansl(l)}(7,:)));hold on
    subplot(2,8,8)
    plot(parval_l_slow{l,posmediansl(l)}(8,:),(parres_l_slow{l,posmediansl(l)}(8,:)));hold on
    subplot(2,8,9)
    plot(parval_l_fast{l,posmediansl(l)}(1,:),(parres_l_fast{l,posmediansl(l)}(1,:)));hold on
    subplot(2,8,10)
    plot(parval_l_fast{l,posmediansl(l)}(2,:),(parres_l_fast{l,posmediansl(l)}(2,:)));hold on
    subplot(2,8,11)
    plot(parval_l_fast{l,posmediansl(l)}(3,:),(parres_l_fast{l,posmediansl(l)}(3,:)));hold on
    subplot(2,8,12)
    plot(parval_l_fast{1,posmediansl(l)}(4,:),(parres_l_fast{l,posmediansl(l)}(4,:)));hold on
    subplot(2,8,13)
    plot(parval_l_fast{l,posmediansl(l)}(5,:),(parres_l_fast{l,posmediansl(l)}(5,:)));hold on
    subplot(2,8,14)
    plot(parval_l_fast{l,posmediansl(l)}(6,:),(parres_l_fast{l,posmediansl(l)}(6,:)));hold on
    subplot(2,8,15)
    plot(parval_l_fast{l,posmediansl(l)}(7,:),(parres_l_fast{l,posmediansl(l)}(7,:)));hold on
    subplot(2,8,16)
    plot(parval_l_fast{l,posmediansl(l)}(8,:),(parres_l_fast{l,posmediansl(l)}(8,:)),'DisplayName',['y = '+string(samplesY(l))]);hold on
    
end


%% Functions

function [par_opt, error_res, exitflag,output,solutions, opt_fun,lb,ub] = llhFits(par_true, lb, ub, tspanX, tspanY, n_starts, par_init)
% Input
%     par_true    True parameters for simulation
%     lb, ub      bounds for parameter estimation
%     days        total days of simulation
%     samples     amount of samples in respectively the X and Y compartment
%     spacing     sample spacing  1 = linear; 2 = log10; 3 = log
%     tail        fixed end of the spacing timepoints
%     n_starts    amount of optimization starts

% Output
%     par_opt     the optimized parameters
%     resnorm     value of the residuals
%     exitflag,output,solutions   multistart outputs
%     opt_fun     function used for optimization, including data


%simulated
[x_true, y_true] = intLea(par_true, tspanX, tspanY);

xnoise = normrnd(0, 0.1, 1, length(x_true));
ynoise = normrnd(0, 0.1, 1, length(y_true));

xdata = log(x_true) + xnoise;
ydata = log(y_true) + ynoise;

%use data information in boundaries for x0 and y0
lb(4) = lb(4) * exp(xdata(1)) ; ub(4) = ub(4) * exp(xdata(1)); par_init(4) = exp(xdata(1));
lb(5) = lb(5) * exp(ydata(1)) ; ub(5) = ub(5) * exp(ydata(1)); par_init(5) = exp(ydata(1));
ub(6) = ub(6) * exp(ydata(1));

% Get latin hypercube sampled starts

[par_starts,~] = lhsdesign_modified(n_starts,log(lb),log(ub));
startpoints = CustomStartPointSet(par_starts);

% optimization for integral
%set inputs for fmincon
opt_fun = @(par)negLLH(par,xdata,ydata,tspanX, tspanY);

problem = createOptimProblem('fmincon','x0',log(par_init),'objective',opt_fun,...
    'lb',log(lb),'ub',log(ub),'xdata',xdata,'ydata',ydata);

ms_llh = MultiStart('UseParallel',true,'Display','iter','FunctionTolerance',1e-8);%,'PlotFcns'@gsplotbestf,'FunctionTolerance',1e-8);

[par_opt,error_res,exitflag,output,solutions] = ...
    run(ms_llh,problem,startpoints)

end


function [spac_days] = setSpacing(days, n_samples, tail, distribution)
%Input:
%   days        amount of days in total
%   n_samples   amount of samples
%   tail        fixed points at end of the spacing
%   distribution distritbution of the days, linear, log or log10

%Output:
%   spac_days   spacing of days

    if distribution == 1 %'linear'
        spacing = [linspace(0, days, n_samples),tail];
    elseif distribution == 2 % 'log10'
        spacing = [0, logspace(0, log10(days), n_samples-1),tail];
    elseif distribution == 3 %'log'
        spacing = [exp(linspace(0,log(days),n_samples)),tail];
        
    end
    
    spac_days = round(spacing);
    
    for i = 1:(length(spac_days)-1)
        if spac_days(i) == spac_days(i+1)
            spac_days(i+1) = spac_days(i+1) + 1;
        end
    end
      
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


function [negllh] = negLLH(par, xdata, ydata, tspanX, tspanY)
%calculates the negative log likelihood for minimization
%outputs:                          
%   sensitivity = output feature, squared difference of both outputs

%inputs:
%   tspan       timespan for integral
%   par =       adapted input parameters, including standard deviation
%   data =      data vectors for x and y respectively

% determine output for adapted parameters
[x_out, y_out] = intLea(exp(par(1:6)), tspanX, tspanY);

sigmX = par(7); 
sigmY = par(8);

%determine value of output of interest
featX = zeros(length(tspanX),1);
featY = zeros(length(tspanY),1);  

featX = (xdata - log(x_out)); %calculate differences
featY = (ydata - log(y_out)); %calculate differences

llhX = - 0.5 * log( 2 * pi * sigmX^2) - 0.5 * (featX.^2/sigmX^2);
llhY = - 0.5 * log( 2 * pi * sigmY^2) - 0.5 * (featY.^2/sigmY^2);
%http://jrmeyer.github.io/machinelearning/2017/08/18/mle.html

negllh = -1*(sum(llhX) + sum(llhY)) ; % possibility to apply weight factors

end 


function [parval, parres] = Profiles_llh(par_true, par_opt, timeframe, num, lb, ub,llhfun)
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
    originalpar = par_opt; %save original parameters

    parnum = length(par_opt);
    
    parres = zeros(parnum,num) ;
    parval = zeros(parnum,num) ;

    for j = 1:parnum
        
        parval(j,:) = linspace(lb(j), ub(j),num);
        
        for i = 1:num
            par_opt(j) = parval(j,i);
            
            if j <= 6
                parres(j,i) = sumofsquares(par_true,par_opt,timeframe);
            else 
                parres(j,i) = llhfun(par_opt); %%! TO DO: MEASURE AANPASSEN
            end
            
        end
        
        par_opt = originalpar; %reset to original input parameters
        
    end 
    
end


function [SSE] = sumofsquares(par_true,par_opt,timeframe)
%Calculates the differences between the data and simulation output
%can be used as input for lsqnonlin/least squares optimizer
%outputs:                          
%   sensitivity = output feature, differences of both outputs

%inputs:
%   tspan       timespan for integral
%   par =       adapted input parameters
%   data =      data vectors for x and y respectively

% determine output for adapted parameters
[x_out, y_out] = intLea(exp(par_opt), timeframe, timeframe);
[x_true, y_true] = intLea(exp(par_true), timeframe, timeframe);

maxX = max(x_true);
maxY = max(y_true);

normx_true = x_true/maxX; normx_out = x_out/maxX;
normy_true = y_true/maxY; normy_out = y_out/maxY;

%determine value of output of interest
featX = zeros(length(timeframe),1);
featY = zeros(length(timeframe),1);  

featX = sum((normx_out - normx_true).^2); %calculate differences
featY = sum((normy_out - normy_true).^2); %calculate differences

SSE = featX + featY;

end