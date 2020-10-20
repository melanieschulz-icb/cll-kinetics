
%patients you want to fit; 12 will be excluded as ibrutinib dose was paused
%(available patients: 1-30)
patient=[1];
runn=2;
%params
%xmulti=[0.060520314	0.02236084	0.003954597	3.29507E+12	6.19671E+11	2.69818E+11]
xmulti=[0.03670422	0.000201142	0.002571473	3.35665E+11	6.14665E+11	3.15997E+12	1.04563E+12];
%xmulti=[0.0185	0.0033	0.0174	1.806e+12	5.794e+11	1.805e+11];

%assumed volume of a CLL cell (in fl)
cll_volume=166;

%tissue measurements to be taken into account (max 3 available, Burger
%takes only one)
tss_measurements=2;

%if one should use the complex computation of the tissue cll burden
use_complex_tss=1;

%if one should normaize when using the least square optimizer
normalize_lls=1;

%weight for the second tissue measurement (if second measurement is taken
%into account)
weight=1;

folder='/Figures/08_22/imp/'

[num, txt, ~] = xlsread('Heavywater-all counts-decoded.xlsx');
[chars, charDesc,~]=xlsread('Heavy Water patients characterists.xlsx');

lb = [1e-3 1e-5 1e-6 0.7 0.7 1e9]; 
ub = [1 1 0.5 1.3 1.3 1e12]; 

n_starts = 2000;

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

lsqfun = @(par)diff_integral_norm_summed(par,xdata,ydata,tspanX, tspanY, weight);

[parval, parres] = Profiles_lsq(lsqfun, xmulti, 10*n_starts, lb, ub);

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

f=gcf;
exportgraphics(f,[pwd,strcat(folder,plotTitle.String,'.png')],'Resolution',300);
savefig(f,[pwd,strcat(folder,plotTitle.String)],'compact')

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