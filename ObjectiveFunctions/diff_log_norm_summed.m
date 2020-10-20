function [sensitivity] = diff_log_norm_summed(par,xdata,ydata,zdata,tspanX, tspanY,tspanZ, model, refill)
%Calculates the differences between the data and simulation output
%can be used as input for lsqnonlin/least squares optimizer
%uses multiplicative noise (logs the data)
%normalizes the logged data
%outputs:                          
%   sensitivity = output feature, differences of both outputs

%inputs:
%   tspan       timespan for integral
%   par =       adapted input parameters
%   data =      data vectors for x and y respectively

%calculate normed differences
[featX, featY, featZ] = normedDifferencesMultiplicativeNoise(par,xdata,ydata,zdata,tspanX, tspanY,tspanZ, model, refill);

%sum squared differences
sensitivity = norm([featX, featY, featZ])^2;

end