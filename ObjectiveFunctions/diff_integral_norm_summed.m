function [sensitivity] = diff_integral_norm_summed(par,xdata,ydata,zdata,tspanX, tspanY,tspanZ, model, refill)
%Calculates the differences between the data and simulation output
%can be used as input for lsqnonlin/least squares optimizer
%normalizes different datasources with l2 norm
%outputs:                          
%   sensitivity = output feature, differences of both outputs

%inputs:
%   tspanX,Y,Z       timespan for measurements
%   par =       adapted input parameters
%   dataX,Y,Z =      data vectors for x,y and z respectively


%calculate normed differences
[featX, featY, featZ] = normedDifferencesAdditiveNoise(par,xdata,ydata,zdata,tspanX, tspanY,tspanZ, model, refill);


%sum the squared differences
sensitivity = norm([featX, featY, featZ])^2;

end