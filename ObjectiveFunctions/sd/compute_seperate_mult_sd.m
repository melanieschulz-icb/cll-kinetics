function [sd_x, sd_y, sd_z] = compute_seperate_mult_sd(par,xdata,ydata,zdata,tspanX, tspanY, tspanZ, model,refill)
%Calculates the sd_s between the data and simulation output
%uses multiplicative noise (logs data)
%outputs:                          
%   sd_x,sd_y, sd_z

%inputs:
%   tspanX,Y,Z =      timespan for integral
%   par =       adapted input parameters
%   dataX,Y,Z =      data vectors for x and y respectively
%   model = model to be used

% determine output for adapted parameters
[featX, featY, featZ] = normedDifferencesMultiplicativeNoise(par,xdata,ydata,zdata,tspanX, tspanY,tspanZ, model,refill);

%compute sd_s
[sd_x, sd_y, sd_z]=compute_seperate_sd_fromFeatures(featX, featY, featZ);

end