function [sensitivity] = diff_integral_summed_sd(par,xdata,ydata,zdata,tspanX, tspanY,tspanZ, model,refill,seperate_sd)
%Calculates the differences between the data and simulation output
%can be used as input for lsqnonlin/least squares optimizer
%Optimizes with additive noise, uses hierarchical optimization
%uses l2norm for normalization
%outputs:                          
%   sensitivity = output feature, differences of both outputs

%inputs:
%   tspanX,Y,Z       timespan for integral
%   par =       adapted input parameters
%   dataX,Y,Z =      data vectors for x,y and z respectively

%compute normed Differences
[featX, featY, featZ] = normedDifferencesAdditiveNoise(par,xdata,ydata,zdata,tspanX, tspanY,tspanZ, model, refill);

%compute sd_s
[sd_x, sd_y, sd_z] = compute_sd_fromFeatures(featX, featY, featZ, seperate_sd);

%compute logged likelihood
sensitivity = norm([featX./sd_x, featY./sd_y, featZ./sd_z])^2+ ...
    +length(xdata)*log(2*pi*(sd_x^2))+...
    +length(ydata)*log(2*pi*(sd_y^2))+...
    +length(zdata)*log(2*pi*(sd_z^2));


end