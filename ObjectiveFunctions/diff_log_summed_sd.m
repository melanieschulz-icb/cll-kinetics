function [sensitivity] = diff_log_summed_sd(par,xdata,ydata,zdata,tspanX, tspanY,tspanZ, model,refill ,seperate_sd)
%Calculates the differences between the data and simulation output
%can be used as input for lsqnonlin/least squares optimizer
%uses hierarchical optimization for optimizing sd
%uses multiplicative noise (logs data)
%outputs:                          
%   sensitivity = output feature, differences of both outputs

%inputs:
%   tspanX,Y,Z =      timespan for integral
%   par =       adapted input parameters
%   dataX,Y,Z =      data vectors for x and y respectively



%compute differences
[featX, featY, featZ] = normedDifferencesMultiplicativeNoise(par,xdata,ydata,zdata,tspanX, tspanY,tspanZ, model, refill);

%compute sd_s
[sd_x, sd_y, sd_z]=compute_sd_fromFeatures(featX, featY, featZ,[], seperate_sd);
sd_x=max(sd_x, 1e-20);

%compute sensitivity
sensitivity = norm([featX./sd_x, featY./sd_y, featZ./sd_z])^2+...
    +length(xdata)*log(2*pi*(sd_x^2))+...
    +length(ydata)*log(2*pi*(sd_y^2))+...
    +length(zdata)*log(2*pi*(sd_z^2));

end