function [sensitivity] = diff_log_summed(par,xdata,ydata,zdata, tspanX, tspanY,tspanZ,model, refill)
%Calculates the differences between the data and simulation output
%can be used as input for lsqnonlin/least squares optimizer
%outputs:                          
%   sensitivity = output feature, differences of both outputs

%inputs:
%   tspanX,Y,Z       timespan for integral
%   par =       adapted input parameters
%   dataX,Y,Z =      data vectors for x,y,z respectively

%compute differences
[featX, featY, featZ]=differencesMultiplicativeNoise(par,xdata,ydata,zdata,tspanX, tspanY,tspanZ, model, refill);

%sum squared differences
sensitivity = norm([featX, featY, featZ])^2;

end