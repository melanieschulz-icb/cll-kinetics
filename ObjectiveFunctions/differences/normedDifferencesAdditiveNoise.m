function [featX, featY, featZ] = normedDifferencesAdditiveNoise(par,xdata,ydata,zdata,tspanX, tspanY,tspanZ, model, refill)
%Calculates the normed differences between the data and simulation output
%normalizes different datasources with l2 norm
%outputs:                          
%   [featX, featY, featZ] = output feature, differences of both outputs

%inputs:
%   tspanX,Y,Z       timespan for measurements
%   par =       adapted input parameters
%   dataX,Y,Z =      data vectors for x,y and z respectively

% determine differences
[featX, featY, featZ] = differencesAdditiveNoise(par,xdata,ydata,zdata,tspanX, tspanY,tspanZ, model, refill);


%normalize
featX = featX./norm(xdata);
featY = featY./norm(ydata); 
featZ = featZ./norm(zdata);

end