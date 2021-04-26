function [featX, featY, featZ, featHB] = normedDifferencesAdditiveNoise(par,xdata,ydata,zdata,tspanX, tspanY,tspanZ, model, refill,hbdata, tspanYPlats, hb_model,model_pb_pre)
%Calculates the normed differences between the data and simulation output
%normalizes different datasources with l2 norm
%outputs:                          
%   [featX, featY, featZ] = output feature, differences of both outputs

%inputs:
%   tspanX,Y,Z       timespan for measurements
%   par =       adapted input parameters
%   dataX,Y,Z =      data vectors for x,y and z respectively

% determine differences
[featX, featY, featZ, featHB] = differencesAdditiveNoise(par,xdata,ydata,zdata,tspanX, tspanY,tspanZ, model, refill,hbdata, tspanYPlats,hb_model,model_pb_pre);


%normalize
% featX =featX./max(xdata);
% featY =featY./max(ydata); 
% featZ =featZ./max(zdata);
% featHB=featHB./max(hbdata);

featX =featX./norm(xdata);
featY =featY./norm(ydata); 
featZ =featZ./norm(zdata);
featHB=featHB./norm(hbdata);
end