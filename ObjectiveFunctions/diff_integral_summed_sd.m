function [sensitivity] = diff_integral_summed_sd(par,xdata,ydata,zdata,tspanX, tspanY,tspanZ, model,refill,seperate_sd,hbdata,tspanYPlats, hb_model,model_pb_pre)
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
[featX, featY, featZ, featHB] = normedDifferencesAdditiveNoise(par,xdata,ydata,zdata,tspanX, tspanY,tspanZ, model, refill, hbdata, tspanYPlats, hb_model,model_pb_pre);

%compute sd_s
% if ~isempty([tspanX, tspanY,tspanZ])
%     [sd_x, sd_y, sd_z] = compute_sd_fromFeatures(featX, featY, featZ, seperate_sd);
% else
%     sd_x=1;
%     sd_y=1;
%     sd_z=1;
% end
% if ~isempty(tspanYPlats)
%     sd_hb=compute_sd_fromFeatures(featHB,[],[], 0);
% else
%     sd_hb=1;
% end

[sd_x, sd_y, sd_z, sd_hb]=compute_sd_fromFeatures(featX,featY,featZ,featHB, seperate_sd);

%compute logged likelihood
sensitivity = norm([featX./sd_x, featY./sd_y, featZ./sd_z,featHB./sd_hb])^2+ ...
    +length(xdata)*log(2*pi*(sd_x^2))+...
    +length(ydata)*log(2*pi*(sd_y^2))+...
    +length(zdata)*log(2*pi*(sd_z^2))+...
    +length(hbdata)*log(2*pi*(sd_hb^2));
% sensitivity=norm([featX, featY, featZ, featHB])^2


end