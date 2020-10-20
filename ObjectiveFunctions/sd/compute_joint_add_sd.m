function [sd_x, sd_y, sd_z] = compute_joint_add_sd(par, xdata, ydata, zdata, tspanX, tspanY,tspanZ,model,refill)

    [featX, featY, featZ] = normedDifferencesAdditiveNoise(par,xdata,ydata,zdata,tspanX, tspanY,tspanZ, model,refill);
    [sd_x, sd_y, sd_z] = compute_joint_sd_fromFeatures(featX, featY, featZ);
end