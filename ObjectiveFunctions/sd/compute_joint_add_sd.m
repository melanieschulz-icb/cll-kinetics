function [sd_x, sd_y, sd_z, sd_hb] = compute_joint_add_sd(par, xdata, ydata, zdata, tspanX, tspanY,tspanZ,model,refill,hbdata, tspanYPlats,hbmodel,model_pb_pre)

    [featX, featY, featZ,featHB] = normedDifferencesAdditiveNoise(par,xdata,ydata,zdata,tspanX, tspanY,tspanZ, model,refill,hbdata, tspanYPlats,hbmodel,model_pb_pre);
    [sd_x, sd_y, sd_z, sd_hb] = compute_joint_sd_fromFeatures(featX, featY, featZ, featHB);
end