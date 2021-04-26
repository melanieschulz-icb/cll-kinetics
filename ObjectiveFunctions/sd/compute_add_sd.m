function [sd_x, sd_y, sd_z, sd_hb] = compute_add_sd(par, xdata, ydata, zdata, tspanX, tspanY,tspanZ,model, refill,seperate_sd, hbdata, tspanYPlats,hbmodel,model_pb_pre)
    if seperate_sd
        [sd_x, sd_y, sd_z, sd_hb]=compute_seperate_add_sd(par, xdata, ydata, zdata, tspanX, tspanY,tspanZ,model,refill,hbdata, tspanYPlats,hbmodel,model_pb_pre);
    else
        [sd_x, sd_y, sd_z, sd_hb]=compute_joint_add_sd(par, xdata, ydata, zdata, tspanX, tspanY,tspanZ,model,refill,hbdata, tspanYPlats,hbmodel,model_pb_pre);
    end
end