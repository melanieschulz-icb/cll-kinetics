function [sd_x, sd_y, sd_z] = compute_mult_sd(par, xdata, ydata, zdata, tspanX, tspanY,tspanZ,model, refill,seperate_sd)
    if seperate_sd
        [sd_x, sd_y, sd_z]=compute_seperate_mult_sd(par, xdata, ydata, zdata, tspanX, tspanY,tspanZ,model,refill);
    else
        [sd_x, sd_y, sd_z]=compute_joint_mult_sd(par, xdata, ydata, zdata, tspanX, tspanY,tspanZ,model,refill);
    end
end