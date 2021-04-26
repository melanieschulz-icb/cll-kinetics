function [sd_x, sd_y, sd_z, sd_hb] = compute_sd_fromFeatures(featX, featY, featZ,featHB, seperate_sd)
    
    if seperate_sd
        [sd_x, sd_y, sd_z, sd_hb] = compute_seperate_sd_fromFeatures(featX, featY, featZ, featHB);
    else
        [sd_x, sd_y, sd_z, sd_hb] = compute_joint_sd_fromFeatures(featX, featY, featZ, featHB);
    end
end