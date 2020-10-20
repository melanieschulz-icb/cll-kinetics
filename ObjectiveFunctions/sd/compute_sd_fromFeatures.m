function [sd_x, sd_y, sd_z] = compute_sd_fromFeatures(featX, featY, featZ, seperate_sd)
    
    if seperate_sd
        [sd_x, sd_y, sd_z] = compute_seperate_sd_fromFeatures(featX, featY, featZ);
    else
        [sd_x, sd_y, sd_z] = compute_joint_sd_fromFeatures(featX, featY, featZ);
    end
end