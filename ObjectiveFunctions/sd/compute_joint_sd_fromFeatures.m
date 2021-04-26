function [sd_x, sd_y, sd_z, sd_hb] = compute_joint_sd_fromFeatures(featX, featY, featZ, featHB)
    
    sd_x=sqrt((1/length([featX, featY, featZ, featHB]))*(norm([featX, featY, featZ, featHB])^2));
    sd_y=sd_x;
    sd_z=sd_x;
    sd_hb=sd_x;
end