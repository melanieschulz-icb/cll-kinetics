function [sd_x, sd_y, sd_z] = compute_joint_sd_fromFeatures(featX, featY, featZ)
    
    sd_x=sqrt((1/length([featX, featY, featZ]))*(norm([featX, featY, featZ])^2));
    sd_y=sd_x;
    sd_z=sd_x;
end