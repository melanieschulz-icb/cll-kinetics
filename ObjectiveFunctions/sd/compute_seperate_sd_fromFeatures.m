function [sd_x, sd_y, sd_z] = compute_seperate_sd_fromFeatures(featX, featY, featZ)
    
    sd_x=sqrt((1/length(featX))*(norm(featX)^2));
    sd_y=sqrt((1/length(featY))*(norm(featY)^2));
    sd_z=sqrt((1/length(featZ))*(norm(featZ)^2));
end
