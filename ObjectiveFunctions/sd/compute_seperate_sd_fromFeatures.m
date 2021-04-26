function [sd_x, sd_y, sd_z, sd_hb] = compute_seperate_sd_fromFeatures(featX, featY, featZ, featHB)
    

    sd_x=sqrt((1/length(featX))*(norm(featX)^2));
    sd_y=sqrt((1/length(featY))*(norm(featY)^2));
    sd_z=sqrt((1/length(featZ))*(norm(featZ)^2));
    sd_hb=sqrt((1/length(featZ))*(norm(featHB)^2));
    
    sd_x=max(isnan(sd_x),sd_x);
    sd_y=max(isnan(sd_y),sd_y);
    sd_z=max(isnan(sd_z),sd_z);
    sd_hb=max(isnan(sd_hb),sd_hb);
    
end
