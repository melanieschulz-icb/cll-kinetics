function [likelihood]=likelihood_add_noise(par,xdata,ydata,zdatanorm,tspanX, tspanY,tspanZ, model, refill,seperate_sd,hbdata,tYspanPlats,hbmodel,model_pb_pre)
    likelihood = exp(-0.5*diff_integral_summed_sd(par,xdata,ydata,zdatanorm,tspanX, tspanY,tspanZ, model,refill, seperate_sd,hbdata,tYspanPlats,hbmodel,model_pb_pre));
end