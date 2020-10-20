function [likelihood]=likelihood_add_noise(par,xdata,ydata,zdatanorm,tspanX, tspanY,tspanZ, model, refill,seperate_sd)
    likelihood = exp(-0.5*diff_integral_summed_sd(par,xdata,ydata,zdatanorm,tspanX, tspanY,tspanZ, model,refill, seperate_sd));
end