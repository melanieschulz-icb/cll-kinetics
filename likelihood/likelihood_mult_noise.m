function [likelihood]=likelihood_mult_noise(par,xdata,ydata,zdatanorm,tspanX, tspanY,tspanZ, model, refill,seperate_sd)
    likelihood = exp(-0.5*diff_log_summed_sd(par,xdata,ydata,zdatanorm,tspanX, tspanY,tspanZ, model, refill,seperate_sd));
end