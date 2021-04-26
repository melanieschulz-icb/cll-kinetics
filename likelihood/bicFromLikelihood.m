function [bic]=bicFromLikelihood(parNumb, dataNumb, likelihood)
    bic=parNumb*log(dataNumb)-2*log(likelihood);
end