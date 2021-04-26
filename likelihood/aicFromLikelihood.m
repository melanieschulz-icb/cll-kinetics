function [aic]=aicFromLikelihood(parNumb, likelihood)
    aic=2*parNumb-2*log(likelihood);
end