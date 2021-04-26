function bm_total = getBMTotalCount(zdata, z0, refill)
if refill
    bm_total=(z0/zdata(1)).*ones(1,length(zdata));
else
    bm_total=((1-zdata(1))/zdata(1))*z0./(1-zdata);
end

end