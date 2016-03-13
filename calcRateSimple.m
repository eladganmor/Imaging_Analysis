%calculate rate given parameters and features
function rate = calcRateSimple(params,fVals)
    rate = exp(params*fVals);
end