%calculate rate with parametric tuning function f
function [rate, dLogRate, d2LogRate] = calcRateParametric(xVal,params,f)
    [logRate, dLogRate, d2LogRate] = f(xVal,params);
    rate = exp(logRate);
end