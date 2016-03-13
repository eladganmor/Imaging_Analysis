%generate a random calcium trace
function [caTrace, r, rate, stim] = genCaTrace(stimMat,tuning,bias)
    global presets
    
    fStim = tuning*stimMat;
    fStim = [zeros(size(fStim,1),1), fStim(:,1:end-1)];%stimulus at time t determines rate at time t+dt
    rate = exp(fStim + bias);
    r = poissrnd(rate);%spikes

    decayFun = presets.a*exp(-(0:100)*presets.dt/presets.tau);
    
    %%convert spikes to calcium
    caTrace = conv(r,decayFun,'valid');
    caTrace = caTrace + normrnd(zeros(size(caTrace)),presets.sig);

    %make sure r and caTrace correspond in indices
    r = r(:,length(decayFun):end);
    stim = stimMat(:,length(decayFun):end);
    rate = rate(:,length(decayFun):end);
end