%%Calculates the log likelihood given the data, parameters, and a
%%parametric tuning function
%%returns log likelihood, derivatives, Hessian and rate
function [ll, deriv, H, lam] = logLikelihoodSumGivenCaParametric(params,data,tuningFun)
    global presets
    persistent pYGivenX x xFact xRep%stuff that doesn't depend on parameters and will be used every time
    maxSpikes = presets.rMax;%maximal number of spike per time bin that we will consider
    
    xVals = data.xVals(:,2:end);
    goodInds = ~isnan(xVals);%I use this to mark blanks because I don't want to try to fit the rate during blanks
    [lam, dLogLam, d2LogLam] = calcRateParametric(xVals,params,tuningFun);
    
    %%prepare some variables and match in size
    lamRep = repmat(lam,maxSpikes+1,1);
    
    calculatePersistent();
    
    pX = exp(-lamRep).*lamRep.^xRep./xFact;
    pXY = pX.*pYGivenX;
    sum_PXY = sum(pXY,1);
    ll = sum(log(sum_PXY(goodInds)),2);

    if nargout>1%%derivatives
        sumX_PXY = x'*pXY;
        sumX_PXY_over_sum_PXY = sumX_PXY./sum_PXY;
        gamma = sumX_PXY_over_sum_PXY - lam;
        
        dLogLam = dLogLam(:,goodInds);%I use this to mark blanks because I don't want to try to fit the rate during blanks
        gamma = gamma(goodInds);
        
        deriv = dLogLam*gamma';
        %%Hessian
        if nargout>2 && ~isempty(d2LogLam)%%hessian
            c = lam - x'.^2*pXY./sum_PXY + (sumX_PXY_over_sum_PXY).^2;
            c = c(goodInds);
            part1 = -bsxfun(@times,c,dLogLam)*dLogLam';%this part is free of second derivatives
            d2LogLam = d2LogLam(:,:,goodInds);
            part2 = squeeze(sum(bsxfun(@times,d2LogLam,reshape(gamma,1,1,[])),3));%this part has second derivatives
            H = part1 + part2;
        else
            H = [];
        end
    end
    
    
    %%%%%%%%%Sub routines
    function calculatePersistent()
        if isempty(pYGivenX)%%This is the first time, let's calculate a bunch of stuff that will be used over and over as is
            if isempty(x)
                x = single((0:maxSpikes)');%%number of spikes (unobserved)
            end
            xFact = factorial(x);%simply x!
            xRep = repmat(x,1,size(lamRep,2));%maxSpikes x T
            xFact = repmat(xFact,1,size(lamRep,2));
            %%Calculate the probability of seeing observed fluorescence given
            %%previously observed fluorescence and guessed number of spikes
            r = repmat(data.response,maxSpikes+1,1);
            yPred = r(:,1:end-1)*exp(-presets.dt/presets.tau) + presets.a*xRep;%the response we expect to see if there were x spikes in a bin
            %         r(:,1:end-1)*exp(-presets.dt/presets.tau) + presets.baseline*(1 - exp(-presets.dt/presets.tau)) + presets.a*xRep;%%in case of non-zero baseline
            if presets.sig==0
                pYGivenX = abs(yPred - r(:,2:end)) < 1e-3;%if there is no error (just for testing) than the observed calcium has to exactly match one (and only one) possibility of number of spikes
            else
                pYGivenX = exp(-0.5*((r(:,2:end) - yPred)/presets.sig).^2)/(sqrt(2*pi)*presets.sig);%Prob(Ca(t)|Ca(t-dt),x). The probability of seeing a certain calcium response at time t given the calcium level at time t-1 and that x spikes occured
            end
            %%Put stuff on the gpu if needed
            if isfield(presets,'useGPU') && presets.useGPU
                x = gpuArray(x);
                xRep = gpuArray(xRep);
                xFact = gpuArray(xFact);
                pYGivenX = gpuArray(pYGivenX);
            end
        end
    end

end