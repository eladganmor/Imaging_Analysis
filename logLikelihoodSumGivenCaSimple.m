%%Calculates the log likelihood given the data and parameters
%%returns log likelihood, derivatives, Hessian and rate
function [ll, deriv, H, lam] = logLikelihoodSumGivenCaSimple(params,data)
    global presets
    persistent pYGivenX x xFact xRep%stuff that doesn't depend on parameters and will be used every time
    maxSpikes = presets.rMax;%maximal number of spike per time bin that we will consider
    
    f = data.fVals(:,2:end);%this is to align features (stimuli) with the response at the next time step
    
    lam = calcRateSimple(params,f);
    
    %%prepare some variables and match in size
    lamRep = repmat(lam,maxSpikes+1,1);
    
    calculatePersistent();%calculate presistent variables that don't depend on parameters
    
    pX = exp(-lamRep).*lamRep.^xRep./xFact;
    pXY = pX.*pYGivenX;
    sum_PXY = sum(pXY,1);
    ll = sum(log(sum_PXY),2);

    if nargout>1%%derivatives
        sumX_PXY = x'*pXY;
        sumX_PXY_over_sum_PXY = sumX_PXY./sum_PXY;
        gamma = sumX_PXY_over_sum_PXY - lam;
        
        dLogLam = calcdLogLam();
        
        deriv = dLogLam*gamma';
        %%Hessian
        if nargout>2%%hessian
            c = lam - x'.^2*pXY./sum_PXY + (sumX_PXY_over_sum_PXY).^2;
            H = -bsxfun(@times,c,dLogLam)*dLogLam';
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

    function d = calcdLogLam()
        d = f;%derivative of log rate relative to parameters as a function of time
    end
    
end