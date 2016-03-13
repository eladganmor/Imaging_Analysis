%%In this version no separable filters, gains or anything like that
function [optParams, paramsErr, lambda] = searchBestParamsSimple(data,varargin)
    
[params0, paramsErr, tuningFun] = initParams(varargin);%initialize parameter values
    minMe();
    
    function ll = minMe()
        options = optimset('GradObj','on','Hessian','on','display','off','UseParallel','always','MaxIter',1e3);
        clear logLikelihoodSumGivenCaSimple logLikelihoodSumGivenCaParametric %clear persistent variables
        
        [optParams,ll,~,~,~,hessian] = fminunc(@innerMin,params0,options);
        paramsErr = sqrt(1./diag(hessian));
        
        %%%%%sub-routines
        function [llInner, deriv, H] = innerMin(params)
            if isempty(tuningFun)%non-parametric
                [llInner, deriv, H, lambda] = logLikelihoodSumGivenCaSimple(params, data);
            else%parametric
                [llInner, deriv, H, lambda] = logLikelihoodSumGivenCaParametric(params, data, tuningFun);
            end
            llInner = -llInner;
            deriv = -deriv;
            H = -H;
            
            %fminunc doesn't like single
            llInner = double(llInner);
            deriv = double(deriv);
            H = double(H);
        end
    end

    function [params, initParamsErr, tunFun] = initParams(inputs)
        params = []; tunFun = [];
        if length(inputs)>=1
            params = inputs{1};%%May be empty
            if length(inputs)==2
                tunFun = inputs{2};
            else
                tunFun = [];
            end
        end
        if isempty(params)
            numOfParams = size(data.fVals,1);%1 + presets.numOfFuncsSelf + (data(1).nCells-1)*presets.numOfFuncsInter;
            params = zeros(1,numOfParams);
        end
        initParamsErr = zeros(size(params));
    end
end