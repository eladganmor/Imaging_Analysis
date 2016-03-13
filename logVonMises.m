%Von mises tuning function y = k*cos(2*(x - mu)) + b;
%input:
%1) x - orientation values
%2) parameter vector - [b, mu, k]

%output:
%1) y - tuning values
%2) d - first derivative
%3) d2 - second derivative
function [y, d, d2] = logVonMises(x,params)
    b = params(1); mu = params(2); k = exp(params(3));
    
    y = k*cos(2*(x - mu)) + b;%I multiply by 2 because all orientations lie in [0 pi] not [0 2*pi]
    
    d = [ones(1,length(x));
        2*k*sin(2*(x - mu));
        k*cos(2*(x - mu))];
    
    d2 = zeros(3,3,length(x));
    d2(2,2,:) = -2*k*cos(2*(x - mu));%derive twice wrt mu
    d2(2,3,:) = 2*k*sin(2*(x - mu));%dmu dk
    d2(3,2,:) = d2(2,3,:);
    d2(3,3,:) = k*cos(2*(x - mu));
    
    %%Handle blank (no stimulus)
    blankInds = isnan(x);
    y(blankInds) = b;
    d(:,blankInds) = 0;
    d2(:,:,blankInds) = 0;
end