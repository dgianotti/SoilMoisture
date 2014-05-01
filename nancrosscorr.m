function [XC, lags] = nancrosscorr(X,Y)

% Make sure X and Y are the same length:
if ~isvector(X) || ~isvector(Y) || (length(X) ~= length(Y))
    error('The inputs to the function nancrosscorr(X,Y) should be two vectors of equal length. Aborting!');
end

N = length(X);

max_lags = 96;

if N < max_lags
    lags = (-N+3):(N-3);
else
    lags = -max_lags:max_lags;
end

XC = zeros(size(lags));

for i = 1:length(lags)
    lag = lags(i);
    
    if lag < 0
        x = X(1:(N+lag));
        y = Y( (1-lag):end );
    elseif lag == 0
        x = X;
        y = Y;
    else % lag > 0
        x = X( (1+lag):end );
        y = Y(1:(N-lag));
    end
    
    
    r = corrcoef( x( ~isnan(x) & ~isnan(y) ), y( ~isnan(x) & ~isnan(y) ) );
    XC(i) = r(2,1);
end



end