function V = AnalyticVanillaCall(T, S, sigma, r, K, ans )

%The Greeks for a vanilla European call as per B-S
%T - time to expiry
%S - Asset Price
%sigma - vol
%r - risk free rate
%K - Strike
%ans - What we want, options are 'value', 'delta'
%'gamma', 'vega', 'vanna'

S  = max(1e-16,S);     % no division by 0
K  = max(1e-16,K);     
d1 = (log(S/K)+(r+(sigma^2)/2)*T)/(sigma * sqrt(T));
d2 = d1 - sigma*(sqrt(T));
%d2 = ( log(S) - log(K) + (r-0.5*sigma^2)*T ) / (sigma*sqrt(T));
switch ans
case 'value'
V = (S *N(d1)) - (K*exp(-r*T)*N(d2));
case 'delta'
V = N(d1);
case 'gamma'
V =  (1./sqrt(2*pi))* (exp(-0.5*d1.^2)/( S*sigma*sqrt(T)) );
%V = exp(-r*T)*(exp(-0.5*d2.^2)./sqrt(2*pi))./(sigma*sqrt(T)*S);
%V = -exp(-r*T)
case 'vega'
V = S*exp(-0.5*d1.^2)*sqrt(T)*(1/sqrt(2*pi));
case 'vanna'
%V = -exp(-r*T)*(exp(-0.5*d2.^2).*(-d2).*(1./(S*sigma*sqrt(T)))) ...
%    + exp(-r*T)*(exp(-0.5*d2.^2)./(sigma*sqrt(T)*S*(sigma+sqrt(T))));
V = -exp(-r*T)*exp(-0.5*d1.^2)*(1/sqrt(2*pi))*(d2./sigma);
otherwise
error('Wrong value field')
end
%
% Normal cumulative distribution function
%
function ncf = N(x)
%ncf = 0.5*(1+erf(x/sqrt(2)));
xr = real(x);
xi = imag(x);
if abs(xi)>1e-10
error 'imag(x) too large in N(x)'
end
ncf = 0.5*(1+erf(xr/sqrt(2))) ...
+ 1i*xi.*exp(-0.5*xr.^2)/sqrt(2*pi);