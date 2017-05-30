% European digital/vanilla call option with asset price following 1-dimensional
% geometric Brownian motion (Euler-Maruyama approximation). Calculation of
% gamma,Vega,vanna,delta and their variances using Vibrato Monte Carlo.
%
% This method is described in https://arxiv.org/pdf/1606.06143.pdf. We use
% a direct implementation here as opposed to using AD. Equations referenced
% are from this paper.
%
% Values in brackets were test values I used.
% mu - drift (0.05)
% sig - volatility (0.1)
% T     - time interval (1)
% S0    - initial asset value (50)
% K     - strike price (55)
% M     - number of Monte Carlo paths (10000)
% d   - number of final Z variables (10)
% N     - number of timesteps (100)
% flavor - type of option payoff (currently valid is digital or vanilla)
% greek - greek desired (currently valid is delta, gamma, vega, vanna)
function [value, variance] =... 
    Vibrato2ndOrder(mu, sig, T, S0, K, M, d, N, flavor, greek)

h = T/N;

S = S0*ones(1,M);
tanS0 = ones(1,M);
tanSig = zeros(1,M);
tanS0Sig = zeros(1,M);

%First lets Euler-Maruyama the process and associated tangent processes.
%Note the second order tangent process with respect to S0 will be
%identically zero, as the intial condition will be zero, and it takes the
%same form as tanS0
for n=1:(N-1)
   dW = sqrt(h)*randn(1,M);
   tanSig = tanSig.*(1+mu*h+sig*dW) + S.*dW;
   tanS0Sig = tanS0Sig.*(1+mu*h+sig*dW) + tanS0.*dW;
   tanS0 = tanS0.*(1+mu*h+sig*dW);
   S = S.*(1+mu*h+sig*dW);
end

Z = randn(d,M);

muw = S*(1+mu*h); %array of drift of final timesteps for each path
sigmaw = S*(sig*sqrt(h)); %array of vols of final timestep

dmuw = tanS0*(1+mu*h); % S0 tangent drifts
dsigmaw = tanS0*(sig*sqrt(h)); % vols

vmuw = tanSig*(1+mu*h); %Same for the sigma tangent
vsigmaw = tanSig*(sig*sqrt(h)) + S*sqrt(h);

vamuw = tanS0Sig*(1+mu*h); %Same for the S0-sigma tangent
vasigmaw = tanS0Sig*(sig*sqrt(h)) + tanS0*sqrt(h);

Sp=ones(d,1)*muw+(ones(d,1)*sigmaw).*Z;
Sm=ones(d,1)*muw-(ones(d,1)*sigmaw).*Z;
Sneutral=ones(d,1)*muw;
%The above quantities are used for the calculation of antithetic variables
%later, to speed up convergence

switch flavor
    case 'vanilla'
        fSp = exp(-mu*T)*(0.5)*(1+ sign(Sp-K)).*(Sp-K);
        fSm = exp(-mu*T)*(0.5)*(1 + sign(Sm - K)).*(Sm-K); %These are payoff values
        fmuw = exp(-mu*T)*(0.5)*(1 + sign(Sneutral - K)).*(Sneutral-K);
    case 'digital'
        fSp = exp(-mu*T)*(0.5)*(1 + sign(Sp - K));
        fSm = exp(-mu*T)*(0.5)*(1 + sign(Sm - K)); 
        fmuw = exp(-mu*T)*(0.5)*(1 + sign(Sneutral - K));
    otherwise
        error('Invalid payoff time. Currently "digital" or "vanilla" is supported')
end
    
oddAntithetic = fSp - fSm; % For terms where the polynomial in Z is odd
evenAntithetic = fSp - 2*fmuw + fSm; % When its even, subtract the 2 to
%stay O(1)

switch greek %Section 2.2-2.3 for explanation of 1st order, 2.4 for 2nd
    case 'delta'
        Ymu = Z.*(ones(d,1)*(1./sigmaw))*(1/2).*(oddAntithetic);
        Ysigma = (Z.^2-1).*(ones(d,1)*(1./sigmaw))*(1/2).*(evenAntithetic);
        deltMat = (ones(d,1)*dmuw).*Ymu + (ones(d,1)*dsigmaw).*Ysigma;
        valsDelta = 1/d*sum(deltMat);
        value = sum(valsDelta)/M;
        variance=(1/M)*(1/M*sum(valsDelta.^2)-(value)^2)+1/(M*d)...
        *1/M*sum(sum(deltMat.^2)/d-valsDelta.^2);
    case 'vega'
        Ymu = Z.*(ones(d,1)*(1./sigmaw))*(1/2).*(oddAntithetic);
        Ysigma = (Z.^2-1).*(ones(d,1)*(1./sigmaw))*(1/2).*(evenAntithetic);
        vegaMat = (ones(d,1)*vmuw).*Ymu + (ones(d,1)*vsigmaw).*Ysigma;
        valsVega = 1/d*sum(vegaMat);
        value = sum(valsVega)/M;
        variance=(1/M)*(1/M*sum(valsVega.^2)-(value)^2)+1/(M*d)...
        *1/M*sum(sum(vegaMat.^2)/d-valsVega.^2);
    case 'gamma'
        Ymusquare = (Z.^2-1).*(ones(d,1)*(1./(sigmaw.^2)))*(1/2).*(evenAntithetic);
        Ysigmasquare = (Z.^4-5.*Z.^2 - 2).*(ones(d,1)*(1./(sigmaw.^2)))*(1/2).*(evenAntithetic);
        Ymixed = (Z.^3-3.*Z).*(ones(d,1)*(1./(sigmaw.^2)))*(1/2).*(oddAntithetic);
        gammaMat = (ones(d,1)*(dmuw.^2)).*Ymusquare + ...
        (ones(d,1)*(dsigmaw.^2)).*Ysigmasquare + ...
        2.*(ones(d,1)*(dmuw.*dsigmaw)).*Ymixed;
        valsGamma = 1/d*sum(gammaMat);
        value = sum(valsGamma)/M;
        variance=(1/M)*(1/M*sum(valsGamma.^2)-(value)^2)+1/(M*d)...
        *1/M*sum(sum(gammaMat.^2)/d-valsGamma.^2);
    case 'vanna'
        Ydoublemu = (Z).*(ones(d,1)*(1./(sigmaw)))*(1/2).*(oddAntithetic);
        Ydoublesigma = (Z.^2-1).*(ones(d,1)*(1./(sigmaw)))*(1/2).*(evenAntithetic);
        Ymusquare = (Z.^2-1).*(ones(d,1)*(1./(sigmaw.^2)))*(1/2).*(evenAntithetic);
        Ysigmasquare = (Z.^4-5.*Z.^2 - 2).*(ones(d,1)*(1./(sigmaw.^2)))*(1/2).*(evenAntithetic);
        Ymixed = (Z.^3-3.*Z).*(ones(d,1)*(1./(sigmaw.^2)))*(1/2).*(oddAntithetic);
        vannaMat = (ones(d,1)*(vamuw)).*Ydoublemu + ...
        (ones(d,1)*(dmuw.*vmuw)).*Ymusquare + ...
        (ones(d,1)*(vasigmaw)).*Ydoublesigma + ...
        (ones(d,1)*(dsigmaw.*vsigmaw)).*Ysigmasquare + ...
        (ones(d,1)*(dmuw.*vsigmaw + vmuw.*dsigmaw)).*Ymixed;
        valsVanna = 1/d*sum(vannaMat);
        value = sum(valsVanna)/M;
        variance=(1/M)*(1/M*sum(valsVanna.^2)-(value)^2)+1/(M*d)...
        *1/M*sum(sum(vannaMat.^2)/d-valsVanna.^2);
    otherwise
        error('Invalid Greek. Delta, gamma, vega and vanna are currently supported');
end
