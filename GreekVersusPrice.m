function GreekVersusPrice(flavor, spots, greek, begin, upper, mu, sig, T, K, M, d, N)

%Plots the evolution of a certain Greek with price. Choices are delta,
%gamma,vega, vanna. Flavor of option is also supported (digital or vanilla)
%begin/upper - interval of spot over which to plot greek
%spots - desired number of spots in interval
%mu - risk-free rate
%sig - vol
%T - time to expiry
%K - strike
%M - number of Monte Carlo paths used in each simulation
%d - number of random variables for each final timestep
%N - number of timesteps in discretisation of path

InitialSpots = linspace(begin,upper,spots);
Analytics = zeros(1, spots);
MonteCarlos = zeros(1, spots);

for i = 1:spots
    Analytics(i) = AnalyticVanillaCall(T,InitialSpots(i), vol,mu,K,greek);
    vibArray = Vibrato2ndOrder(mu, sig, T, InitialSpots(i), K, M, d, N, flavor, greek);
    MonteCarlos(i) = vibArray(1);
end

plot(InitialSpots, Analytics);
hold on;
plot(InitialSpots, MonteCarlos);
legend('Analytic value', 'VMC value');
