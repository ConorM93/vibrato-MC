%Plots the evolution of a certain Greek with price. Choices are delta,
%gamma,vega, vanna


spots = 100;
InitialSpots = linspace(70,130,spots);
Analytics = zeros(1, spots);
MonteCarlos = zeros(1, spots);

greek = 'vanna';
item = 1;
%vanna-1, vega-2 delta-3, gamma-4
for i = 1:spots
    Analytics(i) = AnalyticVanillaCall(1,InitialSpots(i), 0.1,0.05,100,greek);
    vibArray = Vibrato2ndOrder(0.05, 0.1, 1, InitialSpots(i), 100, 100000, 10, 100);
    MonteCarlos(i) = vibArray(1);
end

plot(InitialSpots, Analytics);
hold on;
plot(InitialSpots, MonteCarlos);
legend('Analytic value', 'VMC value');
