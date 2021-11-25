% AMAN PAREKH - 180073
% ME630 - Fall 2021
%% Implicit Euler for 2b
T = 15;

dt = 0.1 ;     % For timestep = 0.1
time = 0:dt:T;
n = T/dt;     % Number of timesteps

phi1 = zeros(ceil(n),1);      % Initializing Solution to 1
phi1_prime = zeros(ceil(n),1);
phi1(1) = 1;           % Assigning Initial Condition
phi1_prime(1) = 0;

for i = 1:n-1
    phi1(i+1) = (phi1(i) + dt*phi1_prime(i))/(1 + dt*dt);
    phi1_prime(i+1) = (phi1_prime(i) - dt*phi1(i))/(1 + dt*dt);
end


t = 0:0.1:T;
plot(t, cos(t), 'r','LineWidth', 2);
hold on;
plot(time(1:length(phi1)), phi1, 'bo-', 'MarkerFaceColor', 'b');
xlabel('x');
xlim([0 10]);
ylabel('\phi');
legend('Exact Solution', 'Implicit Solution for dt=0.1')
grid on;
grid minor;