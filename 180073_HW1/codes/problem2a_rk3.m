% AMAN PAREKH - 180073
% ME630 - Fall 2021
%% RK-3 (Low Storage) for 2a
T = 15;

dt = 0.1 ;     % For timestep = 0.1
time = 0:dt:T;
n = T/dt;     % Number of timesteps

phi1 = zeros(ceil(n),1);      % Initializing Solution to 0
phi1(1) = 0;           % Assigning Initial Condition at x=0

phis = 0;
k1 = 0;

for i=1:n-1
    phis = phi1(i);             % Step 1
    k1 = -phi1(i);             % Step 1
    phis = phis + (dt/3)*k1;             % Step 2
    k1 = (-5/9)*k1 - phis;             % Step 2
    phis = phis + (15/16)*dt*k1;             % Step 3
    k1 = (-153/128)*k1 - phis;             % Step 3
    phi1(i+1) = phis + (8/15)*dt*k1;             % Solution Step
end

t = 0:0.1:T;
plot(t, exp(-t), 'r','LineWidth', 2);
hold on;
plot(time(1:length(phi1)), phi1, 'bo-', 'MarkerFaceColor', 'b');
xlabel('x');
xlim([0 10]);
ylabel('\phi');
legend('Exact Solution', 'RK-3 for dt=0.1')
grid on;
grid minor;