% AA 273, Spring 2018
%
% 5/23/18
%
% Final Project
%
clear variables
close all

dt = 5*60; % 5 minute intervals [s]
tspan = 0:dt:(86400*10); % Run for 10 days
tsteps = length(tspan);
% [t, x] = Truth_sim(tspan);
% x = [rhox, rhoy, rhoz, rhodotx, rhodoty, rhodotz, rchief, thchief]

% for testing:
t = tspan;
x = zeros(8, length(tspan));

% TODO: what do we want our noise to be?
Q = 0.1*eye(8);
R = 0.1; % should measurement noise be different for the 
         % different measurements? e.g. because rchief is much
         % larger than rhox

mu0 = [1, 1, 1, 1, 1, 1, 7000, 0]'; % random initial guess
cov0 = 100.*eye(8)'; % very uncertain

% Initialize EKF
mu_EKF = zeros(8, tsteps);
cov_EKF = zeros(8, 8, tsteps);
mu_EKF(:, 1) = mu0;
cov_EKF(:, :, 1) = cov0;

% set up timers
EKF_time = zeros(1, tsteps-1);


for tstep = 2:tsteps
    t = dt*(tstep-1);
    
    % for now, we have no control
    
    % propogate state to get xt
    xt = x(:, tstep);
    
    % get measurement at t (measurement right now is just noisy truth)
    V = sqrtm(R)*randn(8,1);
    y = xt + V;
    
    % EKF
    tic;
    [mu_EKF(:, tstep), cov_EKF(:, :, tstep)] = ...
        proj_EKF(y, mu_EKF(:, tstep-1), cov_EKF(:, :, tstep-1), Q, R);
    EKF_time(tstep-1) = toc;
end


%% Plots
% TODO: plot both satellites' positions relative to Earth
figure; hold on;
plot3(x(1,:), x(2,:), x(3,:), 'k--');
plot3(mu_EKF(1,:), mu_EKF(2,:), mu_EKF(3,:), 'b')
legend('actual', 'EKF')
xlabel('x');
ylabel('y');
zlabel('z');
title('Relative positions')

figure; hold on;
plot3(x(4,:), x(5,:), x(6,:), 'k--');
plot3(mu_EKF(4,:), mu_EKF(5,:), mu_EKF(6,:), 'b')
legend('actual', 'EKF')
xlabel('x');
ylabel('y');
zlabel('z');
title('Relative velocities')


disp(['Average EKF loop time: ' num2str(mean(EKF_time)) ' sec'])

