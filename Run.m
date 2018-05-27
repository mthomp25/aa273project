% AA 273, Spring 2018
%
% 5/23/18
%
% Final Project
%
clear variables
close all

addpath('functions');

%% Get truth state
dur = 86400*3;% [s] Run for 3 days
dt = 5; % [s]
t = 0:dt:dur;
tsteps = length(t);

x_true = Truth_sim(dur, dt);
% x_true  	truth state vector [10xN] where N is number of time steps
%           [rho;       [km]
%            rhodot;    [km/s]
%            r;         [km]
%            theta;     [rad]
%            rdot;      [km/s]
%            thetadot]  [rad/s]

% TODO: what do we want our noise to be?
nstates = size(x_true,1);
Q = 0.1*eye(nstates);
R = 0.1; % should measurement noise be different for the 
         % different measurements? e.g. because rchief is much
         % larger than rhox

mu0 = [1, 1, 1, 1, 1, 1, 7000, 0, 0, 1]'; % random initial guess
cov0 = 100.*eye(nstates)'; % very uncertain

% Initialize EKF
mu_EKF = zeros(nstates, tsteps);
cov_EKF = zeros(nstates, nstates, tsteps);
mu_EKF(:, 1) = mu0;
cov_EKF(:, :, 1) = cov0;

% set up timers
EKF_time = zeros(1, tsteps-1);


for tstep = 2:tsteps
    t = dt*(tstep-1);
    
    % for now, we have no control
    
    % propogate state to get xt
    xt = x_true(:, tstep);
    
    % get measurement at t (measurement right now is just noisy truth)
    V = sqrtm(R)*randn(nstates,1);
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
plot3(x_true(1,:), x_true(2,:), x_true(3,:), 'k--');
plot3(mu_EKF(1,:), mu_EKF(2,:), mu_EKF(3,:), 'b')
legend('actual', 'EKF')
xlabel('x');
ylabel('y');
zlabel('z');
title('Relative positions')

figure; hold on;
plot3(x_true(4,:), x_true(5,:), x_true(6,:), 'k--');
plot3(mu_EKF(4,:), mu_EKF(5,:), mu_EKF(6,:), 'b')
legend('actual', 'EKF')
xlabel('x');
ylabel('y');
zlabel('z');
title('Relative velocities')


disp(['Average EKF loop time: ' num2str(mean(EKF_time)) ' sec'])
