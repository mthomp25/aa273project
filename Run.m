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
dur = 86400*1;% [s] Run for 1 days
dt = 5; % [s] (this could probably be as large as 30s)
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
R = [0.0001*eye(3) zeros(3,7); % 10 cm noise on pos measurement?
     zeros(3) 1e-6*eye(3) zeros(3,4); % 1 mm/s noise on vel. Wait, shoud this be sqrt?
     zeros(1,6) 0.001 0 0 0; % 1 m error on absolute position
     zeros(1,7) 1e-5 0 0; % not sure how big this should be [rad]
     zeros(1,8) 1e-5 0; % 1 cm/s error on absolute velocity
     zeros(1,9) 1e-7]; % again, idk [rad/s]

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

figure % plot relative position
subplot(2,2,1)
hold on; 
plot(x_true(2,:)*1000, x_true(1,:)*1000, 'k--', 'LineWidth', 2)
plot(mu_EKF(2,:)*1000, mu_EKF(1,:)*1000, 'b', 'LineWidth', 2)
xlabel('\rho_T (m)'); ylabel('\rho_R (m)'); grid on; axis equal;
legend('actual', 'EKF')
subplot(2,2,2)
hold on; 
plot(x_true(3,:)*1000, x_true(1,:)*1000, 'k--', 'LineWidth', 2)
plot(mu_EKF(3,:)*1000, mu_EKF(1,:)*1000, 'b', 'LineWidth', 2)
xlabel('\rho_N (m)'); ylabel('\rho_R (m)'); grid on; axis equal;
subplot(2,2,3)
hold on; 
plot(x_true(2,:)*1000, x_true(3,:)*1000, 'k--', 'LineWidth', 2)
plot(mu_EKF(2,:)*1000, mu_EKF(3,:)*1000, 'b', 'LineWidth', 2)
xlabel('\rho_T (m)'); ylabel('\rho_N (m)'); grid on; axis equal;
subplot(2,2,4)
plot3(x_true(2,:)*1000, x_true(3,:)*1000, x_true(1,:)*1000, 'k--', ...
        'LineWidth', 2)
plot3(mu_EKF(2,:)*1000, mu_EKF(3,:)*1000, mu_EKF(1,:)*1000, 'b', ...
        'LineWidth', 2)
grid on;
zlabel('\rho_R (m)'); xlabel('\rho_T (m)'); ylabel('\rho_N (m)');
view(3);
axis equal;
