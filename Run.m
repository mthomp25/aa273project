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
dur = 86400*0.8;% [s] Run for 1 days

dt = 5; % [s] (this could probably be as large as 30s)
t = 0:dt:dur;
tsteps = length(t);
MEAS_SCHEME = '2cam_switch'; %options: gps_only, 1cam_gps, 2cam_gps, 1cam_switch, 2cam_switch

x_true = Truth_sim(dur, dt);
% x_true  	truth state vector [10xN] where N is number of time steps
%           [rho;       [km]
%            rhodot;    [km/s]
%            r;         [km]
%            theta;     [rad]
%            rdot;      [km/s]
%            thetadot]  [rad/s]

% =========== flip the states =============
rdot = x_true(8,:);
x_true(8,:) = x_true(9,:);
x_true(9,:) = rdot;
x_true(9,:) = wrapTo2Pi(x_true(9,:));
% x_true  	truth state vector [10xN] where N is number of time steps
%           [rho;       [km]
%            rhodot;    [km/s]
%            r;         [km]
%            rdot;      [km/s]
%            theta;     [rad]
%            thetadot]  [rad/s]

nstates = size(x_true,1);
Qdiag = [1e-4, 1e-4, 1e-4, 1e-6, 1e-6, 1e-6, 0.2, 0.0001, 1e-5, 1e-8];
Q = diag(Qdiag.^2);

mu0 = x_true(:,1);
mu0 = mu0 + [(rand(3,1)*4e-3-2e-3); rand(3,1)*2e-5 - 1e-5; zeros(4,1)]; % start with some estimation error
cov0 = diag([5e-3, 5e-3, 5e-3, 0.02e-3, 0.02e-3, 0.02e-3, 1, 0.01e-3, 1e-4, 1e-4]); % taken from Kim et al.

% Initialize EKF
mu_EKF = zeros(nstates, tsteps);
cov_EKF = zeros(nstates, nstates, tsteps);
mu_EKF(:, 1) = mu0;
cov_EKF(:, :, 1) = cov0;

in_view = zeros(2, tsteps); %first row = visnav, second row = gps

for tstep = 2:tsteps
%     t = dt*(tstep-1);
    
    % for now, we have no control
    
    % propogate state to get xt
    xt = x_true(:, tstep);
    
    [y,~,R] = measure(xt); %true measurement
    y = y + sqrtm(R)*randn(length(R(:,1)),1); %add noise
    
    % EKF
    [mu_EKF(:, tstep), cov_EKF(:, :, tstep), in_view(:,tstep)] = ...
        proj_EKF(y, mu_EKF(:, tstep-1), cov_EKF(:, :, tstep-1), Q, dt, MEAS_SCHEME);

end


%% Plots
% TODO: plot both satellites' positions relative to Earth
figure; hold on;
plot3(x_true(1,:), x_true(2,:), x_true(3,:), 'k--');
plot3(mu_EKF(1,2:end), mu_EKF(2,2:end), mu_EKF(3,2:end), 'b')
scatter3(x_true(1,1), x_true(2,1), x_true(3,1), 'o');
legend('actual', 'EKF', 'start')
xlabel('x');
ylabel('y');
zlabel('z');
title('Relative positions')

figure; hold on;
plot3(x_true(4,:), x_true(5,:), x_true(6,:), 'k--');
plot3(mu_EKF(4,2:end), mu_EKF(5,2:end), mu_EKF(6,2:end), 'b')
legend('actual', 'EKF')
xlabel('x');
ylabel('y');
zlabel('z');
title('Relative velocities')

figure % plot relative position
subplot(2,2,1)
hold on; 
plot(x_true(2,:)*1000, x_true(1,:)*1000, 'k--', 'LineWidth', 2)
plot(mu_EKF(2,2:end)*1000, mu_EKF(1,2:end)*1000, 'b', 'LineWidth', 2)
scatter(x_true(2,1)*1000, x_true(1,1)*1000, 'o');
xlabel('\rho_T (m)'); ylabel('\rho_R (m)'); grid on; axis equal;
legend('actual', 'EKF', 'start')
subplot(2,2,2)
hold on; 
plot(x_true(3,:)*1000, x_true(1,:)*1000, 'k--', 'LineWidth', 2)
plot(mu_EKF(3,2:end)*1000, mu_EKF(1,2:end)*1000, 'b', 'LineWidth', 2)
scatter(x_true(3,1)*1000, x_true(1,1)*1000, 'o');
xlabel('\rho_N (m)'); ylabel('\rho_R (m)'); grid on; axis equal;
subplot(2,2,3)
hold on; 
plot(x_true(2,:)*1000, x_true(3,:)*1000, 'k--', 'LineWidth', 2)
plot(mu_EKF(2,2:end)*1000, mu_EKF(3,2:end)*1000, 'b', 'LineWidth', 2)
scatter(x_true(2,1)*1000, x_true(3,1)*1000, 'o');
xlabel('\rho_T (m)'); ylabel('\rho_N (m)'); grid on; axis equal;
subplot(2,2,4)
hold on;
plot3(x_true(2,:)*1000, x_true(3,:)*1000, x_true(1,:)*1000, 'k--', ...
        'LineWidth', 2)
plot3(mu_EKF(2,2:end)*1000, mu_EKF(3,2:end)*1000, mu_EKF(1,2:end)*1000,...
        'b', 'LineWidth', 2)
scatter3(x_true(2,1)*1000, x_true(3,1)*1000, x_true(1,1)*1000, 'o');
grid on;
zlabel('\rho_R (m)'); xlabel('\rho_T (m)'); ylabel('\rho_N (m)');
view(3);
% axis equal;

t = t/3600; % change to hrs
range = vecnorm(x_true(1:3,:),1);
% Plot range
figure
plot(t,range*1000)
ylabel('Separation [m]')
xlabel('Time [hr]')
% TODO: Plot ECI or ECEF?

%% Error plots
figure;  % plot relative position error (truth - EKF)
ax1=subplot(4,1,1);
hold on;
plot(t(2:end), x_true(1,2:end) - mu_EKF(1,2:end));
% TODO: Plot covariance?
plot(t(2:end), 3*sqrt(squeeze(cov_EKF(1,1,2:end))), '--r', 'LineWidth', 1);
plot(t(2:end), -3*sqrt(squeeze(cov_EKF(1,1,2:end))), '--r', 'LineWidth', 1);
grid on;
ylabel('\rho_R error [km]'); %ylim([-1e-3 1e-3])
title('Relative position error')
ax2=subplot(4,1,2);
hold on;
plot(t(2:end), x_true(2,2:end) - mu_EKF(2,2:end));
plot(t(2:end), 3*sqrt(squeeze(cov_EKF(2,2,2:end))), '--r', 'LineWidth', 1);
plot(t(2:end), -3*sqrt(squeeze(cov_EKF(2,2,2:end))), '--r', 'LineWidth', 1);
grid on;
ylabel('\rho_T error [km]'); %ylim([-1e-3 1e-3])
ax3=subplot(4,1,3);
hold on;
plot(t(2:end), x_true(3,2:end) - mu_EKF(3,2:end));
plot(t(2:end), 3*sqrt(squeeze(cov_EKF(3,3,2:end))), '--r', 'LineWidth', 1);
plot(t(2:end), -3*sqrt(squeeze(cov_EKF(3,3,2:end))), '--r', 'LineWidth', 1);
grid on;
ylabel('\rho_N error [km]'); %ylim([-1e-3 1e-3])

ax4=subplot(4,1,4);
plot(t, in_view(1,:), 'k*', t, in_view(2,:), 'r.');
ylim([-.5 1.5])
grid on;
ylabel('Activated');
xlabel('Time [hr]');
legend('VISNAV', 'GPS')
linkaxes([ax1,ax2,ax3,ax4],'x')

figure;  % plot relative velocity error (truth - EKF)
subplot(3,1,1)
plot(t(2:end), x_true(4,2:end) - mu_EKF(4,2:end));
% TODO: Plot covariance?
grid on;
ylabel('v_R error [km/s]');
title('Relative velocity error')
subplot(3,1,2)
plot(t(2:end), x_true(5,2:end) - mu_EKF(5,2:end));
grid on;
ylabel('v_T error [km/s]');
subplot(3,1,3)
plot(t(2:end), x_true(6,2:end) - mu_EKF(6,2:end));
grid on;
ylabel('v_N error [km/s]');
xlabel('Time [hr]')

figure;  % plot relative velocity error (truth - EKF)
subplot(4,1,1)
plot(t(2:end), x_true(7,2:end) - mu_EKF(7,2:end));
% TODO: Plot covariance?
grid on;
ylabel('r error [km]');
title('Absolute position and velocity error')
subplot(4,1,2)
plot(t(2:end), x_true(8,2:end) - mu_EKF(8,2:end));
grid on;
ylabel('\theta error [rad]');
subplot(4,1,3)
plot(t(2:end), x_true(9,2:end) - mu_EKF(9,2:end));
grid on;
ylabel('$\dot{r}$ error [km/s]', 'Interpreter','latex')
subplot(4,1,4)
plot(t(2:end), x_true(10,2:end) - mu_EKF(10,2:end));
grid on;
ylabel('$\dot{\theta}$ error [km/s]', 'Interpreter','latex')
xlabel('Time [hr]')

% Plot covariances
postrace = zeros(1, length(t)-1);
veltrace = zeros(1, length(t)-1);
for ii = 2:length(t)
    postrace(ii-1) = trace(cov_EKF(1:3, 1:3, ii));    
    veltrace(ii-1) = trace(cov_EKF(4:6, 4:6, ii));
end
figure;
subplot(3,1,1)
plot(t(2:end), postrace);
ylabel('Relative Position Covariance Trace');
grid on;
subplot(3,1,2)
plot(t(2:end), veltrace);
ylabel('Relative Velocity Covariance Trace');
grid on;
subplot(3,1,3)
plot(t, in_view(1,:), 'k*', t, in_view(2,:), 'r.');
ylim([-.5 1.5])
grid on;
ylabel('Activated');
xlabel('Time [s]');
legend('VISNAV', 'GPS')
xlabel('Time [s]');

%% animation 
% WARNING: Takes a while to run so uncomment at your own risk
% 
% x = x_true(1,:);
% y = x_true(2,:);
% z = x_true(3,:);
% x_EKF = mu_EKF(1,:);
% y_EKF = mu_EKF(2,:);
% z_EKF = mu_EKF(3,:);
% 
% fig = figure;
% for i = 1:50:length(x)
%     plot3(0,0,0,'*');
%     hold on
%     scatter3(y(i)*1000,z(i)*1000,x(i)*1000, 'o','filled') 
%     scatter3(y_EKF(i)*1000, z_EKF(i)*1000, x_EKF(i)*1000, 'sb')
%     xlabel('\rho_T (m)'); ylabel('\rho_N (m)');  zlabel('\rho_R (m)');
%     grid on; axis equal;
%     legend('Sat 1', 'Sat 2', 'EKF Est', 'location','northwest');
%     hold off
%     xlim([-150 100]); ylim([-50 50]); zlim([-50 50]);
%     
%     drawnow
%     frame = getframe(fig);
%     im{i} = frame2im(frame);
% end
% % Save to .gif
% fname = 'rel_motion.gif';
% for i  = 1:50:length(x)
%     [A, map] = rgb2ind(im{i},256);
%     if i == 1
%         imwrite(A,map,fname,'gif','LoopCount',Inf,'DelayTime',1);
%     else
%         imwrite(A,map,fname,'gif','WriteMode','append','DelayTime',0);
%     end
% end
