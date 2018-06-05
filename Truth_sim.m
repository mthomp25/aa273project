% AA 273, Spring 2018
%
% 5/23/18
%
% Final Project
%
% clear variables
% close all

function x_true = Truth_sim(dur, dt, plotFlag)
% Compute true relative state between 2 spacecraft from predefined orbits.
% The dynamics include 2-body gravity, J2 oblateness effect, and air drag
%
% Inputs
% ------
% dur       sim duration time [s]
% dt        time step [s]
% plotFlag  Option to generate plot of orbit
%           0 (default) - no plots
%           1 - plot
%
% Outputs
% -------
% x_true  	truth state vector [10xN] where N is number of time steps
%           [rho;       [km]
%            rhodot;    [km/s]
%            r;         [km]
%            theta;     [rad]
%            rdot;      [km/s]
%            thetadot]  [rad/s]
%

if nargin < 3
    plotFlag = 0;
end
%% Simulation time
tspan = 0:dt:dur; 

%% Constants
RE = 6378.137;          % [km]
mu = 398600.4415;       % [km^3/s^2] (Montenbruck)
J2 = 1.08263e-3;

m1 = 500;       % [kg]
m2 = 500;       % [kg]
CD = 2.3;
A1 = 5e-6;     % [km^2]
A2 = 7e-6;      % [km^2]
rho0 = 1.225e9; % [kg/km^3]
h0 = 0;         % [km]
H = 10;         % [km]
B1 = CD*A1/m1; % Ballistic coef of spacecraft 1
B2 = CD*A2/m2; % Ballistic coef of spacecraft 2

%% Initial conditions of spacecraft
% --------------- Spacecraft 1 -------------------
% semi-major axis (LEO)
a1 = RE + 586;              % [km]
% eccentricity
e1 = 0.01234567;
% inclination 
i1 = deg2rad(63);           % [rad]
% Right ascension of the ascending node
O1 = deg2rad(4);            % [rad]
% argument of perigee
w1 = deg2rad(125);          % [rad]
% true anomaly
f1 = deg2rad(0);            % [rad]

% Convert to ECI position and velocity 
[r0_1, v0_1] = oe2eci(a1, e1, i1, O1, w1, f1, mu);

% --------------- Spacecraft 2 -------------------
% semi-major axis (slightly lower, will catch up to S/C 1)
a2 = a1-0.002;                % [km]
% eccentricity
% e2 = e1*1.05;
e2 = e1;
% inclination (slight offset to have cross-track separataion)
% i2 = i1 + deg2rad(0.001);   % [rad]
i2 = i1 + deg2rad(0.0001);  
% Right ascension of the ascending node
O2 = O1;                    % [rad]
% argument of perigee
w2 = w1;                    % [rad]
% true anomaly (start it "behind" S/C 1)
f2 = deg2rad(-0.001);         % [rad]

% Convert to ECI position and velocity 
[r0_2, v0_2] = oe2eci(a2, e2, i2, O2, w2, f2, mu);

%% Simulate using numerical integration
options = odeset('RelTol', 1e-9, 'AbsTol', 1e-12);

% Simulate Spacecraft 1
y_init1 = [r0_1, v0_1];
[~, y_out1] = ode113(@(t,y) get_statedot(t, y, mu, rho0, h0, H, B1),...
                    tspan, y_init1, options);
            
% Simulate Spacecraft 2
y_init2 = [r0_2, v0_2];
[~, y_out2] = ode113(@(t,y) get_statedot(t, y, mu, rho0, h0, H, B2), ...
                tspan, y_init2, options);
            
% Extract simulation outputs
r1_eci = y_out1(:,1:3)';
v1_eci = y_out1(:,4:6)';
r2_eci = y_out2(:,1:3)';
v2_eci = y_out2(:,4:6)';

%% Convert to states
% x = [rho, rhodot, r, theta, rdot, thetadot]

N = length(tspan);

% Rotation matrices
R_eci2rtn = ECI2RTN(r1_eci, v1_eci);

rho = zeros(3,N);
rhodot = zeros(3,N);
r = zeros(1,N);
theta = zeros(1,N);
rdot = zeros(1,N);
thetadot = zeros(1,N); 
for i = 1:N
    % Relative position
    rho(:,i) = R_eci2rtn(:,:,i)*(r2_eci(:,i) - r1_eci(:,i));
    
    % Relative velocity
    h = norm(cross(r1_eci(:,i), v1_eci(:,i)));
    omd = h/norm(r1_eci(:,i))^2; 
    
    rhodot(:,i) = R_eci2rtn(:,:,i)*(v2_eci(:,i) - v1_eci(:,i)) - ...
                    cross([0;0;omd], rho(:,i));

    r(:,i) = norm(r1_eci(:,i));
    
    % Compute orbital elements
    [~, ~, ~, ~, w, f] = eci2oe(r1_eci(:,i), v1_eci(:,i), mu);
    theta(:,i) = w + f; % true argument of latitude
    
    % Inertial velocity in the rotating frame
    v_RTN = R_eci2rtn(:,:,i)* v1_eci(:,i); % velocity in rotating frame.
    rdot(:,i) = v_RTN(1);
    
    % theta dot
    thetadot(:,i) = h/norm(r1_eci(:,i))^2; 
end

% Truth state
x_true = [rho; rhodot; r; theta; rdot; thetadot];

%% Plot orbits
if plotFlag
    figure
    grid on; hold on;
    earthPlot
    plot3(r1_eci(1,:), r1_eci(2,:), r1_eci(3,:), '--m', 'LineWidth', 2);
    plot3(r2_eci(1,:), r2_eci(2,:), r2_eci(3,:), '--g', 'LineWidth', 2);
    axis equal;
    xlabel('X [km]'); ylabel('Y [km]'); zlabel('Z [km]');
    legend('S/C 1', 'S/C 2');

%% Plot relative motion
    figure
    grid on; hold on;
    plot3(0,0,0, '*r')
    plot3(rho(2,:), rho(3,:), rho(1,:), 'LineWidth', 2)
    zlabel('\rho_R (m)'); xlabel('\rho_T (m)'); ylabel('\rho_N (m)');
    axis equal; view(3)
    legend('S/C 1', 'S/C 2', 'location', 'best');
    title('Relative Position')
end