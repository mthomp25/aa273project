% AA 273, Spring 2018
%
% 5/23/18
%
% Final Project
%
clear variables
close all

addpath('functions');
%% Constants
RE = 6378.137;          % [km]
mu = 398600.4415;       % [km^3/s^2] (Montenbruck)
J2 = 1.08263e-3;

m1 = 500;       % [kg]
m2 = 300;       % [kg]
CD = 2.3;
A1 = 20e-6;      % [km^2]
A2 = 5e-6;      % [km^2]
rho0 = 1.225e9; % [kg/km^3]
h0 = 0;         % [km]
H = 10;         % [km]
B1 = CD*A1/m1; % Ballistic coef of spacecraft 1
B2 = CD*A2/m2; % Ballistic coef of spacecraft 2

%% Initial conditions of spacecraft
% --------------- Spacecraft 1 -------------------
% semi-major axis
a1 = RE + 586; % [km]
% eccentricity
e1 = 0.01234567;
% inclination 
i1 = deg2rad(63); % [rad]
% Right ascension of the ascending node
O1 = deg2rad(4); % [rad]
% argument of perigee
w1 = deg2rad(125); % [rad]
% true anomaly
f1 = deg2rad(0); % [rad]

% Convert to ECI position and velocity 
[r0_1, v0_1] = oe2eci(a1, e1, i1, O1, w1, f1, mu);

% --------------- Spacecraft 2 -------------------
% semi-major axis
a2 = a1; % [km]
% eccentricity
e2 = e1;
% inclination 
i2 = i1 + deg2rad(10); % [rad]
% Right ascension of the ascending node
O2 = O1; % [rad]
% argument of perigee
w2 = w1; % [rad]
% true anomaly
f2 = deg2rad(2); % [rad]

% Convert to ECI position and velocity 
[r0_2, v0_2] = oe2eci(a2, e2, i2, O2, w2, f2, mu);

%% Simulate using numerical integration
dt = 5; % [s]
tspan = 0:5:(86400*10); % Run for 10 days
options = odeset('RelTol', 1e-6, 'AbsTol', 1e-9);

% Simulate Spacecraft 1
y_init1 = [r0_1, v0_1];
[t_out, y_out1] = ode113(@(t,y) get_statedot(t, y, mu, rho0, h0, H, B1),...
                    tspan, y_init1, options);
            
% Simulate Spacecraft 2
y_init2 = [r0_1, v0_1];
[~, y_out2] = ode113(@(t,y) get_statedot(t, y, mu, rho0, h0, H, B2), ...
                tspan, y_init1, options);
            
% Extract simulation outputs
% r1_eci =
            
%% Plot orbits
figure
grid on; hold on;
earthPlot

axis equal;
