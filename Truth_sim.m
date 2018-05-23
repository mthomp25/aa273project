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

%% Initial conditions of satellites

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
[r0_1, v0_1] = oe2eci(a1, e1, i1, O1, w1, f1);

% --------------- Spacecraft 2 -------------------


%% Simulate using numerical integration

