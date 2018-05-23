function [r, v] = oe2eci(a, e, i, O, w, nu, mu) 

% Compute the position and velocity vectors in ECI frame given orbital
% elements
%
% Inputs
% ------
% a         semi-major axis
% e         eccentricity
% i         inclination [rad]
% O       	right ascension of ascending node [rad]
% w         argument of pariapsis [rad]
% nu        true anomaly [rad]
% mu        gravitation parameter
%
% Outputs
% -------
% r         postition vector
% v         velocity vector

% Compute E
E = nu2E(nu, e);

% Compute n, the mean motion
n = sqrt(mu/a^3);

% Compute position vector in Perifocal coordinates
r_PQW = [a*(cos(E) - e); a*sqrt(1-e^2)*sin(E); 0];

% Compute velocity vector in Perifocal coordinates
v_PQW = a*n/(1 - e*cos(E)) * [-sin(E); sqrt(1-e^2)*cos(E); 0];

% Compute rotation matrix from PQW to ECI (3-1-3 Euler rots)
Rz_O = [cos(-O)     sin(-O)     0;
        -sin(-O)    cos(-O)     0;
        0           0           1;];
Rx_i = [1   0           0;
        0   cos(-i)     sin(-i);
        0   -sin(-i)    cos(-i);];
Rz_w = [cos(-w)     sin(-w)     0;
        -sin(-w)    cos(-w)     0;
        0           0           1;];
R_PWQ2ECI = Rz_O * Rx_i * Rz_w;

% Compute outputs
r = R_PWQ2ECI * r_PQW;
v = R_PWQ2ECI * v_PQW;
