function [a, e, i, O, w, nu] = eci2oe(r, v, mu)

% Compute orbital elements given the position and velocity vectors in ECI 
% frame
%
% Inputs
% ------
% r         postition vector
% v         velocity vector
% mu        gravitation parameter
%
% Outputs
% -------
% a         semi-major axis
% e         eccentricity
% i         inclination [rad]
% O       	right ascension of ascending node [rad]
% w         argument of pariapsis [rad]
% nu        true anomaly [rad]

% Compute eccentricity vector
e_vec = 1/mu * ((norm(v)^2 - mu/norm(r)).*r - dot(r,v) .* v );
e = norm(e_vec);


% Compute W_hat direction
h = cross(r, v);
W = h/norm(h);

% Compute inclination
i = atan2(sqrt(W(1)^2 + W(2)^2), W(3));

% Compute right ascension of ascending node
if i == 0
    O = NaN; % This means equatorial orbit
else
    O = atan2(W(1), -W(2));
end

% Compute specific mechanical energy
Eps = norm(v)^2/2 - mu/norm(r);

% Compute semi-major axis
a = -mu/2/Eps;

% Compute semi-latus rectum
p = norm(h)^2 / mu;

% Compute true anomaly
if e == 0
    nu = NaN;
else
    n = sqrt(mu/a^3);
    E = atan2(dot(r,v)/(a^2*n), 1 - norm(r)/a);
    nu = atan2(sin(E)*sqrt(1-e^2)/(1 - e*cos(E)), (cos(E) - e)/(1-e*cos(E)));
    nu = wrapTo2Pi(nu);
end

% Compute argument of latitude
u = atan2(r(3)/sin(i), r(1)*cos(O) + r(2)*sin(O));

% Compute argument of perigee
w = wrapTo2Pi(u - nu);
