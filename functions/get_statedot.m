function statedot = get_statedot(t, state, mu, rho0, h0, H, B)

% Compute state derivative vector given state and time. This is used by
% MATLAB integration functions.
%
% Inputs
% ------
% t         time
% state     vector [rx ry rz vx vy vz] expressed in ECI frame
% mu        gravitational parameter in desired unit
% rho0      density at reference altitude
% h0        reference altitude
% H         characteristic height
% B         ballistic coefficient CD*A/M
%
% Outputs
% -------
% statedot  Time derivative of the state vector

% Preallocate
statedot = zeros(6,1);

% Call drag force function
a_drag = eci2Fdrag(state(1:3), state(4:6), rho0, h0, H, B);

% Solar Radiation Pressure
% Should we add this? 3 orders of magnitude smaller than J2 effects and
% requires keeping track of time and position of the sun.

% J2 acceleration (Alfriend (4.93) with spherical grav removed)
r_norm = norm(state(1:3));
J2 = 1.08263e-3;
RE = 6378.137;
a_J2 = [state(1)*(5*state(3)^2/r_norm^2 - 1);
        state(2)*(5*state(3)^2/r_norm^2 - 1);
        state(3)*(5*state(3)^2/r_norm^2 - 3)] ...
        * 3/2 * J2 * (RE/r_norm)^2 * mu/r_norm^3;

% Kepler 
statedot(1:3) = state(4:6); % velocities
statedot(4:6) = -mu*state(1:3) ./ norm(state(1:3))^3; % FODE

% Output
statedot(4:6) = statedot(4:6) + a_drag + a_J2;

end