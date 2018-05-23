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

% J2 acceleration


statedot(1:3) = state(4:6); % velocities
statedot(4:6) = -mu*state(1:3) ./ norm(state(1:3))^3; % FODE

end