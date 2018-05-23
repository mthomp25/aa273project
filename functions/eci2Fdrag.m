function a_drag = eci2Fdrag(r, v, rho0, h0, H, B)

% Compute acceleration due to atmospheric drag given the position and 
% velocity vectors in ECI frame
%
% Inputs
% ------
% r         postition vector
% v         velocity vector
% rho0      density at reference altitude
% h0        reference altitude
% H         characteristic height
% B         ballistic coefficient CD*A/M
%
% Outputs
% -------
% a_drag    Acceleration due to drag

% Earth rotation rate
w_E = [0; 0; 2*pi / 86164]; % [rad/s]

% Compute relative velocity
v_rel = v - cross(w_E, r);

% Compute air density
RE = 6378.137; % assume spherical Earth
h = norm(r) - RE;
rho = rho0 * exp(-(h-h0)/H);

% Compute acceleration due to drag
a_drag = -0.5*B*rho*norm(v_rel).*v_rel;
