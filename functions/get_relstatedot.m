function statedot = get_relstatedot(t, state, mu)
% Compute state derivative vector given state and time. This is used by
% MATLAB integration functions.
%
% Inputs
% ------
% state     relative vector [rx ry rz vx vy vz r theta rdot thetadot]
% mu        gravitational parameter in desired unit
%
% Outputs
% -------
% statedot  Time derivative of the state vector

r0 = state(7);
theta_dot = state(10);
x = state(1);
y = state(2);
z = state(3);
xdot = state(4);
ydot = state(5);
num = ((r0+x)^2 + y^2 + z^2)^(3/2);

statedot = zeros(10,1);
statedot(1:3) = state(4:6); % velocities
statedot(7:8) = state(9:10); % r, theta
statedot(9) = r0 * state(10)^2 - mu/r0^2; % r_ddot
statedot(10) = -2*state(9)*state(10)/r0; % theta_ddot

% Alfriend (4.14) - (4.16)
statedot(4) = -mu*(r0 + x)/num + mu/r0^2 + 2*theta_dot*ydot + ...
              statedot(10)*y + theta_dot^2*x ; % x_ddot
statedot(5) = -mu*y/num - 2*theta_dot*xdot - statedot(10)*x + ...
              theta_dot^2*y; % y_ddot
statedot(6) = -mu*z/num; % z_ddot
end