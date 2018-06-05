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
statedot(7) = state(8); % r_dot
statedot(9) = state(10); % theta_dot

p = (6378.137+586)*(1-0.01234567^2); %semi-latus rectum of cheif

x1=state(1); x2=state(2); x3=state(3); x4=state(4); x5=state(5); 
x6=state(6); x7=state(7); x8=state(8); x9=state(9); x10=state(10);

statedot = [x4;
    x5;
    x6;
    x1*x10^2 *(1 + 2*x7/p) + 2*x10 *(x5-x2*x8/x7);
    -2*x10 *(x4-x1*x8/x7) + x2*x10^2 *(1-x7/p); 
    -x7*x10^2 *x3/p;
    x8;
    x7*x10^2 *(1-x7/p);
    x10;
    -2*x8*x10/x7];

% statedot(8) = r0 * state(10)^2 - mu/r0^2; % r_ddot
% statedot(10) = -2*state(8)*state(10)/r0; % theta_ddot
% 
% % Alfriend (4.14) - (4.16)
% statedot(4) = -mu*(r0 + x)/num + mu/r0^2 + 2*theta_dot*ydot + ...
%               statedot(10)*y + theta_dot^2*x ; % x_ddot
% statedot(5) = -mu*y/num - 2*theta_dot*xdot - statedot(10)*x + ...
%               theta_dot^2*y; % y_ddot
% statedot(6) = -mu*z/num; % z_ddot
end