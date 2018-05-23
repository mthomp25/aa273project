function E = nu2E(nu, e)

% Compute eccentric anomaly from true anomaly
%
% Inputs
% ------
% nu        true anomaly [rad]
% e         eccentricity
%
% Outputs
% -------
% E         Eccentric anomaly [rad]

sinE = (sin(nu)*sqrt(1-e^2)) / (1+e*cos(nu));
cosE = (e + cos(nu)) / (1 + e*cos(nu));
E = wrapTo2Pi(atan2(sinE, cosE));