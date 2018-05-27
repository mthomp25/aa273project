function r_ECEF = eci2ecef(t_sim, r_ECI, DD, MM, YYYY, dUT1)

% Transform vector in perifocal to ECI
% elements
%
% Inputs
% ------
% t_sim     simulation time [s]
% r_ECI     position vector in ECI frame
% DD       	UTC day
% MM        month
% YYYY      year
% dUT1      Time correction factor UT1 - UTC [msec]
%
% Outputs
% -------
% r_ECI     position in ECEF coordinates

% Compute MJD
MJD = Cal2MJD(MM, DD, YYYY);
% Add simulation time
MJD = MJD + t_sim / 86400 + dUT1/1000/86400; % convert to UT1

% Compute GMST
GMST = UT12GMST(MJD); % [rad]

% Get rotation matrix from ECI to ECEF
R_ECI2ECEF = DCM_ECI2ECEF_ideal(GMST);

r_ECEF = R_ECI2ECEF * r_ECI;