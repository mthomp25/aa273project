function GMST = UT12GMST(UT1)

% Convert UT1 time in MJD format to Greenwich Mean Sidereal Time
%
% Inputs
% ------
% UT1       UT1 time expressed in MJD format
%
% Outputs
% -------
% GMST      Greenwich Mean Sidereal Time % [rad] [0, 2pi)

% Compute time since Janyary 1, 2000 at 12:00
d = UT1 - 51544.5;

% Compute GMST
GMST = (280.4606 + 360.9856473 * d) * pi/180; % [rad]
% Wrap angle to 2pi
GMST = wrapTo2Pi(GMST);  % [rad]