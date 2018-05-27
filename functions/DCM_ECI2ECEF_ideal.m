function R_ECI2ECEF = DCM_ECI2ECEF_ideal(GMST)

% Compute rotation matrix from ECI to ECEF reference frame using the ideal
% configuration (neglect nutation, precession, polar motion, equinoxes)
%
% Inputs
% ------
% GMST       UT1 time expressed in MJD format
%
% Outputs
% -------
% R_ECI2ECEF    Rotation matrix to convert from ECI to ECEF

% Rotation about z-axis
R_ECI2ECEF = [cos(GMST)     sin(GMST)     0;
              -sin(GMST)    cos(GMST)     0;
              0             0           1;];