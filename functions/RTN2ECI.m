function dcm = RTN2ECI(R, V)
% RTN2ECI Computes rotation matrices from RTN frame to inertial frame
%
% Inputs:
%   R - 3xn inertial position vectors
%   V - 3xn inertial velocity vectors
%
% Outputs:
%   dcm - 3x3xn rotation matrices from inertial to RTN

% Set up the input matrices as a (3 x n) series of column vectors
if size(R, 2) == 3 && size(R, 1) ~= 3
    R = R.';
end

if size(V, 2) == 3 && size(V, 1) ~= 3
    V = V.';
end

% Check input sizes
numPos = size(R, 2);
numVel = size(V, 2);

if size(R, 1) ~= 3
    error('Position vectors must come in the form of 3xn column vectors');
elseif size(V, 1) ~= 3
    error('Velocity vectors must come in the form of 3xn column vectors');
elseif numPos ~= numVel
    error('Must provide the same number of position and velocity vectors');
end

dcm = zeros(3, 3, numPos);

for i = 1:numPos
    rVec = R(:, i)./norm(R(:, i));
    
    H = cross(R(:, i), V(:, i));
    nVec = H./norm(H);
    
    T = cross(nVec, rVec);
    tVec = T./norm(T);
    
    dcm(:, :, i) = [rVec, tVec, nVec];
end

end

