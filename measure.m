function [y, C, R] = measure(x)
% x is the 10-dimensional state vector

% y is the measurement (size depends on x, either 10 or 12)
% C is the measurement derivative (of appropriate size)
% R is the measurement covariance (of appropriate size)


y = x; %GPS
C = eye(10); %GPS
Rdiag = [0.0001*ones(3,1); % 10 cm noise on pos measurement?
     	1e-6*ones(3,1); % 1 mm/s noise on vel. Wait, shoud this be sqrt?
        0.001; % 1 m error on absolute position
        1e-3; % 1 cm/s error on absolute velocity
        1e-5; % not sure how big this should be [rad]
        1e-7]; % again, idk [rad/s]

    
% check to see if in view of camera
%if abs(x(1)/x(2)) < tand(70) && abs(x(3)/x(2)) < tand(70) && x(2) < 0

ycam = [x(1)/x(2); x(3)/x(2)];
yrange = norm(x(1:3));

Ccam = [1/x(2), -x(1)/x(2)^2, 0, zeros(1,7); %camera
    0, -x(3)/x(2)^2, 1/x(2), zeros(1,7)];
Crange = [x(1), x(2), x(3), zeros(1,7)]./yrange; %range

y = [y; ycam; yrange];
C = [C; Ccam; Crange];
Rdiag = [Rdiag; 8.7266e-06*ones(2,1).*ycam*yrange; 1e-5];
%end

R = diag(Rdiag.^2);

end

