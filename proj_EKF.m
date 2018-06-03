function [mu, cov] = proj_EKF(y, mut, covt, Q, R, dt) 

    % EKF predict
    At = A_t(mut, dt); % TODO: replace with A function
    mu = non_lin_dyn(mut, dt); %At*mut; % no control
    cov = At*covt*At' + Q;

    % EKF update
    [yhat, Ct, R] = measure(mu);
    
    if length(y) ~= length(yhat)
        y = y(1:10); %just use GPS measurement 
        yhat = yhat(1:10);
        Ct = Ct(1:10,:);
        R = R(1:10, 1:10);
    end
    
    Kt = (cov*Ct')/(Ct*cov*Ct' + R);
    mu = mu + Kt*(y - yhat);
    cov = cov - Kt*Ct*cov;
end


function x_tp1 = non_lin_dyn(x, dt)
x1=x(1); x2=x(2); x3=x(3); x4=x(4); x5=x(5); 
x6=x(6); x7=x(7); x8=x(8); x9=x(9); x10=x(10); 
% TODO, make this not hard-coded
p = (6378.137+586)*(1-0.01234567^2); %semi-latus rectum of cheif

xdot = [x4;
    x5;
    x6;
    x1*x10^2 *(1 + 2*x7/p) + 2*x10 *(x5-x2*x8/x7);
    -2*x10 *(x4-x1*x8/x7) + x2*x10^2 *(1-x7/p); 
    -x7*x10^2 *x3/p;
    x8;
    x7*x10^2 *(1-x7/p);
    x10;
    -2*x8*x10/x7];

x_tp1 = x + xdot.*dt;

end

function A = A_t(x, dt)
x1=x(1); x2=x(2); x3=x(3); x4=x(4); x5=x(5); 
x6=x(6); x7=x(7); x8=x(8); x9=x(9); x10=x(10); 

% TODO, make this not hard-coded
p = (6378.137+586)*(1-0.01234567^2); %semi-latus rectum of cheif

A = [
[  1,  0,  0,                            dt*x10^2*((2*x7)/p + 1),                                (2*dt*x8*x10)/x7,                   0,  0,                       0,  0,                  0];
[  0,  1,  0,                                  -(2*dt*x8*x10)/x7,                            -dt*x10^2*(x7/p - 1),                   0,  0,                       0,  0,                  0];
[  0,  0,  1,                                                  0,                                               0,    -(dt*x7*x10^2)/p,  0,                       0,  0,                  0];
[ dt,  0,  0,                                                  1,                                       -2*dt*x10,                   0,  0,                       0,  0,                  0];
[  0, dt,  0,                                           2*dt*x10,                                               1,                   0,  0,                       0,  0,                  0];
[  0,  0, dt,                                                  0,                                               0,                   1,  0,                       0,  0,                  0];
[  0,  0,  0,           dt*((2*x1*x10^2)/p + (2*x2*x8*x10)/x7^2),         -dt*((x2*x10^2)/p + (2*x1*x8*x10)/x7^2),    -(dt*x3*x10^2)/p,  1, (dt*x10^2*(p - 2*x7))/p,  0, (2*dt*x8*x10)/x7^2];
[  0,  0,  0,                                  -(2*dt*x2*x10)/x7,                                (2*dt*x1*x10)/x7,                   0, dt,                       1,  0,     -(2*dt*x10)/x7];
[  0,  0,  0,                                                  0,                                               0,                   0,  0,                       0,  1,                  0];
[  0,  0,  0, dt*(2*x5 - (2*x2*x8)/x7 + 2*x1*x10*((2*x7)/p + 1)), -dt*(2*x4 - (2*x1*x8)/x7 + 2*x2*x10*(x7/p - 1)), -(2*dt*x3*x7*x10)/p,  0, -2*dt*x7*x10*(x7/p - 1), dt,   1 - (2*dt*x8)/x7];
];
end













%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% These are replicated in the measure function

function y = g(x)
 
 
ycam = x(1:2)/x(3);
yrange = norm(x(1:3));
 
y = [x];%;ycam;yrange];
 
end
 
function C = C_t(x)
 
C1 = eye(10); %GPS
 
C2 = [1/x(3), 0, -x(1)/x(3)^2, zeros(1,7); %camera
      0, 1/x(3), -x(2)/x(3)^2, zeros(1,7)];    
  
C3 = [x(1), x(2), x(3), zeros(1,7)]./norm(x(1:3)); %range
 
C = [C1];%;C2;C3];
end



