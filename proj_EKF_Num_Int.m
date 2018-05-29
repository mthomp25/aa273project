function [mu, cov] = proj_EKF_Num_Int(y, mut, covt, Q, R, dt) 
    muEarth = 398600.4415;       % [km^3/s^2] (Montenbruck)

    % EKF predict
    At = A_t(mut, dt); % TODO: replace with A function
    options = odeset('RelTol', 1e-4, 'AbsTol', 1e-6);
    [~, yout] = ode113(@(t,y) get_relstatedot(t, y, muEarth), [0 dt], ...
                mut, options);
    mu = yout(end,:)';
    cov = At*covt*At' + Q;

    % EKF update
    Ct = eye(10);
    yhat = mu;
    Kt = (cov*Ct')/(Ct*cov*Ct' + R);
    mu = mu + Kt*(y - yhat);
    cov = cov - Kt*Ct*cov;
end


function x_tp1 = lin_dyn(x, dt)
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


