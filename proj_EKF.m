function [mu, cov, in_view] = proj_EKF(y, mut, covt, Q, dt, MEAS_SCHEME) 

    % EKF predict
    At = A_t(mut, dt); % TODO: replace with A function
    mu = non_lin_dyn(mut, dt); %At*mut; % no control
    cov = At*covt*At' + Q;

    % EKF update
   [yhat, Ct, R] = measure(mu);
    numGpsMeas = 10;

    if strcmp(MEAS_SCHEME,'1cam_gps') || strcmp(MEAS_SCHEME,'1cam_switch')
        visible = abs(mu(1)/mu(2)) < tand(70) && abs(mu(3)/mu(2)) < tand(70) && mu(2) < 0;
    else
        visible = abs(mu(1)/mu(2)) < tand(70) && abs(mu(3)/mu(2)) < tand(70);
    end
    
    if strcmp(MEAS_SCHEME,'gps_only')
        %if length(y) ~= length(yhat)
        y = y(1:numGpsMeas); %just use GPS measurement 
        yhat = yhat(1:numGpsMeas);
        Ct = Ct(1:numGpsMeas,:);
        R = R(1:numGpsMeas, 1:numGpsMeas);
        
        in_view = [0;1];
    elseif strcmp(MEAS_SCHEME,'1cam_gps') || strcmp(MEAS_SCHEME,'2cam_gps')
        if visible
            in_view = [1;1];
        else
            y = y(1:numGpsMeas); %just use GPS measurement 
            yhat = yhat(1:numGpsMeas);
            Ct = Ct(1:numGpsMeas,:);
            R = R(1:numGpsMeas, 1:numGpsMeas);
        
            in_view = [0;1];
        end
    elseif strcmp(MEAS_SCHEME,'1cam_switch') || strcmp(MEAS_SCHEME,'2cam_switch')
        
        if norm(mu(1:3)) < .03
            if ~visible
                in_view = [0;0];
                return
            end
            y = y(numGpsMeas:end); %only use camera
            yhat = yhat(numGpsMeas:end);
            Ct = Ct(numGpsMeas:end,:);
            R = R(numGpsMeas:end, numGpsMeas:end);
            in_view = [1;0];
        else
            y = y(1:numGpsMeas); %just use GPS measurement
            yhat = yhat(1:numGpsMeas);
            Ct = Ct(1:numGpsMeas,:);
            R = R(1:numGpsMeas, 1:numGpsMeas);
            in_view = [0;1];
        end
    else
        error('MEAS_SCHEME not recognized!!!')
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

