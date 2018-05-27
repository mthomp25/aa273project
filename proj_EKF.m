function [mu, cov] = proj_EKF(y, mut, covt, Q, R) 

    % EKF predict
    At = eye(10); % TODO: replace with A function
    mu = At*mut; % no control
    cov = At*covt*At' + Q;

    % EKF update
    Ct = eye(10);
    yhat = mu;
    Kt = (cov*Ct')/(Ct*cov*Ct' + R);
    mu = mu + Kt*(y - yhat);
    cov = cov - Kt*Ct*cov;
end