function MJD = Cal2MJD(MM, DD, YYYY)

% Convert calendar representation into Modified Julian Date
%
% Inputs
% ------
% MM        month number [1-12]
% DD        date, including hours, minutes, and seconds
% YYYY      year
%
% Outputs
% -------
% MJD       Modified Julian Date
%% Check if inputs are valid
switch MM
    case 1
        if DD > 31
            error('Invalid Date Input')
        end
    case 2
        if mod(YYYY, 4) > 0
            % non leap year
            if DD > 28
                error('Invalid Date Input')
            end
        else
            % leap year
            if DD > 29
                error('Invalid Date Input')
            end
        end
    case 3
        if DD > 31
            error('Invalid Date Input')
        end
    case 4
        if DD > 30
            error('Invalid Date Input')
        end
    case 5
        if DD > 31
            error('Invalid Date Input')
        end
    case 6
        if DD > 30
            error('Invalid Date Input')
        end
    case 7
        if DD > 31
            error('Invalid Date Input')
        end
    case 8
        if DD > 31
            error('Invalid Date Input')
        end
    case 9
        if DD > 30
            error('Invalid Date Input')
        end
    case 10
        if DD > 31
            error('Invalid Date Input')
        end
    case 11
        if DD > 30
            error('Invalid Date Input')
        end
    case 12
        if DD > 31
            error('Invalid Date Input')
        end
    otherwise
        error('Invalid Month Input')
end

%% Compute MJD from Montenbruck Appendix A.1.1
if MM <= 2
    y = YYYY - 1;
    m = MM + 12;
else
    y = YYYY;
    m = MM;
end

% Assume we won't deal with time before Oct 10, 1582
B = floor(y/400) - floor(y/100) + floor(y/4);

MJD = 365*y - 679004 + B + floor(30.6001*(m+1)) + DD;
