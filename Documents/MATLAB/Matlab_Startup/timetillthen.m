function timetillthen(date)
% Calculates the number of hours until the specified date, from current
% time

% DO NOT USE ON SUNDAYS

if nargin<1
    error('Please input date or day');
end

current_date = now;
future_date = current_date;
if ischar(date)
    if strcmp('nextmon',date)
        future_date = future_date + (8-(weekday(now)-1));
    elseif strcmp('tom',date)
        future_date = future_date + 1;
    end
    
    future_date = datevec(future_date);
    future_date(4)=10; future_date(5)=0; future_date(6)=0;
    future_date = datenum(future_date);
    timebetween = etime(datevec(future_date),datevec(current_date));
    disp(['Time until ',datestr(datevec(future_date)),' is ',...
          num2str(timebetween/3600),' hours'                    ]   )    
end
