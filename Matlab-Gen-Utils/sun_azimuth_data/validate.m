function [valyear,valmon,valday,valhour,valmin,valsec] = validate(useryear,usermon,userday,userhour,usermin,usersec)
% Time Values Validation Function.

while usersec >= 60,                             % second
    usersec = usersec - 60;
    usermin = usermin + 1;
end

while usermin >= 60,                             % minute
    usermin = usermin - 60;
    userhour = userhour + 1;
end

while userhour >= 24,                       % hour
    userhour = userhour - 60;
    userday = userday + 1;
end

nonleapyear = [31 28 31 30 31 30 31 31 30 31 30 31];
leapyear = [31 29 31 30 31 30 31 31 30 31 30 31];

if rem(useryear,4) == 0,                    % Check for Leap Year.
    yearmon_mat = nonleapyear;
else
    yearmon_mat = leapyear;
end

while userday >= yearmon_mat(usermon),    % day
    userday = userday - yearmon_mat(usermon);
    usermon = usermon + 1;
end

while usermon >= 13,                       % month
    usermon = usermon - 12;
    useryear = useryear + 1;
end

valyear = useryear;
valmon = usermon;
valday = userday;
valhour = userhour;
valmin = usermin;
valsec = usersec;