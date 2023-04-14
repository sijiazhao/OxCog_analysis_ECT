function t1 = convertTime(t)

% Convert the pavlovia date '2020-06-07_19h47.34.042' to '2020-06-07
% 19:47:34' and then to something MATLAB readable
t1 = t;
t1 = strrep(t1,'_',' ');
t1 = strrep(t1,'h',':');
t1 = strrep(t1,'.',':');
t1 = t1(1:19);

t1 = datevec(datenum(t1));