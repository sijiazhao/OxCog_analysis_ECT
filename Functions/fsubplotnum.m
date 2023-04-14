function [x,y]=fsubplotnum(n)
if n<=4
    x=n; y=1;
elseif n<=6
    x = 4; y =2;
elseif n <= 8
    x = 4; y = 2;
elseif n==9
    x = 3; y = 3;
else
    x=ceil(n/3);
    y=3;
end
    
end
