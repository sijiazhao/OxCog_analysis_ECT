function statlimit = findStatLimit(tw,s,timeaxis)

x0 = s(find(timeaxis<=tw(2) & timeaxis>= tw(1)));

i = 1;
thiscluster=1;
while i<numel(x0)
    if isnan(x0(i))
        i = i+1;
    elseif (x0(i+1) == (x0(i)))
        g = i;
        consec = 1;
        while g<numel(x0) && x0(g+1)==x0(g)
            consec = consec+1;
            g=g+1;
        end
        num(thiscluster) = consec;
        thiscluster=thiscluster+1;
        i=g+1;
    else
        i=i+1;
    end
end
if ~exist('num','var')
    statlimit = 0;
else
    statlimit = max(num);
end
end