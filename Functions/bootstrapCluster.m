function [s_boot_corrected,clusters] = bootstrapCluster(s_boot,statlimit,timeaxis)
clusters = [];
numclu = 0;
consec = 1;

i = 1;
while i<numel(s_boot)
    if isnan(s_boot(i))
        i = i+1;
    elseif ( s_boot(i+1) == (s_boot(i)) )
        numclu = numclu + 1;
        g = i;
        consec = 1;
        clusters(numclu,1) = timeaxis(g);
        while g<numel(s_boot) && s_boot(g+1)==s_boot(g)
            consec = consec+1;
            g=g+1;
        end
        if consec<=statlimit
            s_boot(1,i:g)=NaN;
            clusters(numclu,2) = NaN;
            i=g+1;
        else
            i=g+1;
            clusters(numclu,2) = timeaxis(g);
        end
    else
        i=i+1;
    end
end
s_boot_corrected=s_boot;
