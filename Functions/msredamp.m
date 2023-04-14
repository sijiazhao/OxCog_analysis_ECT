function x = msredamp(MS,Xn,method)

x = MS;
switch method
    case 'zeta'
        x.damp(MS.zeta<=1) = exp(-MS.zeta(MS.zeta<=1)*pi./sqrt(1-MS.zeta(MS.zeta<=1).^2));
        x.damp(MS.zeta>1)=0;
    case 'rawdata'
        for k = 1:length(MS.onset)
            x0 = Xn(MS.onset(k));
            fn = Xn(floor((MS.overtime(k)+MS.endtime(k))/2));
            ov = Xn(floor((MS.peaktime(k)+MS.overtime(k))/2));
            if fn-x0~=0
                x.damp(k) = (ov-fn)/(fn-x0);
            else
                x.damp(k)=0;
            end
            if x.damp(k)>0
                x.zeta(k) = 1./sqrt(1+(pi./log(x.damp(k))).^2);
            elseif x.damp(k)<=0 && ((ov-fn)>0)*((ov-x0)>0)==0
                x.zeta(k)=1;
            elseif x.damp(k)<=0 && ((ov-fn)>0)*((ov-x0)>0)==1
                x.zeta(k)=0;
            end
        end
end