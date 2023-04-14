function newMS = msoptim(MS,Xn,zetalim,varargin)

%
%
% varargin{1} -> Threshold for amplitude
% 
if isempty(varargin)==1 || strcmp(varargin{1},'micro')==1
    maxamp = 240;               % -> used in [4]: threshold used in removal of large saccade  [unit:MOA]
elseif strcmp(varargin{1},'large')==1
    maxamp = inf;               % -> used in [4]: threshold used in removal of large saccade  [unit:MOA]
elseif isnumeric(varargin{1})==1
    maxamp = varargin{1};
end
if isempty(varargin)==1 || length(varargin)==1
    px_to_MOA = 21.6/640*60;    % -> used in [5]: convert pixel to MOA (minute of arc)
elseif length(varargin)==2
    px_to_deg = varargin{2};
    px_to_MOA = px_to_deg*60;
end
newMS = MS;
oldCOD = zeros(1,length(MS.onset));
newMS.COD = zeros(length(MS.onset),1);
for a=1:length(MS.onset)
    ydata = px_to_MOA*Xn(MS.onset(a):MS.endtime(a));
    t = 0:MS.endtime(a)-MS.onset(a);
%     ydata = px_to_MOA*Xn(MS.onset(a):floor((MS.overtime(a)+MS.endtime(a))/2));
%     t = 0:floor((MS.overtime(a)+MS.endtime(a))/2)-MS.onset(a);
    Fparam = [MS.amp(a) MS.zeta(a) MS.wn(a)];
    [Eparam,model] = MSparam(t,ydata,Fparam);
    [sse,~] = model(Eparam);
    newMS.amp(a) = Eparam(1);
    newMS.zeta(a) = Eparam(2);
    newMS.wn(a) = Eparam(3);
    newMS.COD(a) = 1-sse/sum((ydata-mean(ydata)).^2);
    [sse,~] = model(Fparam);
    oldCOD(a) = 1-sse/sum((ydata-mean(ydata)).^2);
end

index = find(abs(newMS.amp)>maxamp | newMS.zeta>zetalim(2) | newMS.zeta<zetalim(1) | newMS.wn==Inf)';
disp(['[msopt] Removed ',num2str(length(index)),' / ',num2str(length(MS.onset)),' trials'])
newMS = removeMS(newMS,index,1);
oldCOD(index)=[];
disp(['[msopt] Coefficient of Determination: ',num2str(mean(oldCOD)),' --> ',num2str(mean(newMS.COD)),' trials'])
newMS.Ka = newMS.wn./(2*newMS.zeta);
% newMS.Tp = pi./newMS.wn./sqrt(1-newMS.zeta.^2);

end

function [estimates, model] = MSparam(t, ydata, Fparam) %???
ydata = ydata-ydata(1);
start_point = Fparam;
model = @expfun;
estimates = fminsearch(model, start_point);
    function [sse, FittedCurve] = expfun(params)
        A = params(1);
        Z = params(2);
        Wn = params(3);
        FittedCurve = A*(1-exp(-Z*Wn*t).*(Z/sqrt(1-Z^2)*sin(sqrt(1-Z^2)*Wn*t)+cos(sqrt(1-Z^2)*Wn*t)));
        ErrorVector = FittedCurve - ydata;
        sse = sum(ErrorVector.^2);
    end
end

function newMS = removeMS(MS,remindex,properties)
    MS.onset(remindex)=[];
    MS.offset(remindex)=[];
    if properties == 1
        MS.duration(remindex)=[];
        MS.amp(remindex)=[];
        MS.peakvel(remindex)=[];
        MS.peaktime(remindex)=[];
        MS.double(remindex)=[];
        MS.overtime(remindex)=[];
        MS.endtime(remindex)=[];
        MS.starttime(remindex)=[];
        MS.overshoot(remindex)=[];
        MS.endvel(remindex)=[];
        MS.zerotime(remindex)=[];
        MS.tremor(remindex)=[];
        MS.damp(remindex)=[];
        MS.zeta(remindex)=[];
        MS.Tp(remindex)=[];
        MS.wn(remindex)=[];
        MS.Ka(remindex)=[];
        MS.preISI(remindex)=[];
        MS.preISPC(remindex)=[];
        MS.postISI(remindex)=[];
        MS.postISPC(remindex)=[];
        MS.COD(remindex)=[];
    end
    newMS = MS;
end