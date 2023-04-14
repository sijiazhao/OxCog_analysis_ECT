function MS = msdetect(Xn, smpfreq, blinkrem, lamda, varargin)

%MSDETECT Microsaccades detection.
%
%   -INPUTS
%       Xn       : time series of EYE POSITION                        <1xN vector> [unit:pixel]
%       smpfreq  : SAMPLING FREQUENCY                                 <scalar>     [unit:Hz]
%       blinkrem : threshold for BLINK detection                      <scalar>     [unit:pixel]
%       lamda    : threshold parameter for MICROSACCADE detection(#)  <scalar>     [unit:none]
%                  -----
%                  (#) MICROSACCADE_detection_threshold = lamda * MS.sigma
%                  -----
%
%   -OUTPUTS
%       MS <1x1 struct>
%       - MS.onset    : ONSET of each microsaccade          <1xK vector> [unit:sample point]
%       - MS.offset   : OFFSET of each microsaccade         <1xK vector> [unit:sample point]
%       - MS.duration : DURATION of each microsaccade       <1xK vector> [unit:msec]
%       - MS.amp      : AMPLITUDE of each microsaccade      <1xK vector> [unit:MOA]
%       - MS.peakvel  : PEAK VELOCITY of each microsaccade  <1xK vector> [unit:deg]
%       - MS.vel      : time series of EYE VELOCITY         <1xN vector> [unit:pixel]
%       - MS.sigma    : STANDARD DEVIATION of eye velocity  <scalar>     [unit:pixel]
%       - MS.xxx...

% [0] Initial Setting
ridnum = floor(100*smpfreq/1000);               % -> used in [1]: remove time zone (100ms) immediately before blink onsets or after offsets [unit:samplepoint]
MAnum = 4;                  % -> used in [2]: number of points used in moving average                           [unit:samplepoint]
Tcorrection = 1;            % -> used in [3]: time correction for microsaccade onset/offset                     [unit:samplepoint]
tooclose = 3;               % -> used in [4]: threshold for removing inappropriate trials [unit:samplepoint]
longDur = 100;
adjaThr = 50;
if isempty(varargin)==1 || strcmp(varargin{1},'micro')==1
    maxamp = 90;               % -> used in [4]: threshold used in removal of large saccade  [unit:MOA]
    maxvel = 1000;               % -> used in [4]: threshold used in removal of large saccade  [unit:deg]
elseif strcmp(varargin{1},'large')==1
    maxamp = inf;               % -> used in [4]: threshold used in removal of large saccade  [unit:MOA]
    maxvel = inf;               % -> used in [4]: threshold used in removal of large saccade  [unit:deg]
elseif strcmp(varargin{1},'EMR9')==1 %???
    ridnum = 10;
    longDur = inf;
    adjaThr = 10;
    maxamp = inf;
    maxvel = inf;
elseif isnumeric(varargin{1})==1
    maxamp = varargin{1}(1);
    maxvel = varargin{1}(2);
end

if isempty(varargin)==1 || length(varargin)==1
    px_to_MOA = 21.6/640*60;    % -> used in [5]: convert pixel to MOA (minute of arc)
    px_to_deg = 21.6/640;       % -> used in [5]: convert pixel to degree
elseif length(varargin)==2
    px_to_deg = varargin{2};
    px_to_MOA = px_to_deg*60;
end
% [1] Trial Definition
%      Remove time zone around blinks
%       OUTPUT
%         -> 'MSdomain': Time domain for microsaccades onset/offset <1xM vector> [unit:sample point]   
%         -> 'Vndomain': Time domain for velocity calculation       <1xK vector> [unit:sample point]
%
eyeline       = [0 diff(Xn<blinkrem)];
% eyeline       = [0 isnan(Xn)];
open_to_close = find(eyeline==1);
close_to_open = find(eyeline==-1);
if isempty(open_to_close)==0 && isempty(close_to_open)==0
    startwith     = close_to_open(1) > open_to_close(1);            % Whole trial starts with eyes 1:open, 0:closed
    endwith       = close_to_open(end) > open_to_close(end);        % Whole trial ends with eyes 1:open, 0:closed
    MSdomain      = [];
    Vndomain      = [];
    minlength     = ridnum*2 + 2*MAnum + smpfreq/10;
    tmpdomain = [[ones(1,startwith) close_to_open]' [open_to_close length(Xn)*ones(1,endwith)]'];
    tmpdomain = tmpdomain(diff(tmpdomain')>minlength,:);
    MSdomain = tmpdomain + ridnum * [ones(size(tmpdomain,1),1) -ones(size(tmpdomain,1),1)];
    Vndomain = tmpdomain + MAnum  * [ones(size(tmpdomain,1),1) -ones(size(tmpdomain,1),1)];
    clear tmpdomain
elseif isempty(open_to_close)==0 && isempty(close_to_open)==1
    MSdomain = [ridnum+1 open_to_close(1)-ridnum];
    Vndomain = [MAnum+1 open_to_close(1)-MAnum];
elseif isempty(open_to_close)==1 && isempty(close_to_open)==0
    MSdomain = [close_to_open(1)+ridnum length(Xn)-ridnum];
    Vndomain = [close_to_open(1)+MAnum length(Xn)-MAnum];
elseif isempty(open_to_close)==1 && isempty(close_to_open)==1
    MSdomain = [ridnum+1 length(Xn)-ridnum];
    Vndomain = [MAnum+1 length(Xn)-MAnum];
end

% [2] Caluculate Eye Velocity
%       OUTPUT
%         -> 'Vn': velocity of eye position <1xM vector> [unit:pixel]   
%
% MAnum/2 ... centre point
% 

Vn = zeros(length(Xn),1);
for a=1:size(Vndomain,1)
    for i=Vndomain(a,1):Vndomain(a,2)
%         Vn(i) = (sum(Xn(i+1:i+MAnum))-sum(Xn(i-MAnum:i-1)))/(1/smpfreq*2*MAnum);
        Vn(i) = (sum(Xn(i:i+MAnum))-sum(Xn(i-MAnum:i)))/MAnum/(MAnum+1)*smpfreq;
    end
end

% [3] Detect Microsaccades (MAIN PART)
%  SELECT DETECTION METHOD BY ILLUSTRATING 'VARARGIN' 
%   (i) 'threshold': detect microsaccades as over-threshold points
%   (ii) 'similarity': detect microsaccades by calculating the similarity
%   to the model [REQUIRE MODEL] (not included here)
%
%         -> 'sigma': standard deviation of Vn <scalar> [unit:pixel]   
%            (#) MICROSACCADE_detection_threshold = lamda * sigma
%
%         -> 'MSindex': Index of microsaccades included in MSdomain (domain of definition)
%
sigma = sqrt(median(nonzeros(Vn).^2) - median(nonzeros(Vn))^2);
overthreshold = (Vn>lamda*sigma+median(nonzeros(Vn)))+(Vn<-lamda*sigma+median(nonzeros(Vn)));
MS.onset  = 1+find(diff(overthreshold)== 1) - Tcorrection;
MS.offset = 1+find(diff(overthreshold)==-1) + Tcorrection;
tmpindex = unique([includeCheck(MSdomain,MS.onset) includeCheck(MSdomain,MS.offset)]);
MS.onset = MS.onset(tmpindex);
MS.offset = MS.offset(tmpindex);


% [4] Remove Inappropriate Trials
%  (2) Trials which occur just after the previous trials
    remindex = 1+find(MS.onset(2:end)-MS.offset(1:end-1)<=adjaThr);
    MS = removeMS(MS,remindex,0);  
%  (1) Trials with short duration
    remindex = find(MS.offset-MS.onset<=tooclose+2*Tcorrection | MS.offset-MS.onset>=longDur);
    MS = removeMS(MS,remindex,0);

% [5] Calculate Microsaccades Properties 
    MS.duration = zeros(length(MS.onset),1);
    MS.amp = zeros(length(MS.onset),1);
    MS.peakvel = zeros(length(MS.onset),1);
    MS.overshoot = zeros(length(MS.onset),1);
    MS.peaktime = zeros(length(MS.onset),1);
    MS.overtime = zeros(length(MS.onset),1);
    MS.starttime = zeros(length(MS.onset),1);
    MS.endtime = zeros(length(MS.onset),1);
    MS.endvel = zeros(length(MS.onset),1);
    MS.double = zeros(length(MS.onset),1);
    MS.damp = zeros(length(MS.onset),1);
    MS.tremor = zeros(length(MS.onset),1);
    MS.zerotime = zeros(length(MS.onset),1);
    MS.sigma = sigma;
    MS.vel = px_to_deg*Vn;
    for i=1:length(MS.onset)
        MS.duration(i) = (MS.offset(i)-MS.onset(i))/smpfreq*1000;
%         MS.amp(i) = px_to_MOA*(mean(Xn(MS.offset(i):min([length(Xn) MS.offset(i)+2*MAnum-1])))-mean(Xn(max([1 MS.onset(i)-2*MAnum+1]):MS.onset(i))));
        MS.amp(i) = 0;
        MS.peakvel(i) = px_to_deg*Vn(MS.onset(i)-1+argmax(abs(Vn(MS.onset(i):MS.offset(i)))));
        MS.peaktime(i) = MS.onset(i) -1 + argmax(abs(Vn(MS.onset(i):MS.offset(i))));
        tmppdx = MSfindpeak(MS.vel,MS.peaktime(i),2*(0.5-(MS.peakvel(i)>0)),1,3);
        if isnan(tmppdx)==1
            MS.double(i)=inf;
            continue
        end
        MS.double(i) = argmax(abs(MS.peakvel(i)-MS.vel(tmppdx)));
        MS.overtime(i) = tmppdx(MS.double(i));
        try
            MS.endtime(i) = MSfindpeak(MS.vel,MS.overtime(i),-2*(0.5-(MS.peakvel(i)>0)),1,1);
        catch myerr
            MS.double(i)=inf;
            continue
        end
        MS.starttime(i) = MSfindpeak(MS.vel,MS.peaktime(i),2*(0.5-(MS.peakvel(i)>0)),-1,1);
        MS.overshoot(i) = MS.vel(MS.overtime(i));
        MS.endvel(i)    = MS.vel(MS.endtime(i));
        MS.damp(i) = abs(MS.endvel(i)-MS.overshoot(i))/abs(MS.peakvel(i)-MS.overshoot(i));        % Damping Ratio (not Damping Factor)
%         MS.damp(i) = abs(MS.endvel(i)-MS.overshoot(i))/abs(MS.peakvel(i)-MS.endvel(i));
        MS.tremor(i) = (MS.damp(i)*MS.peakvel(i)+MS.overshoot(i))/(1+MS.damp(i));
        MS.zerotime(i) = MS.peaktime(i) + argmax(-abs(MS.vel(MS.peaktime(i)+1:MS.overtime(i)-1)-MS.tremor(i)));
        MS.amp(i) = px_to_MOA*(Xn(MS.endtime(i))-Xn(MS.onset(i)));
    end

% [4] Remove Large Saccades
    remindex = find(MS.double==inf);
    MS = removeMS(MS,remindex,1);
    remindex = find(abs(MS.amp)>maxamp | abs(MS.peakvel)>maxvel | abs(MS.overshoot)>maxvel);
    MS = removeMS(MS,remindex,1);
    
% [6] servo properties
    MS.zeta   = 1./sqrt(1+(pi./log(MS.damp)).^2);                   % Damping Factor
    MS.Tp = MS.overtime - MS.peaktime;                              % Duration
    MS.wn = pi./(MS.Tp.*sqrt(1-MS.zeta.^2));                        % Natural Frequency
    MS.Ka = MS.wn./(2*MS.zeta);                                     % Characteristic Coefficient of Servomechanism (?)
    if ~isempty(MS.onset)
        MS.preISI = [MS.onset(1) diff(MS.onset)']';                     % Inter-(micro)saccadic Interval (pre)
        MS.preISPC = Xn(MS.onset')'-Xn([1 MS.offset(1:end-1)']')';      % Inter-(micro)saccadic Position Change (pre)
        MS.postISI = [diff(MS.onset') length(Xn)-MS.onset(end)]';       % Inter-(micro)saccadic Interval (post): Cannot define for the last trial
        MS.postISPC = [Xn(MS.onset(2:end)) Xn(end)]' -Xn(MS.offset')';  % Inter-(micro)saccadic Position Change (post): Cannot define for the last trial
    end
    MS = rmfield(MS,{'vel','sigma'});
end

function m = argmax(x)
    [tmp,m]=max(x); %#ok<ASGLU>
    clear tmp
end
function index = includeCheck(domain,trial)
%   check if time domain include each trial
%   domain: <Nx2> [starttime endtime]
%   trial : <Kx1> [time]
    if size(domain,1)>1
        index = find(sum(((repmat(domain(:,1),1,length(trial))-repmat(trial',size(domain,1),1))>0) ~= ((repmat(domain(:,2),1,length(trial))-repmat(trial',size(domain,1),1))>0))==1);
    elseif size(domain,1)==1
        index = find(trial > domain(1) & trial < domain(2));
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
    end
    newMS = MS;
end
function pdx = MSfindpeak(vel,onset,sign,findvel,peaknum)
    % sign:  [-1] Search Local Minimum Point(s), [+1] Search Local Maximum Point(s)
    % delta: Sweep width
    % findvel:  [-1] Search Backward, [+1] Search Forward
    % peaknum: Number of the Peaks to Output
    delta = 100;
    [~,pdx]=findpeaks(sign*vel(onset+findvel:findvel:max([1 min([length(vel) onset+findvel*delta])])),'MINPEAKDISTANCE',5);
    loop = 0;
    nocorr = 0;
    while isempty(pdx)==1
        loop=loop+1;
        if onset+findvel*(loop+1)*delta > length(vel)
            pdx = NaN;
            nocorr = 1;
            break
        elseif onset+findvel*(loop+1)*delta < 1
            pdx = NaN;
            nocorr = 1;
            break
        end
        [~,pdx]=findpeaks(sign*vel(onset+findvel*(loop*delta+1):findvel:onset+findvel*(loop+1)*delta),'MINPEAKDISTANCE',5);
    end
    if nocorr==0
        pdx = pdx + loop*delta;
        pdx = onset + findvel*pdx(1:min([peaknum length(pdx)]));
    end
end