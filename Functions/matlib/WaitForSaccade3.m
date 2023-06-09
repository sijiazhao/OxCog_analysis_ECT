function results=WaitForSaccade3(el, maxtime, minsize)
% results = WaitForSaccade3( el, maxtime, minsize )
% This uses the eyelink parser to determine when the next saccade occurs
% returns array [end-saccade-time, end-saccade-x, end-saccade-y]
% This version dequeues data from the link.
%
% waits until saccade end, or time elapsed, or if saccade is late then
% until saccade has completed. If saccade has already started, then
% behaviour depends on the flag el.disallowEarlySaccades.

startTime = GetSecs;
time=0;
started=0;ended=0;
% drain queue
oldevents=[];
drained=0; while(~drained)
    [samples, events,drained]=Eyelink('GetQueuedData');
    oldevents=[oldevents;events];
end
s=Eyelink('NewestFloatSample');
p0=[s.time s.gx(el.eye) s.gy(el.eye)];
timeout=0;
debug=isfield(el, 'debugSaccades');
disallowEarlySaccades = isfield(el, 'disallowEarlySaccades');


WAIT_FOR_ENDING = false; % wait for the end of the saccade before returning?


while ~ended % keep waiting until end signal reached

    if(GetSecs>startTime+maxtime)  % check if we have timed out
      ended=1;
      results=0;                   % return a blank, if so
    end

    drained=0;
    while(~drained)                % keep reading until the queue is empty
        [samples, events,drained]=Eyelink('GetQueuedData');
        if(~isempty(oldevents) && size(oldevents,2)==size(events,2))
            events=[oldevents;events];
            oldevents=[];
        end
        if(isempty(events))  continue;  end;
        endsaccs = find(events(2,:) == el.ENDSACC);  % find ENDSACC events
        if(length(endsaccs)>0)  % for each one,
            for(i=1:length(endsaccs))
              results = [events([1 14 15],endsaccs(i))]';  % add the saccade to the result
              dist    = norm(events([14,15],endsaccs(i)) - events([9,10], endsaccs(i)));
              if(dist>minsize)  % check if it was big enough, and exit immediately with results
                return; 
              end;
            end
        end;
    end;

    if WAIT_FOR_ENDING   % old version of WaitForSaccade2 - reads values manually
      % this allows checking for the beginnings and ends of saccades
      st=eyelink('getnextdatatype');
      if st==el.STARTSACC
          started=1;
      end;
      
      if st==el.ENDSACC
          s=eyelink('getfloatdata',st);
          p1 = [s.sttime s.gstx s.gsty];
          p2 = [s.entime s.genx s.geny];
          if (~started && s.sttime<p0(1)) & disallowEarlySaccades continue;end;
          if norm(p1(2:3)-p2(2:3))<minsize & norm(p2(2:3)-p0(2:3))<minsize
              started=0; ended=0;
          else
              ended=1;
              results = p2;
          end;
      end;
      if eyelink('NewFloatSampleAvailable')
          s2=eyelink('NewestFloatSample');
          time = s2.time - p0(1);
          if time>maxtime & ~started % if saccade already begun, allow it to complete
              ended=1;
              results = 0;
          end;
          if debug
              p=[s2.gx(el.eye) s2.gy(el.eye)]
              screen('fillrect', el.window, 0,[p(2) p(3) p(2)+10 p(3)+10]);
              screen('flip', el.window);
          end;
      end;
    end

    % check keyboard - and exit if a key is pressed.
    [z z keys] = KbCheck;
    if any(keys) ended=1; results=0; end;
end;


