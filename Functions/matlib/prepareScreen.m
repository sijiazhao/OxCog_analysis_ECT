function scr=prepareScreen(ex)
% setup Screen (returns 'scr' struct) on display 
% creates Screen window
% preloads images and sounds
%
%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% returns:
%  scr.w        = window number
%  xhairCoords  = a crosshair as coords of two rects, at Screen centre
%
%  scr.soundData{i}  - sound file data          }
%  scr.soundFs{i}    - sound file format info   }
%  scr.imageData{i}  - image file data          }  of the given files
%  scr.imageMap{i}   - image file colour map    }  
%  scr.imageSize{i}  - dimensions of image (px) }
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% params:
%  ex.bgColour
%  ex.targetPos = [dx,dy] - returns scr.xhairCoords(:,:,2) and (:,:,3) as
%                          rects for crosshairs at these points relative to
%                          centre: [+dx,+dy], [-dx,-dy]
%  ex.imageFiles = {'f1'} - returns scr.imageData and scr.imageMap as data
%                          from image files image, and a openGL texture
%                          texture.
%  ex.imageAlpha = n      - colour index n is converted to the background
%                          colour index specified by ex.bgColourIndex
%  ex.soundFiles = {'f1'} - returns scr.soundData and scr.soundFs which 
%  ex.displayNumber = n   - use display number to open experiment Screen,
%                          default = display 0
%  ex.useCedrus           - try to connect to a 'cedrus' button box on COM6
% 
%  ex.useAlpha            - load up alpha (transparency) for each image and
%                           use that. also sets blendfunction to allow
%                           transparency
% 
%  ex.usePPA              - if using PsychPortAudio audio device, do not
%                           use audioplayer to load soundData into
%                           soundPlayer as this can cause errors with
%                           competing audio devices
% 
% Sanjay Manohar 2008

% suppress the PTB 'Welcome' screen
Screen('Preference', 'VisualDebugLevel',1);
% check if any screens already open? close them.
if prod(size(Screen('Screens'))), Screen('CloseAll'); end;
% whether to skip screen check
if isfield(ex,'skipScreenCheck')
    Screen('Preference', 'SkipSyncTests',ex.skipScreenCheck);
end
% use display 0 if none specified.
if ~isfield(ex,'displayNumber') ex.displayNumber=0; end;
% Open the window
[scr.w, scr.sszrect] = Screen('OpenWindow', ex.displayNumber, ex.bgColour);% [601 401 1400 1000]);

if isfield(ex, 'useAlpha') && ex.useAlpha % set for transparency
    Screen('BlendFunction', scr.w, 'GL_SRC_ALPHA', 'GL_ONE_MINUS_SRC_ALPHA');
end

Screen('Flip',scr.w);            % blank screen
scr.ssz=scr.sszrect(:,3:4);      % get screen dimensions in px
scr.centre=scr.ssz/2;            % centre of screen
if(isfield(ex,'xhairSize'))       % calculate crosshair coordinates if needed
    xh=repmat(ex.xhairSize,1,2) .* [-1 -1 1 1];
    scr.xhairCoords(:,:,1)=[ repmat(scr.centre,1,2)+xh([1 2 3 4]);...
                             repmat(scr.centre,1,2)+xh([2 1 4 3]) ]; %centre
    if(isfield(ex,'targetPos'))   % target pos = [x,y] distance from centre of screen
        scr.xhairCoords(:,:,2)=[ repmat(scr.centre+ex.targetPos,1,2)+xh([1 2 3 4]);...
                                 repmat(scr.centre+ex.targetPos,1,2)+xh([2 1 4 3]) ]; %right
        scr.xhairCoords(:,:,3)=[ repmat(scr.centre-ex.targetPos,1,2)+xh([1 2 3 4]);...
                                 repmat(scr.centre-ex.targetPos,1,2)+xh([2 1 4 3]) ]; %left
    end
end
if isfield(ex,'imageFiles')       % load images if requested
    for i=1:length(ex.imageFiles) % for each image
        fn=ex.imageFiles{i};      % get file name
        type = fn(end-2:end);     % get extension (last 3 chars)
        if fn(end-3)~='.' 
            if fn(end-4)=='.'
                type=fn(end-3:end);
            else
                type='gif';      % if no extension, then assume gif format.
            end
        end
        % load image from file - the pixel data and the colour map.
        [scr.imageData{i},scr.imageMap{i},alpha] = imread(ex.imageFiles{i},type);
        % now make the PsychToolbox texture.
        if(isfield(ex,'imageAlpha')) % check if 'imageAlpha' is specified
            discrimImage1(scr.imageData{i}==ex.imageAlpha) = ex.bgColourIndex;
            scr.imageTexture(i)=Screen('MakeTexture', scr.w, ...
               cat(3, scr.imageData{i}, scr.imageData{i}==ex.bgColourIndex));
        elseif isfield(ex, 'useAlpha') && ex.useAlpha
            scr.imageTexture(i)=Screen('MakeTexture', scr.w, ...
               cat(3, scr.imageData{i}, alpha));
        else
            scr.imageTexture(i)=Screen('MakeTexture', scr.w, ...
                scr.imageData{i});
        end
        % store the image size too.
        scr.imageSize{i} = [size(scr.imageData{i},2) size(scr.imageData{i},1)];
    end;
end;

if isfield(ex,'useCedrus')       % open the Cedrus buttonbox COM port?
  cedrus=serial('COM6');
  fopen(cedrus);
end

if isfield(ex,'soundFiles')           % preload sound files if required
    for i = 1:length(ex.soundFiles)     % for each sound file
      if isempty(ex.soundFiles{i}), continue; end % ignore blanks
      try                       
        if exist('audioread','file')  % try to read the file using audioread
            [scr.soundData{i}, scr.soundFs{i}]=audioread(ex.soundFiles{i});
            if ~isfield(ex, 'usePPA') || ~ex.usePPA % don't load this if using PPA
                scr.soundPlayer{i} = audioplayer(scr.soundData{i}, scr.soundFs{i});
            end
        else                          % older versions of matlab:
            [scr.soundData{i}, scr.soundFs{i}] = wavread(ex.soundFiles{i});
        end
      catch exc                       % explain if a file couldn't be read
        fprintf('Error reading %s\n',ex.soundFiles{i});
        rethrow(exc)
      end
    end
end
Screen(scr.w, 'TextSize', 32);        % default text font size 
Screen(scr.w, 'TextFont',  'Arial');  % and font.
HideCursor;                           % hide mouse pointer
FlushEvents '';                       % clear any keypresses in the buffer
