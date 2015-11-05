function varargout = MyCameraGUI(varargin)

%% GUI Initialization Code
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @MyCameraGUI_OpeningFcn, ...
    'gui_OutputFcn',  @MyCameraGUI_OutputFcn, ...
    'gui_LayoutFcn',  [] , ...
    'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


%% Executes just before MyCameraGUI is made visible.
function MyCameraGUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to MyCameraGUI (see VARARGIN)

% Choose default command line output for MyCameraGUI
handles.output = hObject;

%% Parameters Initialization
% Image Parameters
handles.sceneSize = [440,798];
handles.frameSize = round(handles.sceneSize/2);
handles.frameCenter = round(handles.frameSize/2);
handles.ROI = [154,248];
handles.returnedColorSpace = 'RGB';
handles.doFlip = 1;

% Edge Detection Filter Parameters
handles.EdgeFilter = 'prewitt';
handles.FilterThreshold = 0.03;             % Range [0.0 ~ 1.0]
handles.BWThreshold = 40;                   % Approx. Range [0 ~ 100]

% Angle range for Left and Right lanes
handles.minThetaL = 10;
handles.maxThetaL = 85;
handles.minThetaR = -85;
handles.maxThetaR = -5;

% Frame division factor for processing
handles.frameProcFactor = 2;                % Frame Processing Rate is improved with the increase in this number

% Sliding Window Filter Size for Hough Transform
handles.HoughSlidingWindowSize = 3;

% Hough Lane Detection Common Parameters
handles.noOfHoughPeaks = 10;                % Number of peaks returned by HoughPeaks function on passing Hough Transform
% Left Lane Specifics
handles.LeftNHood = [5 5];                  % Left Neighbourhood Size   % [Rho Theta], both of them must be odd
handles.LeftLineMinLength = 30;             % Left Line Min. Length     % Size is measured in terms of pixel length
handles.LeftGapFill = 4;                    % Gap Filling               % Size is measured in terms of pixel length
% Right Lane Specifics
handles.RightNHood = [21 7];                % Right Neighbourhood Size  % [Rho Theta], both of them must be odd
handles.RightLineMinLength = 30;            % Right Line Min. Length    % Size is measured in terms of pixel length
handles.RightGapFill = 4;                   % Gap Filling               % Size is measured in terms of pixel length

% Sample Processing Delay
SampleProcessingDelay = 0.02;               % In seconds

% This point helps in selecting the line closest to it
% For Left Lane
handles.x0L = 0.25*handles.frameSize(2);    % Note: Point location is measured from top-left
handles.y0L = handles.frameSize(1);
% For Right Lane
handles.x0R = 0.8*handles.frameSize(2);
handles.y0R = handles.frameSize(1);

% Length of Marker
handles.markerSize = 150;                   % in pixels

%% Masks

% BW Mask
handles.maskBW = zeros(handles.frameSize);  % The masking is applied to the edge detected image
handles.maskBW([140:round(handles.frameCenter(1)*7/4)],:) = 1;

% Window size of Rho and Theta to maintain to the previous location
handles.wRho = 3;
handles.wTheta = 3;
handles.AverageWeight = 3;

%% Initialization

%% Parameters initialization
handles.noOfRun = 0;
handles.myTime = 0;
handles.posBuffLength = 3;
handles.prevH = zeros(911,handles.maxThetaL-handles.minThetaL+handles.maxThetaR-handles.minThetaR+2,handles.HoughSlidingWindowSize);
handles.posBuffL = zeros(4,handles.posBuffLength);
handles.posBuffR = zeros(4,handles.posBuffLength);
handles.prevLinesL = struct([]);
handles.prevLinesR = struct([]);
handles.cam = [];
handles.slopeL = 0;
handles.xInterceptL = 0;
handles.slopeR = 0;
handles.xInterceptR = 0;
handles.posBuffLength = 3;

%% Code Initialization

axes(handles.cameraAxes);
white = uint8(255*ones(handles.sceneSize));
imshow(white);

tempObj = imaqhwinfo('winvideo');
tempVar = table2cell(struct2table(tempObj.DeviceInfo));
set(handles.camListPMenu,'string',tempVar(:,3));
handles.buf1 = get(handles.cameraAxes,'children');

handles.t = timer('TimerFcn',{@timelapse_timer,handles},'StartDelay',1,'Period',SampleProcessingDelay,'ExecutionMode','fixedSpacing');

%% GUI Specific Section
guidata(hObject, handles);
% UIWAIT makes MyCameraGUI wait for user response (see UIRESUME)
uiwait(handles.MyCameraGUI);

function timelapse_timer(hObject, eventdata, handles)

handles = guidata(gcf);

if(~isempty(gco))
    
    snap = getsnapshot(handles.cam);
    if handles.doFlip
        snap = flip(flip(snap,1),2);
    end;
    
    frame = rgb2gray(snap(1:handles.frameProcFactor:end,1:handles.frameProcFactor:end,:));
    
    BW = edge(frame,handles.EdgeFilter,handles.FilterThreshold);
    BW = handles.maskBW.*BW;
    
    [H,T,R] = hough(BW,'RhoResolution',1,'Theta',[handles.minThetaR:handles.maxThetaR,handles.minThetaL:handles.maxThetaL]);
    
    H(H < handles.BWThreshold) = 0;
    
    % --- Uncomment for Hough Transform Smoothing
    %         handles.prevH(:,:,2:end) = handles.prevH(:,:,1:end-1);
    %         handles.prevH(:,:,1) = H;
    %         H = (sum(handles.prevH,3) + H) / (handles.HoughSlidingWindowSize + 1);
    
    %% Detect left lane
    
    T1 = zeros([size(T,1),180]);
    T1(:,90+handles.minThetaL:90+handles.maxThetaL) = T(:,handles.maxThetaR-handles.minThetaR+2:end);
    H1 = zeros([size(H,1),180]);
    H1(:,90+handles.minThetaL:90+handles.maxThetaL) = H(:,handles.maxThetaR-handles.minThetaR+2:end);
    
    peaks = houghpeaks(round(H1), handles.noOfHoughPeaks, 'Threshold', handles.BWThreshold,'NHoodSize',handles.LeftNHood);
    
    if size(peaks) ~= 0
        lines = houghlines(BW, T1, R, peaks,'MinLength',handles.LeftLineMinLength,'FillGap',handles.LeftGapFill);
        
        [handles.posBuffL,handles.slopeL,handles.xInterceptL,handles.prevLinesL] = ...
            chooseXY(lines,handles.posBuffL,handles.frameCenter,handles.x0L,handles.y0L,handles.slopeL,handles.xInterceptL,H1,'left',...
            handles.prevLinesR,handles.wRho,handles.wTheta,handles.AverageWeight);
        if numel(handles.prevLinesL) ~= 0
            snap = drawLane(handles.slopeL, handles.xInterceptL, snap,'left',handles.markerSize);
        end;
    end;
    
    %% Detect right lane
    
    H2 = zeros([size(H,1),180]); H2(:,90+handles.minThetaR:90+handles.maxThetaR)=H(:,1:handles.maxThetaR-handles.minThetaR+1);
    
    T2 = zeros([size(T,1),180]);
    T2(:,90+handles.minThetaR:90+handles.maxThetaR)=T(:,1:handles.maxThetaR-handles.minThetaR+1);
    
    peaks = houghpeaks(round(H2), handles.noOfHoughPeaks, 'Threshold', handles.BWThreshold,'NHoodSize',handles.RightNHood);
    
    if numel(peaks) ~= 0
        lines = houghlines(BW, T2, R, peaks,'MinLength',handles.RightLineMinLength,'FillGap',handles.RightGapFill);
        
        [handles.posBuffR,handles.slopeR,handles.xInterceptR,handles.prevLinesR] = ...
            chooseXY(lines,handles.posBuffR,handles.frameCenter,handles.x0R,handles.y0R,handles.slopeR,handles.xInterceptR,H2,'right',...
            handles.prevLinesR,handles.wRho,handles.wTheta,handles.AverageWeight);
        if numel(handles.prevLinesR) ~= 0
            snap = drawLane(handles.slopeR, handles.xInterceptR, snap,'right',handles.markerSize);
        end;
    end;
    
    %%
    handles.buf1 = get(handles.cameraAxes,'children');
    set(handles.buf1,'cData',snap);
    set(handles.cameraAxes,'children',handles.buf1)
    
    guidata(gcf, handles);
end;


% --- Outputs from this function are returned to the command line.
function varargout = MyCameraGUI_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
% handles = guidata(gcf);
handles.output = hObject;
varargout{1} = handles.output;


% --- Executes when user attempts to close MyCameraGUI.
function MyCameraGUI_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to MyCameraGUI (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% handles = guidata(gcf);
delete('handles.cam');
stop(handles.t);
% Hint: delete(hObject) closes the figure
delete(hObject);
delete(imaqfind);
% guidata(gcf, handles);


% --- Executes on button press in startStopCamera.
function startStopCamera_Callback(hObject, eventdata, handles)
% hObject    handle to startStopCamera (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles = guidata(gcf);
if strcmp(get(handles.startStopCamera,'String'),'Start Camera') && isempty(get(handles.camListPMenu,'string')) == 0
    % Camera is off. Change button string and start camera.
    set(handles.startStopCamera,'String','Stop Camera');
    
    handles.cam = videoinput('winvideo', get(handles.camListPMenu,'value'));
    handles.cam.ROIPosition = [handles.ROI(2) handles.ROI(1) handles.sceneSize(2) handles.sceneSize(1)];
    handles.cam.returnedColorSpace = handles.returnedColorSpace;
    triggerconfig(handles.cam, 'manual');
    start(handles.cam);
    guidata(gcf, handles);
    start(handles.t);
else
    % Camera is on. Stop camera and change button string.
    set(handles.startStopCamera,'String','Start Camera')
    stop(handles.t);
    delete('handles.cam');
    guidata(gcf, handles);
    stop(handles.t);
end

% --- Executes on selection change in camListPMenu.
function camListPMenu_Callback(hObject, eventdata, handles)
% hObject    handle to camListPMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns camListPMenu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from camListPMenu

% --- Executes during object creation, after setting all properties.
function camListPMenu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to camListPMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
