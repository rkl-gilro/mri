function varargout = view3d(varargin)
% VIEW3D MATLAB code for view3d.fig
%      VIEW3D, by itself, creates a new VIEW3D or raises the existing
%      singleton*.
%
%      H = VIEW3D returns the handle to a new VIEW3D or the handle to
%      the existing singleton*.
%
%      VIEW3D('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in VIEW3D.M with the given input arguments.
%
%      VIEW3D('Property','Value',...) creates a new VIEW3D or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before view3d_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to view3d_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help view3d

% Last Modified by GUIDE v2.5 06-Mar-2014 12:51:09

% Begin initialization code - DO NOT EDIT
gui_Singleton = 0;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @view3d_OpeningFcn, ...
                   'gui_OutputFcn',  @view3d_OutputFcn, ...
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


% --- Executes just before view3d is made visible.
function view3d_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to view3d (see VARARGIN)

%handles.orient = imread('resources/orient2.png');

% Choose default command line output for view3d
handles.output = hObject;

handles.views = varargin{1};
handles.triangulation = [];

if length(varargin) > 1
    handles.triangulation = varargin{2};
end

handles.plot_slices = 1;
handles.plot_vol    = 1;
%% Axis orientation 3D
handles.az = -37.5;
handles.el = 30;

handles.vol_axial = handles.views.axial;
handles.r_axial = size(handles.vol_axial,1);
handles.c_axial = size(handles.vol_axial,2);
handles.s_axial = size(handles.vol_axial,3); %% number of slices
handles.slider_axial_val = 1;%% slider value


handles.vol_sagittal = handles.views.sagittal;
handles.r_sagittal = size(handles.vol_sagittal,1);
handles.c_sagittal = size(handles.vol_sagittal,2);
handles.s_sagittal = size(handles.vol_sagittal,3); %% number of slices
handles.slider_sagittal_val = 1;%% slider value

handles.vol_coronal = handles.views.coronal;
handles.r_coronal = size(handles.vol_coronal,1);
handles.c_coronal = size(handles.vol_coronal,2);
handles.s_coronal = size(handles.vol_coronal,3); %% number of slices
handles.slider_coronal_val = 1; %% slider value

handles.view3d = 0; % flag to determine if we should show the 3dview each time the slider changes

axes(handles.axes_axial);
imshow(handles.vol_axial(:,:,1),[]);

set(handles.edit_AXIAL,'String',['AXIAL:     #Slices ',...
   num2str(handles.slider_axial_val), '/',num2str(handles.s_axial)]);

set(handles.slider_axial,'min',1);
set(handles.slider_axial,'max',handles.s_axial);
set(handles.slider_axial, 'SliderStep', [1/handles.s_axial , 10/handles.s_axial ]);



axes(handles.axes_sagittal);
imshow(handles.vol_sagittal(:,:,1),[]);

set(handles.edit_SAG,'String',['SAGITTAL:    #Slices ',...
   num2str(handles.slider_sagittal_val), '/',num2str(handles.s_sagittal)]);

set(handles.slider_sagittal,'min',1);
set(handles.slider_sagittal,'max',handles.s_sagittal);
set(handles.slider_sagittal, 'SliderStep', [1/handles.s_sagittal, 10/handles.s_sagittal ]);



axes(handles.axes_coronal);
imshow(handles.vol_coronal(:,:,1),[]);

set(handles.edit_COR,'String',['CORONAL:     #Slices ',...
   num2str(handles.slider_coronal_val), '/',num2str(handles.s_coronal)]);

set(handles.slider_coronal,'min',1);
set(handles.slider_coronal,'max',handles.s_coronal);
set(handles.slider_coronal, 'SliderStep', [1/handles.s_coronal , 10/handles.s_coronal ]);



axes(handles.axes_3dview);
zoom on
axis off

axes(handles.axes_orient);
imshow(imread('pictures/icon/logo.png'),[]);
 
% Update handles structure
guidata(hObject, handles);

% UIWAIT makes view3d wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = view3d_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;




% --- Executes on button press in pushbutton_3dviews.
function pushbutton_3dviews_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_3dviews (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Flag 
handles.view3d = 1;

axes(handles.axes_3dview);
axis on
cla

% set(cla,'CameraViewAngle',handles.camera);
colors = [0 1 1];
%% Plot volume if it selected and it exists
if handles.plot_vol
    
    if ~isempty(handles.triangulation)
        
        for i = 1:length(handles.triangulation)
               
            colors = [i/length(handles.triangulation) 0 mod(i+1,2)]; 
            trisurf(handles.triangulation(i).faces, handles.triangulation(i).vertices(:,1),handles.triangulation(i).vertices(:,2),handles.triangulation(i).vertices(:,3),'facecolor',colors,'edgecolor','none','facelighting','flat');hold on
            alpha(.1*i)
        
        end
        
    end
end

if handles.plot_slices
    %% Axial %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ax = handles.slider_axial_val;
    
    aa(:,:,1) = convert2u8(handles.vol_axial(:,:,ax));
    aa(:,:,2) = convert2u8(handles.vol_axial(:,:,ax));
    aa(:,:,3) = convert2u8(handles.vol_axial(:,:,ax));
    
    [X, Y, Z, triTexture] = compute_triTexture(handles.views.axial_info{ax},aa,[0 handles.r_axial-1], [0 handles.c_axial-1]);
    surf(X,Y,Z,triTexture,'FaceColor','texturemap','EdgeColor','none');hold on
    axis equal
    
    
    %% Sagittal %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    sag = handles.slider_sagittal_val;
    
    ss(:,:,2) = convert2u8(handles.vol_sagittal(:,:,sag));
    ss(:,:,1) = convert2u8(handles.vol_sagittal(:,:,sag));
    ss(:,:,3) = convert2u8(handles.vol_sagittal(:,:,sag));
    
    [X, Y, Z, triTexture] = compute_triTexture(handles.views.sagittal_info{sag},ss,[0 handles.r_sagittal-1], [0 handles.c_sagittal-1]);
    surf(X,Y,Z,triTexture,'FaceColor','texturemap','EdgeColor','none');hold on
    axis equal
    
%     M = zeros(4,3);
%     
%     M(1:3,3) = handles.views.axial_info{ax}.ImagePositionPatient;
%     M(4,3) = 1;
%     M(1:3,1) = handles.views.axial_info{ax}.ImageOrientationPatient(1:3).*handles.views.axial_info{ax}.PixelSpacing(1);
%     M(1:3,2) = handles.views.axial_info{ax}.ImageOrientationPatient(4:6).*handles.views.axial_info{ax}.PixelSpacing(2);
%     
%     p1 = handles.views.axial_info{ax} * [300-1 200-1 1]';
%     p2 = handles.views.axial_info{ax} * [100-1 200-1 1]';
%     
%     plot3(p1(1), p1(2), p1(3),'m*');hold on
%     plot3(p2(1), p2(2), p2(3),'g*');hold on
    
    %% Coronal %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    cor = handles.slider_coronal_val;
    
    cc(:,:,1) = convert2u8(handles.vol_coronal(:,:,cor));
    cc(:,:,2) = convert2u8(handles.vol_coronal(:,:,cor));
    cc(:,:,3) = convert2u8(handles.vol_coronal(:,:,cor));
    
    [X, Y, Z, triTexture] = compute_triTexture(handles.views.coronal_info{cor},cc,[0 handles.r_coronal-1], [0 handles.c_coronal-1]);
    
    surf(X,Y,Z,triTexture,'FaceColor','texturemap','EdgeColor','none');hold on

end

axis equal
axis off


view([handles.az handles.el])



% handles.camera = get(cla,'CameraViewAngle');

% Update handles structure
guidata(hObject, handles);

% --- Executes on slider movement.
function slider_coronal_Callback(hObject, eventdata, handles)
% hObject    handle to slider_coronal (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
handles.slider_coronal_val = int32(get(handles.slider_coronal,'Value'));

set(handles.edit_COR,'String',['CORONAL:     #Slices ',...
   num2str(handles.slider_coronal_val), '/',num2str(handles.s_coronal)]);

axes(handles.axes_coronal);
imshow(handles.vol_coronal(:,:,handles.slider_coronal_val),[]);

if handles.view3d == 1
    [handles.az,handles.el] = view(handles.axes_3dview);
    pushbutton_3dviews_Callback(hObject, eventdata, handles);
end

% Update handles structure
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function slider_coronal_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider_coronal (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function slider_sagittal_Callback(hObject, eventdata, handles)
% hObject    handle to slider_sagittal (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
handles.slider_sagittal_val = int32(get(handles.slider_sagittal,'Value'));

set(handles.edit_SAG,'String',['SAGITTAL:    #Slices ',...
   num2str(handles.slider_sagittal_val), '/',num2str(handles.s_sagittal)]);

axes(handles.axes_sagittal);
imshow(handles.vol_sagittal(:,:,handles.slider_sagittal_val),[]);

if handles.view3d == 1
    [handles.az,handles.el] = view(handles.axes_3dview);
    pushbutton_3dviews_Callback(hObject, eventdata, handles);
end

% Update handles structure
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function slider_sagittal_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider_sagittal (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function slider_axial_Callback(hObject, eventdata, handles)
% hObject    handle to slider_axial (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
handles.slider_axial_val = int32(get(handles.slider_axial,'Value'));

set(handles.edit_AXIAL,'String',['AXIAL:     #Slices ',...
   num2str(handles.slider_axial_val), '/',num2str(handles.s_axial)]);

axes(handles.axes_axial);
imshow(handles.vol_axial(:,:,handles.slider_axial_val),[]);

if handles.view3d == 1
    [handles.az,handles.el] = view(handles.axes_3dview);
    pushbutton_3dviews_Callback(hObject, eventdata, handles);
end

% Update handles structure
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function slider_axial_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider_axial (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end



function edit_axial_Callback(hObject, eventdata, handles)
% hObject    handle to edit_axial (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_axial as text
%        str2double(get(hObject,'String')) returns contents of edit_axial as a double


% --- Executes during object creation, after setting all properties.
function edit_axial_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_axial (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_sagittal_Callback(hObject, eventdata, handles)
% hObject    handle to edit_sagittal (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_sagittal as text
%        str2double(get(hObject,'String')) returns contents of edit_sagittal as a double


% --- Executes during object creation, after setting all properties.
function edit_sagittal_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_sagittal (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_coronal_Callback(hObject, eventdata, handles)
% hObject    handle to edit_coronal (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_coronal as text
%        str2double(get(hObject,'String')) returns contents of edit_coronal as a double


% --- Executes during object creation, after setting all properties.
function edit_coronal_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_coronal (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_AXIAL_Callback(hObject, eventdata, handles)
% hObject    handle to edit_AXIAL (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_AXIAL as text
%        str2double(get(hObject,'String')) returns contents of edit_AXIAL as a double


% --- Executes during object creation, after setting all properties.
function edit_AXIAL_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_AXIAL (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_SAG_Callback(hObject, eventdata, handles)
% hObject    handle to edit_SAG (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_SAG as text
%        str2double(get(hObject,'String')) returns contents of edit_SAG as a double


% --- Executes during object creation, after setting all properties.
function edit_SAG_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_SAG (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_COR_Callback(hObject, eventdata, handles)
% hObject    handle to edit_COR (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_COR as text
%        str2double(get(hObject,'String')) returns contents of edit_COR as a double


% --- Executes during object creation, after setting all properties.
function edit_COR_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_COR (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton_close.
function pushbutton_close_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_close (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
close; 


% --- Executes on button press in radiobutton_vol.
function radiobutton_vol_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton_vol (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton_vol
if (get(hObject,'Value') == get(hObject,'Max'))
	% Radio button is selected-take appropriate action
    handles.plot_vol = 1;
    [handles.az,handles.el] = view(handles.axes_3dview);
    pushbutton_3dviews_Callback(hObject, eventdata, handles);
else
	% Radio button is not selected-take appropriate action
    handles.plot_vol = 0;
end
% Update handles structure
guidata(hObject, handles);

% --- Executes on button press in radiobutton_slices.
function radiobutton_slices_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton_slices (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton_slices

if (get(hObject,'Value') == get(hObject,'Max'))
	% Radio button is selected-take appropriate action
    handles.plot_slices = 1;
    [handles.az,handles.el] = view(handles.axes_3dview);
    pushbutton_3dviews_Callback(hObject, eventdata, handles);
else
	% Radio button is not selected-take appropriate action
    handles.plot_slices = 0;
end
% Update handles structure
guidata(hObject, handles);