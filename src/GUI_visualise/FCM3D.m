function varargout = FCM3D(varargin)
% FCM3D MATLAB code for FCM3D.fig
%      FCM3D, by itself, creates a new FCM3D or raises the existing
%      singleton*.
%
%      H = FCM3D returns the handle to a new FCM3D or the handle to
%      the existing singleton*.
%
%      FCM3D('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in FCM3D.M with the given input arguments.
%
%      FCM3D('Property','Value',...) creates a new FCM3D or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before FCM3D_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to FCM3D_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help FCM3D

% Last Modified by GUIDE v2.5 22-Apr-2013 11:21:12

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @FCM3D_OpeningFcn, ...
                   'gui_OutputFcn',  @FCM3D_OutputFcn, ...
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


% --- Executes just before FCM3D is made visible.
function FCM3D_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to FCM3D (see VARARGIN)

% Choose default command line output for FCM3D
handles.output = hObject;
handles.dbt = [];
handles.mask = [];
handles.result = [];
handles.current_dbt = 1;
handles.current_result = 1;
handles.slices_dbt = 0;
handles.file_names = '';
handles.mean_dense = 0;
handles.mean_fatty = 0;
handles.v = zeros(1,3);

axes(handles.axes_load);
axis off

axes(handles.axes_result);
axis off

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes FCM3D wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = FCM3D_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in pushbutton_fcm3d.
function pushbutton_fcm3d_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_fcm3d (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Do the fuzzy clustering

%% Compile the c-code
mex src/3Dsegmentation/FCM/BCFCM3D.c -v

%% B: estimated biasfield
%% U: 3 classes of the partition matrix
[B,U] = BCFCM3D(handles.dbt.*handles.mask,handles.v,struct('maxit',5,'epsilon',1e-5,'sigma',1));
  
handles.result = classes2gray(U);

set(handles.edit_result,'String',['FCM 3D ',...
    ', Total: ',num2str(handles.current_result), '/',num2str(handles.slices_dbt)]);

axes(handles.axes_result)
imshow(handles.result(:,:,handles.current_result))
axis off

% Update handles structure
guidata(hObject, handles);

% --- Executes on button press in pushbutton_previous_result.
function pushbutton_previous_result_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_previous_result (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.current_result = handles.current_result - 1;

if handles.current_result < 1
    handles.current_result = 1;
end

set(handles.edit_result,'String',['FCM 3D ',...
    ', Total: ',num2str(handles.current_result), '/',num2str(handles.slices_dbt)]);

axes(handles.axes_result)
imshow(handles.result(:,:,handles.current_result),[])
axis off

% Update handles structure
guidata(hObject, handles);

% --- Executes on button press in pushbutton_next_result.
function pushbutton_next_result_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_next_result (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.current_result = handles.current_result + 1;

if handles.current_result > handles.slices_dbt
    handles.current_result = handles.slices_dbt;
end

set(handles.edit_result,'String',['FCM 3D ',...
    ', Total: ',num2str(handles.current_result), '/',num2str(handles.slices_dbt)]);

axes(handles.axes_result)
imshow(handles.result(:,:,handles.current_result),[])
axis off

% Update handles structure
guidata(hObject, handles);


function edit_result_Callback(hObject, eventdata, handles)
% hObject    handle to edit_result (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_result as text
%        str2double(get(hObject,'String')) returns contents of edit_result as a double


% --- Executes during object creation, after setting all properties.
function edit_result_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_result (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton_save.
function pushbutton_save_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_save (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%Path = uigetdir('/home/raquel/MT/thesis_tomosynthesis/resources/', 'Select folder to save ...');

seg_vol = handles.result;
filename = [Path,'fcm3d_',num2str(handles.v(2)),'_',num2str(handles.v(3)),'.mat'];

save(['fcm3d_',num2str(handles.v(2)),'_',num2str(handles.v(3)),'.mat'],'seg_vol');

% for i=1:handles.slices_dbt
%     tmp = handles.result(:,:,i);
%     imwrite(tmp,[Path,numPad(i,2),'.png']);
% end

% --- Executes on button press in pushbutton_load.
function pushbutton_load_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_load (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
Path = uigetfile(' ', 'Select a .mat file dbt images ...');
Path1 = uigetfile(' ', 'Select a .mat file for mask ...');

% read the images
%[D file_names]= read_dbt(Path,.2,0);
D = load(Path);
D = imresize(D.images,.5,'bicubic');

D_mask = load(Path1);
D_mask = imresize(D_mask.mask,.5,'bicubic');

% Convert to single
D = im2single(D);
handles.dbt = D;
handles.mask = D_mask;
handles.slices_dbt = size(D,3);
%handles.file_names = file_names;

axes(handles.axes_load)
imshow(handles.dbt(:,:,1),[])
axis off

%set(handles.edit_dbtimage,'String',[handles.file_names(handles.current_dbt+2).name,...
%    ', Total: ',num2str(handles.current_dbt), '/',num2str(handles.slices_dbt)]);

% Update handles structure
guidata(hObject, handles);

% --- Executes on button press in pushbutton_previous_load.
function pushbutton_previous_load_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_previous_load (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.current_dbt = handles.current_dbt - 1;

if handles.current_dbt < 1
    handles.current_dbt = 1;
end

set(handles.edit_dbtimage,'String',[handles.file_names(handles.current_dbt+2).name,...
    ', Total: ',num2str(handles.current_dbt), '/',num2str(handles.slices_dbt)]);

axes(handles.axes_load)
imshow(handles.dbt(:,:,handles.current_dbt),[])
axis off

% Update handles structure
guidata(hObject, handles);

% --- Executes on button press in pushbutton_next_load.
function pushbutton_next_load_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_next_load (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.current_dbt = handles.current_dbt + 1;

if handles.current_dbt > handles.slices_dbt
    handles.current_dbt = handles.slices_dbt;
end

set(handles.edit_dbtimage,'String',[handles.file_names(handles.current_dbt+2).name,...
    ', Total: ',num2str(handles.current_dbt), '/',num2str(handles.slices_dbt)]);

axes(handles.axes_load)
imshow(handles.dbt(:,:,handles.current_dbt),[])
axis off

% Update handles structure
guidata(hObject, handles);

function edit_dbtimage_Callback(hObject, eventdata, handles)
% hObject    handle to edit_dbtimage (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_dbtimage as text
%        str2double(get(hObject,'String')) returns contents of edit_dbtimage as a double


% --- Executes during object creation, after setting all properties.
function edit_dbtimage_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_dbtimage (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton_dense.
function pushbutton_dense_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_dense (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

figure('name','Select polygon dense tissue');
BW = roipoly(handles.dbt(:,:,handles.current_dbt));
tmp = handles.dbt(:,:,handles.current_dbt);
close;

handles.mean_dense = mean2(tmp(BW~=0));
set(handles.edit_mean_dense,'String',['Mean: ',num2str(handles.mean_dense)]);
handles.v(3) = handles.mean_dense;

% Update handles structure
guidata(hObject, handles);

function edit_mean_dense_Callback(hObject, eventdata, handles)
% hObject    handle to edit_mean_dense (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_mean_dense as text
%        str2double(get(hObject,'String')) returns contents of edit_mean_dense as a double


% --- Executes during object creation, after setting all properties.
function edit_mean_dense_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_mean_dense (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton_fatty.
function pushbutton_fatty_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_fatty (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
figure('name',' Select polygon fatty tissue');
BW = roipoly(handles.dbt(:,:,handles.current_dbt));
tmp = handles.dbt(:,:,handles.current_dbt);
close;

handles.mean_fatty = mean2(tmp(BW~=0));
set(handles.edit_mean_fatty,'String',['Mean: ',num2str(handles.mean_fatty)]);
handles.v(2) = handles.mean_fatty;
% Update handles structure
guidata(hObject, handles);


function edit_mean_fatty_Callback(hObject, eventdata, handles)
% hObject    handle to edit_mean_fatty (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_mean_fatty as text
%        str2double(get(hObject,'String')) returns contents of edit_mean_fatty as a double


% --- Executes during object creation, after setting all properties.
function edit_mean_fatty_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_mean_fatty (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
