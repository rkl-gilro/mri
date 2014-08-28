function varargout = choose_ax_sag_cor(varargin)
% CHOOSE_AX_SAG_COR MATLAB code for choose_ax_sag_cor.fig
%      CHOOSE_AX_SAG_COR, by itself, creates a new CHOOSE_AX_SAG_COR or raises the existing
%      singleton*.
%
%      H = CHOOSE_AX_SAG_COR returns the handle to a new CHOOSE_AX_SAG_COR or the handle to
%      the existing singleton*.
%
%      CHOOSE_AX_SAG_COR('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in CHOOSE_AX_SAG_COR.M with the given input arguments.
%
%      CHOOSE_AX_SAG_COR('Property','Value',...) creates a new CHOOSE_AX_SAG_COR or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before choose_ax_sag_cor_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to choose_ax_sag_cor_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help choose_ax_sag_cor

% Last Modified by GUIDE v2.5 20-May-2014 11:37:34

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @choose_ax_sag_cor_OpeningFcn, ...
                   'gui_OutputFcn',  @choose_ax_sag_cor_OutputFcn, ...
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


% --- Executes just before choose_ax_sag_cor is made visible.
function choose_ax_sag_cor_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to choose_ax_sag_cor (see VARARGIN)

% Choose default command line output for choose_ax_sag_cor
handles.output = hObject;

handles.series      = varargin{1}; % struct containing all the series (images & info)
handles.name_series = varargin{2}; % cell containing string with the series names

%% Initialize the list with the first series description for all
%% the 3 views
for i=1:length(handles.series)
    text{i} = handles.series{i}.info{1}.SeriesDescription;
end

set(handles.popupmenu_axial,'String',text);
set(handles.popupmenu_sagittal,'String',text);
set(handles.popupmenu_coronal,'String',text);

%% The selected series for each view (axial, sagittal and coronal)
handles.axial_series    = [];
handles.sagittal_series = [];
handles.coronal_series  = [];

%% Output struct containing the images and dicom info from the 3 orthogonal
%% views
handles.views = [];

%% Initialize axes
axes(handles.axes_image);
zoom on
axis off

axes(handles.axes_logo);
imshow(imread('pictures/icon/logo.png'),[]);
 
% Update handles structure
guidata(hObject, handles);

% UIWAIT makes choose_ax_sag_cor wait for user response (see UIRESUME)
uiwait(hObject);


% --- Outputs from this function are returned to the command line.
function varargout = choose_ax_sag_cor_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.views;


% --- Executes on selection change in popupmenu_coronal.
function popupmenu_coronal_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_coronal (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu_coronal contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_coronal

contents = cellstr(get(hObject,'String'));
handles.coronal_series = contents{get(handles.popupmenu_coronal,'Value')};

handles.views.coronal      = cellimages2mat(handles.series{get(handles.popupmenu_coronal,'Value')}.images);
handles.views.coronal_info = handles.series{get(handles.popupmenu_coronal,'Value')}.info;

num_slices = size(handles.views.coronal,3);
axes(handles.axes_image);
imshow(handles.views.coronal(:, :, round(num_slices / 2)), []);

% Update handles structure
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function popupmenu_coronal_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_coronal (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




% --- Executes on selection change in popupmenu_sagittal.
function popupmenu_sagittal_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_sagittal (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu_sagittal contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_sagittal

contents = cellstr(get(hObject,'String'));
handles.sagittal_series = contents{get(handles.popupmenu_sagittal,'Value')};

handles.views.sagittal      = cellimages2mat(handles.series{get(handles.popupmenu_sagittal,'Value')}.images);
handles.views.sagittal_info = handles.series{get(handles.popupmenu_sagittal,'Value')}.info;

num_slices = size(handles.views.sagittal,3);
axes(handles.axes_image);
imshow(handles.views.sagittal(:, :, round(num_slices / 2)), []);

% Update handles structure
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function popupmenu_sagittal_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_sagittal (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenu_axial.
function popupmenu_axial_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_axial (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu_axial contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_axial

contents = cellstr(get(hObject,'String'));
handles.axial_series = contents{get(handles.popupmenu_axial,'Value')};

handles.views.axial      = cellimages2mat(handles.series{get(handles.popupmenu_axial,'Value')}.images);
handles.views.axial_info = handles.series{get(handles.popupmenu_axial,'Value')}.info;

num_slices = size(handles.views.axial,3);
axes(handles.axes_image);
imshow(handles.views.axial(:, :, round(num_slices / 2)), []);

% Update handles structure
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function popupmenu_axial_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_axial (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton_done.
function pushbutton_done_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_done (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargout{1} = handles.output;

guidata(hObject,handles);

uiresume;

