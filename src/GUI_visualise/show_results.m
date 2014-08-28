function varargout = show_results(varargin)
%SHOW_RESULTS M-file for show_results.fig
%      SHOW_RESULTS, by itself, creates a new SHOW_RESULTS or raises the existing
%      singleton*.
%
%      H = SHOW_RESULTS returns the handle to a new SHOW_RESULTS or the handle to
%      the existing singleton*.
%
%      SHOW_RESULTS('Property','Value',...) creates a new SHOW_RESULTS using the
%      given property value pairs. Unrecognized properties are passed via
%      varargin to show_results_OpeningFcn.  This calling syntax produces a
%      warning when there is an existing singleton*.
%
%      SHOW_RESULTS('CALLBACK') and SHOW_RESULTS('CALLBACK',hObject,...) call the
%      local function named CALLBACK in SHOW_RESULTS.M with the given input
%      arguments.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help show_results

% Last Modified by GUIDE v2.5 21-Nov-2013 12:12:14

% Begin initialization code - DO NOT EDIT
gui_Singleton = 0;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @show_results_OpeningFcn, ...
                   'gui_OutputFcn',  @show_results_OutputFcn, ...
                   'gui_LayoutFcn',  [], ...
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


% --- Executes just before show_results is made visible.
function show_results_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   unrecognized PropertyName/PropertyValue pairs from the
%            command line (see VARARGIN)

% Choose default command line output for show_results
handles.output = hObject;
handles.seg = cell2mat(varargin(1));
handles.current_result = 1;
handles.slices_dbt = size(handles.seg,3);

handles.max = max(max(max(handles.seg)));
handles.min = min(min(min(handles.seg)));

% 
% axes(handles.axes_result);
% axis off

axes(handles.axes_result);
imshow(handles.seg(:,:,1),[handles.min handles.max]);

set(handles.edit_result,'String',['MRI slices',...
    ', Total: ',num2str(handles.current_result), '/',num2str(handles.slices_dbt)]);

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes show_results wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = show_results_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in pushbutton_previous_result.
function pushbutton_previous_result_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_previous_result (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.current_result = handles.current_result - 1;

if handles.current_result < 1
    handles.current_result = 1;
end

set(handles.edit_result,'String',['MRI slices',...
    ', Total: ',num2str(handles.current_result), '/',num2str(handles.slices_dbt)]);

axes(handles.axes_result)
imshow(handles.seg(:,:,handles.current_result),[handles.min handles.max])
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

set(handles.edit_result,'String',['MRI slices',...
    ', Total: ',num2str(handles.current_result), '/',num2str(handles.slices_dbt)]);

axes(handles.axes_result)
imshow(handles.seg(:,:,handles.current_result),[handles.min handles.max])
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


% --- Executes on button press in pushbutton_close.
function pushbutton_close_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_close (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%Path = uigetdir('/home/raquel/MT/thesis_tomosynthesis/resources/', 'Select folder to save ...');

close all;

% seg_vol = handles.seg;
% 
% save(['otsu3d_','.mat'],'seg_vol');

% for i=1:handles.slices_dbt
%     tmp = handles.result(:,:,i);
%     imwrite(tmp,[Path,numPad(i,2),'.png']);
% end
