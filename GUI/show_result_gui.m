function varargout = show_result_gui(varargin)
% SHOW_RESULT_GUI MATLAB code for show_result_gui.fig
%      SHOW_RESULT_GUI, by itself, creates a new SHOW_RESULT_GUI or raises the existing
%      singleton*.
%
%      H = SHOW_RESULT_GUI returns the handle to a new SHOW_RESULT_GUI or the handle to
%      the existing singleton*.
%
%      SHOW_RESULT_GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SHOW_RESULT_GUI.M with the given input arguments.
%
%      SHOW_RESULT_GUI('Property','Value',...) creates a new SHOW_RESULT_GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before show_result_gui_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to show_result_gui_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help show_result_gui

% Last Modified by GUIDE v2.5 31-Jul-2013 04:51:12

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @show_result_gui_OpeningFcn, ...
                   'gui_OutputFcn',  @show_result_gui_OutputFcn, ...
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


% --- Executes just before show_result_gui is made visible.
function show_result_gui_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to show_result_gui (see VARARGIN)


handles.trace_file = '';
handles.node_file = '';
% Choose default command line output for show_result_gui
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes show_result_gui wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = show_result_gui_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on selection change in popup_map.
function popup_map_Callback(hObject, eventdata, handles)
% hObject    handle to popup_map (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popup_map contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popup_map

cla reset;
hold on; grid on; axis auto;
map = cellstr(get(hObject,'String'));
gui_mapCreate(map{get(hObject,'Value')});

if strcmp(handles.trace_file,'')|| strcmp(handles.node_file,'')
else
    gui_plot(handles.trace_file, handles.node_file,get(get(handles.panel_selectPlot, 'SelectedObject'),'Tag'));
end


% --- Executes during object creation, after setting all properties.
function popup_map_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popup_map (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in browse_trace_file.
function browse_trace_file_Callback(hObject, eventdata, handles)
% hObject    handle to browse_trace_file (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[filename pathname] = uigetfile({'*.txt'},'Select Trajectory File');
set(handles.trace_file_text,'String',[pathname filename]);
handles.trace_file = [pathname filename];

%Update handles structure
guidata(hObject, handles);

% --- Executes on button press in browse_node_file.
function browse_node_file_Callback(hObject, eventdata, handles)
% hObject    handle to browse_node_file (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[filename pathname] = uigetfile({'*.txt'},'Select Trajectory File');
set(handles.node_file_text,'String',[pathname filename]);
handles.node_file = [pathname filename];

%Update handles structure
guidata(hObject, handles);


% --- Executes on button press in button_calc_error.
function button_calc_error_Callback(hObject, eventdata, handles)
% hObject    handle to button_calc_error (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in button_plot.
function button_plot_Callback(hObject, eventdata, handles)
% hObject    handle to button_plot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
cla reset;
hold on; grid on; axis auto;
map = cellstr(get(handles.popup_map,'String'));
gui_mapCreate(map{get(handles.popup_map,'Value')});
if strcmp(handles.trace_file,'')|| strcmp(handles.node_file,'')
    msgbox('File(s) Missing!! Please browse the file','Error','error');
else
    gui_plot(handles.trace_file, handles.node_file,get(get(handles.panel_selectPlot, 'SelectedObject'),'Tag'));
end


% --- Executes during object creation, after setting all properties.
function axes1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to axes1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: place code in OpeningFcn to populate axes1
hold on; 
grid on; 
axis auto; 
xlabel('X (meter)');
ylabel('Y (meter)');


% --- Executes during object creation, after setting all properties.
function panel_selectPlot_CreateFcn(hObject, eventdata, handles)
% hObject    handle to panel_selectPlot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes when selected object is changed in panel_selectPlot.
function panel_selectPlot_SelectionChangeFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in panel_selectPlot 
% eventdata  structure with the following fields (see UIBUTTONGROUP)
%	EventName: string 'SelectionChanged' (read only)
%	OldValue: handle of the previously selected object or empty if none was selected
%	NewValue: handle of the currently selected object
% handles    structure with handles and user data (see GUIDATA)
guidata(hObject, handles);
