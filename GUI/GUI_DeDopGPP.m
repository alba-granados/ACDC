function varargout = GUI_DeDopGPP(varargin)
% GUI_DEDOPGPP MATLAB code for GUI_DeDopGPP.fig
%      GUI_DEDOPGPP, by itself, creates a new GUI_DEDOPGPP or raises the existing
%      singleton*.
%
%      H = GUI_DEDOPGPP returns the handle to a new GUI_DEDOPGPP or the handle to
%      the existing singleton*.
%
%      GUI_DEDOPGPP('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GUI_DEDOPGPP.M with the given input arguments.
%
%      GUI_DEDOPGPP('Property','Value',...) creates a new GUI_DEDOPGPP or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before GUI_DeDopGPP_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to GUI_DeDopGPP_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help GUI_DeDopGPP

% Last Modified by GUIDE v2.5 28-Nov-2016 15:53:42

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @GUI_DeDopGPP_OpeningFcn, ...
                   'gui_OutputFcn',  @GUI_DeDopGPP_OutputFcn, ...
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


% --- Executes just before GUI_DeDopGPP is made visible.
function GUI_DeDopGPP_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to GUI_DeDopGPP (see VARARGIN)

% Choose default command line output for GUI_DeDopGPP
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);


% UIWAIT makes GUI_DeDopGPP wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = GUI_DeDopGPP_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in start_processing.
function start_processing_Callback(hObject, eventdata, handles)
% hObject    handle to start_processing (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if(handles.cr2_button.Value)
    mission='CR2';
elseif(handles.s3_button.Value)
    mission='S3';
elseif(handles.s6_button.Value)
    mission='S6';
end

input_Path              = handles.input_folder_txt.String;
output_Path             = handles.output_folder_txt.String;
options.writting_flag   = [handles.l1bs.Value handles.l1b.Value handles.l1b_pLRM.Value handles.kml_check.Value];
options.plotting_flag   = [handles.plot_stacks_check.Value handles.plot_waveforms_check.Value handles.plotting_check.Value];
options.axes            = [];
options.axes            = handles.plotting_axes;
options.wd_axes         = [];
options.wd_axes         = handles.wd_axes;
options.kml_button      =[];
options.kml_button      = handles.open_kml_button;
options.GUI_flag        = 1;
options.pause_button    = [];
options.pause_button    = handles.pause_button;
options.draw_mask_check = [];
options.draw_mask_check = handles.draw_mask_check;
options.text_display    = [];
options.text_display    = handles.text_display;

hObject.ForegroundColor = [0.8 0.2 0.2];
hObject.BackgroundColor = [0.1 0.1 0.1];
hObject.String          = 'RUNNING';

% GPP_bulk_processing_paralelization('S3', './inputs/', './results/',1,0,[0 0 0])

GPP_bulk_processing_paralelization(mission, input_Path, output_Path,1,0,options);
if(handles.kml_check.Value)
    set(handles.open_kml_button, 'Enable', 'On');
end

hObject.String = 'DONE';
hObject.BackgroundColor = [0.94 0.94 0.94];

% --- Executes on button press in input_folder_button.
function input_folder_button_Callback(hObject, eventdata, handles)
% hObject    handle to input_folder_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
folder = uigetdir();
set(handles.input_folder_txt, 'String', folder);


% --- Executes on button press in output_folder_button.
function output_folder_button_Callback(hObject, eventdata, handles)
% hObject    handle to output_folder_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
folder = uigetdir();
set(handles.output_folder_txt, 'String', folder);



function input_folder_txt_Callback(hObject, eventdata, handles)
% hObject    handle to input_folder_txt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of input_folder_txt as text
%        str2double(get(hObject,'String')) returns contents of input_folder_txt as a double


% --- Executes during object creation, after setting all properties.
function input_folder_txt_CreateFcn(hObject, eventdata, handles)
% hObject    handle to input_folder_txt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function output_folder_txt_Callback(hObject, eventdata, handles)
% hObject    handle to output_folder_txt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of output_folder_txt as text
%        str2double(get(hObject,'String')) returns contents of output_folder_txt as a double


% --- Executes during object creation, after setting all properties.
function output_folder_txt_CreateFcn(hObject, eventdata, handles)
% hObject    handle to output_folder_txt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in plotting_check.
function plotting_check_Callback(hObject, eventdata, handles)
% hObject    handle to plotting_check (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if((hObject.Value==1))
    set(handles.plotting_axes, 'Visible', 'On');
    set(handles.wd_axes, 'Visible', 'On');
    set(handles.draw_mask_check, 'Enable', 'on');
    
else
    set(handles.plotting_axes, 'Visible', 'Off');
    set(handles.wd_axes, 'Visible', 'Off');
    set(handles.draw_mask_check, 'Enable', 'off');
end
% Hint: get(hObject,'Value') returns toggle state of plotting_check


% --- Executes on button press in plot_stacks_check.
function plot_stacks_check_Callback(hObject, eventdata, handles)
% hObject    handle to plot_stacks_check (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of plot_stacks_check


% --- Executes on button press in plot_waveforms_check.
function plot_waveforms_check_Callback(hObject, eventdata, handles)
% hObject    handle to plot_waveforms_check (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of plot_waveforms_check


% --- Executes on button press in kml_check.
function kml_check_Callback(hObject, eventdata, handles)
% hObject    handle to kml_check (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if((hObject.Value==1))
    set(handles.open_kml_button, 'Visible', 'On');
    
else
    set(handles.open_kml_button, 'Visible', 'Off');
    
end
% Hint: get(hObject,'Value') returns toggle state of kml_check


% --- Executes on button press in pause_button.
function pause_button_Callback(hObject, eventdata, handles)
% hObject    handle to pause_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
hObject.ForegroundColor = [0.8 0.2 0.2];
hObject.BackgroundColor = [0.3 0.3 0.3];
hObject.String          = 'PAUSING';

% --- Executes on button press in open_kml_button.
function open_kml_button_Callback(hObject, eventdata, handles)
kmlfile=dir([handles.output_folder_txt.String '/*kml']);
winopen(kmlfile.name);
% hObject    handle to open_kml_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in draw_mask_check.
function draw_mask_check_Callback(hObject, eventdata, handles)
% hObject    handle to draw_mask_check (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of draw_mask_check


% --------------------------------------------------------------------
function draw_mb_Callback(hObject, eventdata, handles)
% hObject    handle to draw_mb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function mb_draw_rect_Callback(hObject, eventdata, handles)
% hObject    handle to mb_draw_rect (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function Draw_cm_Callback(hObject, eventdata, handles)
% hObject    handle to Draw_cm (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function cm_draw_rect_Callback(hObject, eventdata, handles)
% hObject    handle to cm_draw_rect (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



% --- Executes during object creation, after setting all properties.
function plotting_axes_CreateFcn(hObject, eventdata, handles)
 
        
% hObject    handle to plotting_axes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: place code in OpeningFcn to populate plotting_axes


% --- Executes during object creation, after setting all properties.
function wd_axes_CreateFcn(hObject, eventdata, handles)
    figlabels('Elevation [m]','Latitude','','Window Elevation',10,hObject);
    set(hObject,'YAxisLocation','right');
    set(hObject,'YLim',[20 60]);
    
    


% hObject    handle to wd_axes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: place code in OpeningFcn to populate wd_axes


% --- Executes during object creation, after setting all properties.
function logo_dedop_CreateFcn(hObject, eventdata, handles)
try
    [dedop_logo]=imread('https://pbs.twimg.com/profile_images/723085760677724160/uf29s56c.jpg','jpg');
    imagesc(dedop_logo);
    axis off;    
catch
end
% hObject    handle to logo_dedop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: place code in OpeningFcn to populate logo_dedop


% --- Executes during object creation, after setting all properties.
function logo_isard_CreateFcn(hObject, eventdata, handles)
try
    [isard_logo]= imread('https://pbs.twimg.com/profile_images/378800000362712682/1dab150f95cb4b36f6ea9e579e2da783.jpeg','jpg');
    imagesc(isard_logo);
    axis off
catch
end
% hObject    handle to logo_isard (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: place code in OpeningFcn to populate logo_isard
