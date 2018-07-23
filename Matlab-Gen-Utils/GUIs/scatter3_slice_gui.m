function varargout = scatter3_slice_gui(varargin)
% SCATTER3_SLICE_GUI create a GUI to display slices of 3D scattered data
%
%   SCATTER3_SLICE_GUI( X, Y, Z ) create the 3D scatter plot with only X,
%   Y, and Z data (no size or color).
%
%   SCATTER3_SLICE_GUI( X, Y, Z, S ) create the 3D scatter plot with size
%   determined by S.
%
%   SCATTER3_SLICE_GUI( X, Y, Z, S, C ) create the 3D scatter plot with
%   size determined by S and color by C. A colorbar is automatically added.
%   To provide color data but not size data, pass [] as S.
%
%   SCATTER3_SLICE_GUI( ___, 'scatter2' ) will use a 2D scatter plot
%   instead. Here, Z will only be used in figuring out what subset of the
%   data to plot for each slice. Essentially, this is like a top-down view
%   of the normal 3D scatter plot.
%
%   Parameters:
%       'xlabel', 'ylabel', 'zlabel', 'clabel', 'title' - each expects a
%       string, and each sets the X, Y, Z, or colorbar label of the scatter
%       plot or sets its title.


% SCATTER3_SLICE_GUI MATLAB code for scatter3_slice_gui.fig
%      SCATTER3_SLICE_GUI, by itself, creates a new SCATTER3_SLICE_GUI or raises the existing
%      singleton*.
%
%      H = SCATTER3_SLICE_GUI returns the handle to a new SCATTER3_SLICE_GUI or the handle to
%      the existing singleton*.
%
%      SCATTER3_SLICE_GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SCATTER3_SLICE_GUI.M with the given input arguments.
%
%      SCATTER3_SLICE_GUI('Property','Value',...) creates a new SCATTER3_SLICE_GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before scatter3_slice_gui_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to scatter3_slice_gui_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help scatter3_slice_gui

% Last Modified by GUIDE v2.5 15-Jan-2018 14:21:33

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @scatter3_slice_gui_OpeningFcn, ...
                   'gui_OutputFcn',  @scatter3_slice_gui_OutputFcn, ...
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


% --- Executes just before scatter3_slice_gui is made visible.
function scatter3_slice_gui_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to scatter3_slice_gui (see VARARGIN)

E = JLLErrors;

% Choose default command line output for scatter3_slice_gui
handles.output = hObject;

% A constant used across several functions
handles.internal.number_format = '%.4g';

% First look for look for flag strings in the input. "scatter2" indicates
% that a 2D scatter plot should be used instead of 3D. 
use_scatter2 = strcmpi(varargin, 'scatter2');
if any(use_scatter2)
    varargin(use_scatter2) = [];
    handles.internal.plot_type = 'scatter2';
else
    handles.internal.plot_type = 'scatter3';
end

% Now get the data, axis labels, and title
p = inputParser;
p.addRequired('X',@isnumeric);
p.addRequired('Y',@isnumeric);
p.addRequired('Z',@isnumeric);
p.addOptional('S',[],@isnumeric);
p.addOptional('C',[],@isnumeric);
p.addParameter('xlabel','',@ischar);
p.addParameter('ylabel','',@ischar);
p.addParameter('zlabel','',@ischar);
p.addParameter('clabel','',@ischar);
p.addParameter('title','',@ischar);

p.parse(varargin{:});
pout = p.Results;

% Ensure all inputs are transformed into vectors
handles.input.x_vals = pout.X(:);
handles.input.y_vals = pout.Y(:);
handles.input.z_vals = pout.Z(:);
if isempty(pout.S)
    % The default size is usually 36 for scatter plots
    handles.input.size_vals = 36*ones(size(handles.input.x_vals));
else
    handles.input.size_vals = pout.S(:);
end

% If no color is given, we also don't want to add a colorbar, since it will
% take up screen space for no reason
if isempty(pout.C)
    handles.input.color_vals = zeros(size(handles.input.x_vals));
    handles.options.add_colorbar = false;
else
    handles.input.color_vals = pout.C(:);
    handles.options.add_colorbar = true;
end

% Store the labels
handles.internal.xlabel = pout.xlabel;
handles.internal.ylabel = pout.ylabel;
handles.internal.zlabel = pout.zlabel;
handles.internal.clabel = pout.clabel;

% The title only needs set one
handles.text_title.String = pout.title;

% Use this to set the x and y limits to always stay fixed
handles.internal.xlim = calc_plot_limits(handles.input.x_vals,1,'pow10');
handles.internal.ylim = calc_plot_limits(handles.input.y_vals,1,'pow10');

% For now, we'll just default to choosing an initial z-bin size that will
% give 10 bins. Later I might add options that will let us choose the
% desired number of bins or the initial bin.
z_range = calc_plot_limits(handles.input.z_vals,1,'pow10');
handles.internal.zlim = z_range;
z_range(2) = z_range(2) / 10;
handles.internal.z_bin_range = z_range;

% By default, fix the z limit and color limit
handles.internal.clim = calc_plot_limits(handles.input.color_vals,1,'pow10');
handles.internal.fix_z = true;
handles.internal.fix_c = true;

handles = update_limit_controls(handles);
handles = update_plot(handles);

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes scatter3_slice_gui wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = scatter3_slice_gui_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function edit_zdata_min_Callback(hObject, eventdata, handles)
% hObject    handle to edit_zdata_min (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_zdata_min as text
%        str2double(get(hObject,'String')) returns contents of edit_zdata_min as a double
handles = text_update_callback(handles, hObject, 'z_bin_range', 1);
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function edit_zdata_min_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_zdata_min (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_zdata_max_Callback(hObject, eventdata, handles)
% hObject    handle to edit_zdata_max (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_zdata_max as text
%        str2double(get(hObject,'String')) returns contents of edit_zdata_max as a double
handles = text_update_callback(handles, hObject, 'z_bin_range', 2);
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function edit_zdata_max_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_zdata_max (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in button_zjump_down.
function button_zjump_down_Callback(hObject, eventdata, handles)
% hObject    handle to button_zjump_down (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Move the z bin range down so that the bin has the same width and the
% current bottom is the new top
limits = handles.internal.z_bin_range;
new_bottom = min(limits) - (limits(2) - limits(1));
handles.internal.z_bin_range = [new_bottom, limits(1)];
handles = update_limit_controls(handles);
handles = update_plot(handles);

guidata(hObject, handles);

% --- Executes on button press in button_zjump_up.
function button_zjump_up_Callback(hObject, eventdata, handles)
% hObject    handle to button_zjump_up (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Move the z bin range up so that the bin has the same width and the
% current top is the new bottom
limits = handles.internal.z_bin_range;
new_top = max(limits) + (limits(2) - limits(1));
handles.internal.z_bin_range = [limits(2), new_top];
handles = update_limit_controls(handles);
handles = update_plot(handles);

guidata(hObject, handles);

% --- Executes on button press in checkbox_fix_zlim.
function checkbox_fix_zlim_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_fix_zlim (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_fix_zlim
handles.internal.fix_z = get(hObject,'Value');
handles = update_limit_controls(handles);
handles = update_plot(handles);
guidata(hObject, handles);


function edit_zlim_min_Callback(hObject, eventdata, handles)
% hObject    handle to edit_zlim_min (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_zlim_min as text
%        str2double(get(hObject,'String')) returns contents of edit_zlim_min as a double
handles = text_update_callback(handles, hObject, 'zlim', 1);
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function edit_zlim_min_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_zlim_min (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_zlim_max_Callback(hObject, eventdata, handles)
% hObject    handle to edit_zlim_max (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_zlim_max as text
%        str2double(get(hObject,'String')) returns contents of edit_zlim_max as a double
handles = text_update_callback(handles, hObject, 'zlim', 2);
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function edit_zlim_max_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_zlim_max (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkbox_fix_clim.
function checkbox_fix_clim_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_fix_clim (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_fix_clim
handles.internal.fix_c = get(hObject,'Value');
handles = update_limit_controls(handles);
handles = update_plot(handles);
guidata(hObject, handles);


function edit_clim_min_Callback(hObject, eventdata, handles)
% hObject    handle to edit_clim_min (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_clim_min as text
%        str2double(get(hObject,'String')) returns contents of edit_clim_min as a double
handles = text_update_callback(handles, hObject, 'clim', 1);
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function edit_clim_min_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_clim_min (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function edit_clim_max_Callback(hObject, eventdata, handles)
% hObject    handle to edit_clim_max (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_clim_max as text
%        str2double(get(hObject,'String')) returns contents of edit_clim_max as a double
handles = text_update_callback(handles, hObject, 'clim', 2);
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function edit_clim_max_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_clim_max (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%%%%%%%%%%%%%%%%%%%%%%%%
% ADDITIONAL FUNCTIONS %
%%%%%%%%%%%%%%%%%%%%%%%%
function handles = update_limit_controls(handles)
% Updates the z range, z limits, and color limits text boxes to match
% what's recorded in the handles structure. Also updates the checkboxes for
% whether the color and z limits are fixed

% Ensure that the ranges are always [min, max]
handles.internal.z_bin_range = [min(handles.internal.z_bin_range), max(handles.internal.z_bin_range)];
handles.internal.zlim = [min(handles.internal.zlim), max(handles.internal.zlim)];
handles.internal.clim = [min(handles.internal.clim), max(handles.internal.clim)];

num_fmt = handles.internal.number_format;
handles.edit_zdata_min.String = sprintf(num_fmt, handles.internal.z_bin_range(1));
handles.edit_zdata_max.String = sprintf(num_fmt, handles.internal.z_bin_range(2));
handles.edit_zlim_min.String = sprintf(num_fmt, handles.internal.zlim(1));
handles.edit_zlim_max.String = sprintf(num_fmt, handles.internal.zlim(2));
handles.edit_clim_min.String = sprintf(num_fmt, handles.internal.clim(1));
handles.edit_clim_max.String = sprintf(num_fmt, handles.internal.clim(2));

handles.checkbox_fix_zlim.Value = handles.internal.fix_z;
handles.checkbox_fix_clim.Value = handles.internal.fix_c;


function handles = text_update_callback(handles, this_handle, limit_field, limit_index)
% Given the structure of handles, the specific handle to read from, which
% limit field to modify, and which index in the limit field to modify, this
% will handle parsing the input text and updating the plot. If an invalid
% number is given, it will reset the modified field to the old value. The
% value modified must be given by:
%   HANDLES.internal.(LIMIT_FIELD)(LIMIT_INDEX)

new_value = str2double(get(this_handle,'String'));
if isnan(new_value)
    % Bad input: need to reset the string to the old value
    old_string = sprintf(handles.internal.number_format, handles.internal.(limit_field)(limit_index));
    set(this_handle, 'String', old_string);
else
    handles.internal.(limit_field)(limit_index) = new_value;
    handles = update_limit_controls(handles);
    handles = update_plot(handles);
end


function handles = update_plot(handles)
E = JLLErrors;

X = handles.input.x_vals;
Y = handles.input.y_vals;
Z = handles.input.z_vals;
S = handles.input.size_vals;
C = handles.input.color_vals;

z_range = handles.internal.z_bin_range;
zz = Z >= min(z_range) & Z < max(z_range);

if strcmpi(handles.internal.plot_type, 'scatter2')
    wstate = warning('off','all'); % avoid warnings if C(zz) is an empty matrix
    scatter(handles.axes_scatter, X(zz), Y(zz), S(zz), C(zz));
    warning(wstate);
elseif strcmpi(handles.internal.plot_type, 'scatter3')
    wstate = warning('off','all'); % avoid warnings if C(zz) is an empty matrix
    scatter3(handles.axes_scatter, X(zz), Y(zz), Z(zz), S(zz), C(zz));
    warning(wstate)
else
    E.notimplemented('No plotting implemented for handles.internal.plot_type == "%s"', handles.internal.plot_type)
end

xlim(handles.axes_scatter, handles.internal.xlim);
ylim(handles.axes_scatter, handles.internal.ylim);

if handles.internal.fix_z
    zlim(handles.axes_scatter, handles.internal.zlim);
end
if handles.internal.fix_c
    caxis(handles.axes_scatter, handles.internal.clim)
end

xlabel(handles.axes_scatter, handles.internal.xlabel);
ylabel(handles.axes_scatter, handles.internal.ylabel);
zlabel(handles.axes_scatter, handles.internal.zlabel);

if handles.options.add_colorbar
    cb=colorbar;
    cb.Label.String = handles.internal.clabel;
end
