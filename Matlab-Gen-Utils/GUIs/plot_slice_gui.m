function varargout = plot_slice_gui(varargin)
% PLOT_SLICE_GUI creates a GUI interface to look at 2D slices in a 3D array
%   PLOT_SLICE_GUI( ARRAY ) will create a GUI to page through slices of the
%   3D ARRAY using default variables for the x and y coordinates using
%   pcolor to plot the values.
%
%   PLOT_SLICE_GUI( U, V ) will create a quiver plot to page through using
%   slices of U and V respectively as the x and y components.
%
%   PLOT_SLICE_GUI( ..., X, Y ) will use the values in the 2D arrays X and
%   Y to set the x and y coordinates for either of the syntaxes above.
%
%   PLOT_SLICE_GUI( ..., X, Y, LINE_X, LINE_Y ) will additionally plot a
%   line with x and y values given by LINE_X and LINE_Y on top of either
%   the pcolor or quiver plots.
%
%   PLOT_SLICE_GUI( ..., TITLE ) will include the string TITLE on the GUI
%   to help keep track of what is being plotted, which is especially
%   useful if you have multiple of these open. TITLE can be added to any of
%   the first syntaxes, but must be a string.
%
%   Josh Laughner <joshlaugh5@gmail.com> 


% Last Modified by GUIDE v2.5 12-Aug-2015 11:48:40

% Begin initialization code - DO NOT EDIT
gui_Singleton = 0;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @plot_slice_gui_OpeningFcn, ...
                   'gui_OutputFcn',  @plot_slice_gui_OutputFcn, ...
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


% --- Executes just before plot_slice_gui is made visible.
function plot_slice_gui_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to plot_slice_gui (see VARARGIN)

% Choose default command line output for plot_slice_gui
handles.output = hObject;

% First look for a string in the input - this will be the title. Find it,
% read it, remove it, and set the correct element's string.
iamstring = iscellcontents(varargin,'ischar');
if any(iamstring)
    if sum(iamstring) == 1
        gui_title = varargin{iamstring};
        varargin(iamstring) = [];
        set(handles.text_title,'String',gui_title);
    else
        E.badinput('Too many strings given - only 1 allowed')
    end
end


% Get the input matrix or matrices, checking that they are formatted
% correctly. Input can be 1 3D matrix of data or a 3D data matrix plus 2D
% matrices of x & y data
E=JLLErrors;
if numel(varargin) < 1
    E.badinput('Must input at least one matrix')
end

data = varargin{1};
handles.input.data = data;
sz = size(data);

% Will reset to quiver if 2 data arrays are given
handles.array_type = 'pcolor';

if mod(numel(varargin),2) == 0 
    data2 = varargin{2};
    sz2 = size(data2);
    handles.input.data2 = data2;
    
    if ndims(data) ~= ndims(data2) || ~all(sz == sz2)
        E.badinput('data1 and data2 must be the same size');
    end
    
    handles.array_type = 'quiver';
end

if ndims(data) ~= 3
    E.badinput('The input data (1st input) is expected to be 3D. This cannot handle 4D, and is a waste for 2D')
end

if numel(varargin)>2
    % Handle if user input x & y values
    if mod(numel(varargin),2)==1
        xvals = varargin{2};
        yvals = varargin{3};
        if numel(varargin)>3
            linex = varargin{4};
            liney = varargin{5};
        end
    elseif mod(numel(varargin),2)==0
        xvals = varargin{3};
        yvals = varargin{4};
        if numel(varargin)>4
            linex = varargin{5};
            liney = varargin{6};
        end
    end
    
    if size(xvals,1) ~= size(data,1) || size(xvals,2) ~= size(data,2)
        E.badinput('x values must have same first two dimension sizes as the data')
    elseif ~ismatrix(xvals)
        E.badinput('x values must be a 2D matrix')
    elseif size(yvals,1) ~= size(data,1) || size(yvals,2) ~= size(data,2)
        E.badinput('x values must have same first two dimension sizes as the data')
    elseif ~ismatrix(yvals)
        E.badinput('x values must be a 2D matrix')
    end
    
    
else
    % otherwise create default values
    [xvals, yvals] = meshgrid(1:sz(2), 1:sz(1));
end

% Add inputs to the handles structure
handles.input.xvals = xvals;
handles.input.yvals = yvals;
if exist('linex','var')
    handles.input.linex = linex;
    handles.input.liney = liney;
end
% Calculate the limits of the data - we'll use it to keep the colorbar with
% consistent limits
handles.input.datalims = [min(data(:)), max(data(:))];
set(handles.edit_cbmin,'String',num2str(handles.input.datalims(1)));
set(handles.edit_cbmax,'String',num2str(handles.input.datalims(2)));

% Set the initial value of the index (which slice in the third dimension
% we're plotting) plus the maximum value this index can take on.

handles.index.value = 1;
handles.index.max = sz(3);

% Default to faceted shading, allow faceted, flat, or interp.
handles.shading.curr_setting = 'faceted';
handles.shading.allowed = {'faceted', 'flat', 'interp'};
set(handles.menu_shading,'String', handles.shading.allowed);

% By default we will NOT fix the color range to the maximum and minimum of
% the data. This can be toggled
handles.input.fix_crange = 0;

% Disable controls not needed for the quiver plot
if strcmp(handles.array_type,'quiver')
    set(handles.menu_shading,'Enable','off');
    set(handles.checkbox_crange,'Enable','off');
end

% Make the initial plot
handles = update_plot(handles);

% Update handles structure
guidata(hObject, handles);
    
% UIWAIT makes plot_slice_gui wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = plot_slice_gui_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in button_next.
function button_next_Callback(hObject, eventdata, handles)
% hObject    handle to button_next (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Increment the index and update the plot - update_plot() will handle if
% the index is out-of-range
handles.index.value = handles.index.value + 1;
handles = update_plot(handles);
guidata(hObject,handles);


% --- Executes on button press in button_prev.
function button_prev_Callback(hObject, eventdata, handles)
% hObject    handle to button_prev (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Decrement the index and update the plot - update_plot() will handle if
% the index is out-of-range
handles.index.value = handles.index.value - 1;
handles = update_plot(handles);
guidata(hObject,handles);

% --- Executes on selection change in menu_shading.
function menu_shading_Callback(hObject, eventdata, handles)
% hObject    handle to menu_shading (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns menu_shading contents as cell array
%        contents{get(hObject,'Value')} returns selected item from menu_shading

if strcmp(handles.array_type,'pcolor')
    contents = cellstr(get(hObject,'String'));
    sel_shading = contents{get(hObject,'Value')};
    shading(sel_shading);
    handles.shading.curr_setting = sel_shading;
    
    guidata(hObject,handles);
end


% --- Executes during object creation, after setting all properties.
function menu_shading_CreateFcn(hObject, eventdata, handles)
% hObject    handle to menu_shading (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in checkbox_crange.
function checkbox_crange_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_crange (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_crange

if strcmp(handles.array_type,'pcolor')
    handles.input.fix_crange = get(hObject,'Value');
    handles = update_plot(handles);
    guidata(hObject,handles);
end
    
function edit_cbmin_Callback(hObject, eventdata, handles)
% hObject    handle to edit_cbmin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_cbmin as text
%        str2double(get(hObject,'String')) returns contents of edit_cbmin as a double

% Get the new user entered limit. If str2double returns a NaN, the user
% entered something that could not be parsed into a number, so reset the
% limit. Otherwise, save it in handles and update the plot with the new
% color range.
new_min = str2double(get(hObject,'String'));
if isnan(new_min)
    set(hObject,'String',num2str(handles.input.datalims(1)));
else
    handles.input.datalims(1) = new_min;
    update_plot(handles);
end

guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function edit_cbmin_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_cbmin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_cbmax_Callback(hObject, eventdata, handles)
% hObject    handle to edit_cbmax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_cbmax as text
%        str2double(get(hObject,'String')) returns contents of edit_cbmax as a double

% Get the new user entered limit. If str2double returns a NaN, the user
% entered something that could not be parsed into a number, so reset the
% limit. Otherwise, save it in handles and update the plot with the new
% color range.
new_max = str2double(get(hObject,'String'));
if isnan(new_max)
    set(hObject,'String',num2str(handles.input.datalims(2)));
else
    handles.input.datalims(2) = new_max;
    update_plot(handles);
end

guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function edit_cbmax_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_cbmax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% ADDITIONAL FUNCTIONS %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function handles = update_plot(handles)
% This function should be called in the callback function of any control
% that needs the entire plot to update (i.e. changing the slice, or
% possibly switching from pcolor to contour in the future). This isn't
% necessary for changing the pcolor shading. It will update the plot and
% the index display string.

% First check if the requested index to plot is allowed. If not, assume
% that no change is needed b/c we should only ever be changing the index by
% 1. Reset to an allowed value and return.
ind = handles.index.value;
indmax = handles.index.max;
if ind < 1
    handles.index.value = 1;
    return
elseif ind > indmax
    handles.index.value = indmax;
    return
end

% Otherwise we need to do 2 things: replot the slice and set the correct
% index in the display string.

data = handles.input.data(:,:,ind);
X = handles.input.xvals;
Y = handles.input.yvals;

% Set for pcolor or quiver
if strcmp(handles.array_type,'pcolor')
    pcolor(handles.axes_plot, X, Y, data);
    colorbar;
    
    % Reset the limits to the min and max of the data if we want to fix them to
    % be the same for all indicies.
    if handles.input.fix_crange
        caxis(handles.input.datalims);
    end
    shading(handles.shading.curr_setting); % pcolor will reset this
elseif strcmp(handles.array_type,'quiver')
    data2 = handles.input.data2(:,:,ind);
    quiver(handles.axes_plot,X,Y,data,data2);
end

if isfield(handles.input, 'linex')
    if numel(handles.input.linex) == 1
        markerspec = 'p';
    else
        markerspec = 'none';
    end
    line(handles.input.linex,handles.input.liney,'color','k','marker',markerspec,'markersize',12,'linewidth',1,'parent',handles.axes_plot);
end


handles.text_index.String = sprintf('%d/%d',ind,indmax);




