function varargout = select_BL_heights(varargin)
% SELECT_BL_HEIGHTS_GUI Allows use to identify boundary layer heights from
% vertical profiles
%
%   [ struct_out ] = SELECT_BL_HEIGHTS_GUI(Merge, UTC_range, altfield_name,
%   field1_name, [opt] field2_name, field3_name)
%
%   Inputs required are a Merge data structure, a UTC time range as a 1x2
%   vector (seconds after midnight UTC), the name of the field containing
%   the altitude field and optionally names for each of the three fields.
%   This GUI will display rolling-binned vertical profiles for up to three
%   fields. When this GUI is closed, it will return a structure of the form:
%
%       struct_out --> field{1,2,3} --> height (in km)
%                                   --> median utc time (sec after midnight)
%                                   --> quality flags (uint16)
%
%   The height returned is the average of the boundary layer heights for
%   that field.  The median UTC time is essentially the middle of the UTC
%   range passed in.
%
%   The quality flags describe the way the BL height was determined, only
%   bits 1-9 are used.  They are grouped in threes, bits 1-3 referring to
%   field 1, 4-6 field 2, and 7-9 field 3.  When set to 1:
%       - Bits 1,4,7 indicate that the boundary layer height from this
%       field was not used.
%       - Bits 2,5,8 indicate that the boundary layer height was modified
%       manually (restricted to a range of bins)
%       - Bits 3,6,9 indicate that multiple heights were averaged for this
%       species
%
% Last Modified by GUIDE v2.5 17-Jul-2014 17:11:05

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @select_BL_heights_OpeningFcn, ...
                   'gui_OutputFcn',  @select_BL_heights_OutputFcn, ...
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


% --- Executes just before select_bl_heights_gui is made visible.
function select_BL_heights_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to select_bl_heights_gui (see VARARGIN)

% Choose default command line output for select_bl_heights_gui
%handles.output = hObject;

% Check the input to the GUI function: there must be between 3 and 5
% inputs, the first must be a structure, the second a vector of UTC times,
% and the 3rd-5th must be field names in the structure
if numel(varargin) < 4; error('BLH_GUI:inputError','Too few inputs, requires a Merge structure, UTC range, altitude field name, and at least one field name'); end
if numel(varargin) > 6; error('BLH_GUI:inputError','Too many inputs, 5 maximum'); end
if ~isstruct(varargin{1}); error('BLH_GUI:inputError','First arugment must be a Merge structure'); end
if ~ismatrix(varargin{2}) || numel(varargin{2})~=2; error('BLH_GUI:inputError','Second argument must be a 2 element row containing beginning and end UTC times'); end
if ~isfield(varargin{1}.Data, varargin{3}); error('BLH_GUI:inputError','Third argument must be a field in the first argument'); end

handles.Merge = varargin{1};
handles.Range = varargin{2};
handles.alt = eval(sprintf('handles.Merge.Data.%s.Values',varargin{3}));
handles.utc = handles.Merge.Data.UTC.Values;
handles.utc_xx = handles.utc >= handles.Range(1) & handles.utc <= handles.Range(2);
draw_map(handles);

% For each of the fields input, check that it is a field of the structure,
% import the data into the handles structure, and prepare the output
% structure within handles. Deactivate any controls relating to unused
% fields, and set the section titles to the appropriate strings
if numel(varargin) >= 4;
    handles.output.field1.name = varargin{4};
    if ~isfield(handles.Merge.Data,handles.output.field1.name); error('BLH_GUI:inputError','%s is not a valid field in the Merge structure',handles.output.field1.name); end
    handles.output.field1.height = NaN;
    handles.output.field1.medianUTC = nanmedian(handles.Range);
    handles.output.field1.qualityFlag = uint8(0);
    set(handles.no2_static,'String',handles.output.field1.name);
    
    field1 = remove_merge_fills(handles.Merge, handles.output.field1.name);
    [fieldbins, fieldalt] = bin_rolling_vertical_profile(handles.alt(handles.utc_xx), field1(handles.utc_xx),0.5,0.1);
    line(fieldbins, fieldalt, 'Parent', handles.no2_axes, 'linewidth', 2, 'color','b', 'marker', 's');
    handles.field1.bins = fieldbins;
    handles.field1.alt = fieldalt;
    handles.field1.include = 1;
    handles.field1.alt_ranges = [];
    handles.field1.heights = [];
    handles.field1.axes_h = handles.no2_axes; % I want the axis handle to be under the fieldx substructure for the plot_bl_heights function
    handles.field1.line_h = [];
    handles.field1.fills = [];
    handles.field1.range_index = [];
    handles.field1.quality_flag = uint8(0);
    
    handles.field1 = plot_bl_heights(handles.field1, handles.output.field1.name);
    
end
if numel(varargin) >= 5;
    handles.output.field2.name = varargin{5};
    if ~isfield(handles.Merge.Data, handles.output.field2.name); error('BLH_GUI:inputError','%s is not a valid field in the Merge structure',handles.output.field2.name); end
    handles.output.field2.height = NaN;
    handles.output.field2.medianUTC = nanmedian(handles.Range);
    handles.output.field2.qualityFlag = uint8(0);
    set(handles.h2o_static,'String',handles.output.field2.name);
    
    field2 = remove_merge_fills(handles.Merge, handles.output.field2.name);
    [fieldbins, fieldalt] = bin_rolling_vertical_profile(handles.alt(handles.utc_xx), field2(handles.utc_xx),0.5,0.1);
    line(fieldbins, fieldalt, 'Parent', handles.h2o_axes, 'linewidth', 2, 'color',[0 0.6 0], 'marker', '^');
    handles.field2.bins = fieldbins;
    handles.field2.alt = fieldalt;
    handles.field2.include = 1;
    handles.field2.alt_ranges = [];
    handles.field2.heights = [];
    handles.field2.axes_h = handles.h2o_axes; % I want the axis handle to be under the fieldx substructure for the plot_bl_heights function
    handles.field2.line_h = [];
    handles.field2.fills = [];
    handles.field2.range_index = [];
    handles.field2.quality_flag = uint8(0);
    
    handles.field2 = plot_bl_heights(handles.field2, handles.output.field2.name);
end
if numel(varargin) >= 6;
    handles.output.field3.name = varargin{6};
    if ~isfield(handles.Merge.Data, handles.output.field3.name); error('BLH_GUI:inputError','%s is not a valid field in the Merge structure',handles.output.field3.name); end
    handles.output.field3.height = NaN;
    handles.output.field3.medianUTC = nanmedian(handles.Range);
    handles.output.field3.qualityFlag = uint8(0);
    set(handles.theta_static,'String',handles.output.field3.name);
    
    field3 = remove_merge_fills(handles.Merge, handles.output.field3.name);
    [fieldbins, fieldalt] = bin_rolling_vertical_profile(handles.alt(handles.utc_xx), field3(handles.utc_xx),0.5,0.1);
    line(fieldbins, fieldalt, 'Parent', handles.theta_axes, 'linewidth', 2, 'color','r', 'marker', 'o');
    handles.field3.bins = fieldbins;
    handles.field3.alt = fieldalt;
    handles.field3.include = 1;
    handles.field3.alt_ranges = [];
    handles.field3.heights = [];
    handles.field3.axes_h = handles.theta_axes; % I want the axis handle to be under the fieldx substructure for the plot_bl_heights function
    handles.field3.line_h = [];
    handles.field3.fills = [];
    handles.field3.range_index = [];
    handles.field3.quality_flag = uint8(0);
    
    handles.field3 = plot_bl_heights(handles.field3, handles.output.field3.name);
end

if numel(varargin) < 5; handles = deactivate_field(handles,2); end
if numel(varargin) < 6; handles = deactivate_field(handles,3); end

% Record the number of fields to be averaged
handles.numfield = numel(varargin)-3;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes select_bl_heights_gui wait for user response (see UIRESUME)
uiwait(handles.select_bl_heights_GUI);


% --- Outputs from this function are returned to the command line.
function varargout = select_BL_heights_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;
delete(hObject);


% --- Executes on button press in no2reject_button.
function no2reject_button_Callback(hObject, eventdata, handles)
% hObject    handle to no2reject_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.field1 = flip_state(handles.field1, hObject);
guidata(hObject, handles);

% --- Executes on button press in no2blh_button.
function no2blh_button_Callback(hObject, eventdata, handles)
% hObject    handle to no2blh_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
fieldname = handles.output.field1.name;
if strcmpi(fieldname,'THETA'); fieldname = 'THETALITE'; end
handles.field1 = plot_bl_heights(handles.field1,fieldname);
guidata(hObject,handles);

% --- Executes on button press in no2AddRange_button.
function no2AddRange_button_Callback(hObject, eventdata, handles)
% hObject    handle to no2AddRange_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.field1 = add_range(handles.field1);
guidata(hObject,handles)

% --- Executes on button press in no2DeleteRange_button.
function no2DeleteRange_button_Callback(hObject, eventdata, handles)
% hObject    handle to no2DeleteRange_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.field1 = delete_range(handles.field1);
guidata(hObject,handles);

% --- Executes on button press in no2MoveUp_button.
function no2MoveUp_button_Callback(hObject, eventdata, handles)
% hObject    handle to no2MoveUp_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.field1 = move_range_up(handles.field1);
guidata(hObject, handles);

% --- Executes on button press in no2Widen_button.
function no2Widen_button_Callback(hObject, eventdata, handles)
% hObject    handle to no2Widen_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.field1 = widen_range(handles.field1);
guidata(hObject,handles);

% --- Executes on button press in no2Shrink_button.
function no2Shrink_button_Callback(hObject, eventdata, handles)
% hObject    handle to no2Shrink_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.field1 = shrink_range(handles.field1);
guidata(hObject,handles);

% --- Executes on button press in no2MoveDown_button.
function no2MoveDown_button_Callback(hObject, eventdata, handles)
% hObject    handle to no2MoveDown_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.field1 = move_range_down(handles.field1);
guidata(hObject,handles);

% --- Executes on button press in h2oAddRange_button.
function h2oAddRange_button_Callback(hObject, eventdata, handles)
% hObject    handle to h2oAddRange_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.field2 = add_range(handles.field2);
guidata(hObject,handles)

% --- Executes on button press in h2oDeleteRange_button.
function h2oDeleteRange_button_Callback(hObject, eventdata, handles)
% hObject    handle to h2oDeleteRange_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.field2 = delete_range(handles.field2);
guidata(hObject,handles);

% --- Executes on button press in h2oMoveUp_button.
function h2oMoveUp_button_Callback(hObject, eventdata, handles)
% hObject    handle to h2oMoveUp_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.field2 = move_range_up(handles.field2);
guidata(hObject,handles);

% --- Executes on button press in h2oWiden_button.
function h2oWiden_button_Callback(hObject, eventdata, handles)
% hObject    handle to h2oWiden_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.field2 = widen_range(handles.field2);
guidata(hObject,handles)

% --- Executes on button press in h2oShrink_button.
function h2oShrink_button_Callback(hObject, eventdata, handles)
% hObject    handle to h2oShrink_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.field2 = shrink_range(handles.field2);
guidata(hObject,handles);

% --- Executes on button press in h2oMoveDown_button.
function h2oMoveDown_button_Callback(hObject, eventdata, handles)
% hObject    handle to h2oMoveDown_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.field2 = move_range_down(handles.field2);
guidata(hObject, handles);

% --- Executes on button press in thetaAddRange_button.
function thetaAddRange_button_Callback(hObject, eventdata, handles)
% hObject    handle to thetaAddRange_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.field3 = add_range(handles.field3);
guidata(hObject,handles)

% --- Executes on button press in thetaDeleteRange_button.
function thetaDeleteRange_button_Callback(hObject, eventdata, handles)
% hObject    handle to thetaDeleteRange_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.field3 = delete_range(handles.field3);
guidata(hObject,handles);

% --- Executes on button press in thetaMoveUp_button.
function thetaMoveUp_button_Callback(hObject, eventdata, handles)
% hObject    handle to thetaMoveUp_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.field3 = move_range_up(handles.field3);
guidata(hObject, handles);

% --- Executes on button press in thetaWiden_button.
function thetaWiden_button_Callback(hObject, eventdata, handles)
% hObject    handle to thetaWiden_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.field3 = widen_range(handles.field3);
guidata(hObject, handles);

% --- Executes on button press in thetaShrink_button.
function thetaShrink_button_Callback(hObject, eventdata, handles)
% hObject    handle to thetaShrink_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.field3 = shrink_range(handles.field3);
guidata(hObject, handles);

% --- Executes on button press in thetaMoveDown_button.
function thetaMoveDown_button_Callback(hObject, eventdata, handles)
% hObject    handle to thetaMoveDown_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.field3 = move_range_down(handles.field3);
guidata(hObject,handles);

% --- Executes on button press in thetareject_button.
function thetareject_button_Callback(hObject, eventdata, handles)
% hObject    handle to thetareject_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.field3 = flip_state(handles.field3, hObject);
guidata(hObject,handles);

% --- Executes on button press in thetablh_button.
function thetablh_button_Callback(hObject, eventdata, handles)
% hObject    handle to thetablh_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
fieldname = handles.output.field3.name;
if strcmpi(fieldname,'THETA'); fieldname = 'THETALITE'; end
handles.field3 = plot_bl_heights(handles.field3,fieldname);
guidata(hObject,handles);

% --- Executes on button press in h2oreject_button.
function h2oreject_button_Callback(hObject, eventdata, handles)
% hObject    handle to h2oreject_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.field2 = flip_state(handles.field2, hObject);
guidata(hObject,handles);

% --- Executes on button press in h2oblh_button.
function h2oblh_button_Callback(hObject, eventdata, handles)
% hObject    handle to h2oblh_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
fieldname = handles.output.field2.name;
if strcmpi(fieldname,'THETA'); fieldname = 'THETALITE'; end
handles.field2 = plot_bl_heights(handles.field2,fieldname);
guidata(hObject,handles);

% --- Executes on mouse press over axes background.
function no2_axes_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to no2_axes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

pt = get(handles.field1.axes_h,'CurrentPoint');
% Find the range that was just clicked on
if ~isempty(handles.field1.alt_ranges); ri = find(pt(1,2) >= handles.field1.alt_ranges(:,1) & pt(1,2) <= handles.field1.alt_ranges(:,2)); 
else ri = []; 
end
% Multiple ranges might meet the above criteria; select the first one
if numel(ri) > 1; ri = ri(1); end
handles.field1.range_index = ri;
% Refresh the drawing of ranges
handles.field1 = redraw_ranges(handles.field1);
% update the handles structure
guidata(hObject,handles);


% --- Executes on mouse press over axes background.
function h2o_axes_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to h2o_axes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

pt = get(handles.field2.axes_h,'CurrentPoint');
% Find the range that was just clicked on
if ~isempty(handles.field2.alt_ranges); ri = find(pt(1,2) >= handles.field2.alt_ranges(:,1) & pt(1,2) <= handles.field2.alt_ranges(:,2)); 
else ri = []; 
end
% Multiple ranges might meet the above criteria; select the first one
if numel(ri) > 1; ri = ri(1); end
handles.field2.range_index = ri;
% Refresh the drawing of ranges
handles.field2 = redraw_ranges(handles.field2);
% update the handles structure
guidata(hObject,handles);

% --- Executes on mouse press over axes background.
function theta_axes_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to theta_axes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

pt = get(handles.field3.axes_h,'CurrentPoint');
% Find the range that was just clicked on
if ~isempty(handles.field3.alt_ranges); ri = find(pt(1,2) >= handles.field3.alt_ranges(:,1) & pt(1,2) <= handles.field3.alt_ranges(:,2)); 
else ri = []; 
end
% Multiple ranges might meet the above criteria; select the first one
if numel(ri) > 1; ri = ri(1); end
handles.field3.range_index = ri;
% Refresh the drawing of ranges
handles.field3 = redraw_ranges(handles.field3);
% update the handles structure
guidata(hObject,handles);

% --- Executes on button press in done_button.
function done_button_Callback(hObject, eventdata, handles)
% hObject    handle to done_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Output field-specific values for whatever fields exist
handles.output.field1.height = nanmean(handles.field1.heights);
handles.output.field1.qualityFlag = handles.field1.quality_flag;
handles.output.field1.include = handles.field1.include;
if handles.output.field1.include 
    allheights(1) = handles.output.field1.height;
else
    allheights(1) = NaN;
end
allquality(1) = uint16(handles.output.field1.qualityFlag);

if handles.numfield > 1;
    handles.output.field2.height = nanmean(handles.field2.heights);
    handles.output.field2.qualityFlag = handles.field2.quality_flag;
    handles.output.field2.include = handles.field2.include;
    if handles.output.field2.include 
        allheights(2) = handles.output.field2.height;
    else
        allheights(2) = NaN;
    end
    allquality(2) = bitshift(uint16(handles.output.field2.qualityFlag),3);
else
    allquality(2) = uin16(0);
end
if handles.numfield > 2;
    handles.output.field3.height = nanmean(handles.field3.heights);
    handles.output.field3.qualityFlag = handles.field3.quality_flag;
    handles.output.field3.include = handles.field3.include;
    if handles.output.field3.include 
        allheights(3) = handles.output.field3.height;
    else
        allheights(3) = NaN;
    end
    allquality(3) = bitshift(uint16(handles.output.field3.qualityFlag),6);
else
    allquality(3) = uint16(0);
end

handles.output.Overall.height = nanmean(allheights);
handles.output.Overall.medianUTC = handles.output.field1.medianUTC;
handles.output.Overall.qualityFlag = bitor(allquality(3), bitor(allquality(2), allquality(1)));


guidata(hObject, handles);
select_bl_heights_GUI_CloseRequestFcn(handles.select_bl_heights_GUI, eventdata, handles);

% --- Executes on button press in help_button.
function help_button_Callback(hObject, eventdata, handles)
% hObject    handle to help_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

homedir = userpath; homedir = homedir(1:end-1);
helptext = fileread(fullfile(homedir,'NO2 Profiles','GUIs','select_BL_heights_help.txt'));
helpdlg(helptext,'Select BL Heights GUI Help');

% --- Executes when user attempts to close select_bl_heights_GUI.
function select_bl_heights_GUI_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to select_bl_heights_GUI (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: delete(hObject) closes the figure
if isequal(get(handles.select_bl_heights_GUI,'waitstatus'),'waiting')
uiresume(handles.select_bl_heights_GUI);
else
delete(hObject);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% ADDITIONAL FCNS %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%
% Make unused field unusable
%%%%%
function S = deactivate_field(S,fieldnum)
names = {'no2','h2o','theta'};
fieldname = names{fieldnum}; h = zeros(1,10);
h(1) = eval(sprintf('S.%s_axes',fieldname));
h(2) = eval(sprintf('S.%sAddRange_button',fieldname));
h(3) = eval(sprintf('S.%sDeleteRange_button',fieldname));
h(4) = eval(sprintf('S.%sWiden_button',fieldname));
h(5) = eval(sprintf('S.%sShrink_button',fieldname));
h(6) = eval(sprintf('S.%sMoveUp_button',fieldname));
h(7) = eval(sprintf('S.%sMoveDown_button',fieldname));
h(8) = eval(sprintf('S.%sAddRange_button',fieldname));
h(9) = eval(sprintf('S.%sreject_button',fieldname));
h(10) = eval(sprintf('S.%sblh_button',fieldname));

set(h(1),'Color',[0.3 0.3 0.3]);
for b=2:10
    set(h(b),'Enable','off');
end

%%%%%
% Plot the boundary layer heights for the given axis
%%%%%
function S_field = plot_bl_heights(S_field, fieldname)
% S_field is handles.fieldx, alt_ranges are the altitude ranges (n x 2
% matrix) to look in for boundary layers, fieldname is usually the data
% field name in the Merge structure, but will actually be what is passed to
% the "find_bdy_layer_height" function as the method to use.
bins = S_field.bins; alts = S_field.alt; 
alt_ranges = S_field.alt_ranges; axes_handle = S_field.axes_h;
% Clear any existing BL heights in the data structure and the axes
if ~isempty(S_field.line_h)
    delete(S_field.line_h);
    S_field.line_h = [];
    S_field.heights = [];
end
heights = -127*ones(size(alt_ranges,1),1);
lines = gobjects(size(alt_ranges,1),1);
if isempty(alt_ranges); 
    alt_ranges = [-Inf, Inf]; 
    S_field.quality_flag = bitset(S_field.quality_flag,2,0);
else
    S_field.quality_flag = bitset(S_field.quality_flag,2,1);
end
for b=1:size(alt_ranges,1)
    rr = alts >= alt_ranges(b,1) & alts <= alt_ranges(b,2);
    heights(b) = find_bdy_layer_height(bins(rr), alts(rr), fieldname);
    xlims = get(axes_handle,'xlim');
    lines(b) = line(xlims,[heights(b), heights(b)],'Parent',axes_handle, 'linewidth',2,'color','k','linestyle','--');
end
S_field.line_h = lines; S_field.heights = heights;
if numel(S_field.heights) > 1; %If more than 1 height is evaluated, set the warning bit
    S_field.quality_flag = bitset(S_field.quality_flag,3,1);
else %If only one (or none) heights, reset the warning bit to 0
    S_field.quality_flag = bitset(S_field.quality_flag,3,0);
end
    

%%%%%
% Change the state of the profile to be included or not
%%%%%
function S_field = flip_state(S_field, button_handle)
btn_text = {'Accept','Reject'};
btn_tooltip = {'Include this profile to find a BL height','Do not use this profile to find a BL height'};
S_field.include = ~S_field.include;
ind = S_field.include + 1;
set(button_handle,'String',btn_text{ind});
set(button_handle,'TooltipString',btn_tooltip{ind});
S_field.quality_flag = bitset(S_field.quality_flag,1,~S_field.include);

%%%%%
% Redraw the ranges for the given axes
%%%%%
function S_field = redraw_ranges(S_field)

% Remove all existing fills
if ~isempty(S_field.fills);
    delete(S_field.fills)
    S_field.fills = [];
end

ranges = S_field.alt_ranges;
plot_xlim = get(S_field.axes_h,'xlim');
for a=1:size(ranges,1)
    if a == S_field.range_index
        fillcol = 'b';
    else
        fillcol = 'r';
    end
    fill_y = [ranges(a,1), ranges(a,2), ranges(a,2), ranges(a,1)];
    fill_x = [plot_xlim(1), plot_xlim(1), plot_xlim(2), plot_xlim(2)];
    fh = patch(fill_x, fill_y, fillcol,'FaceAlpha',0.4,'PickableParts','none','Parent',S_field.axes_h);
    S_field.fills = [S_field.fills, fh];
end

%%%%%
% Add a new range to the axes
%%%%%
function S_field = add_range(S_field)
new_range = [min(S_field.alt(:))-0.05,min(S_field.alt(:))+0.25];
S_field.alt_ranges = [S_field.alt_ranges; new_range];
S_field.range_index = size(S_field.alt_ranges,1);
S_field = redraw_ranges(S_field);


%%%%%
% Delete the currently selected range
%%%%%
function S_field = delete_range(S_field)
x=5;
if ~isempty(S_field.range_index);
    ri = S_field.range_index;
    delete(S_field.fills(ri));
    S_field.fills(ri) = [];
    S_field.alt_ranges(ri,:) = [];
    S_field.range_index = [];
end

%%%%%
% Expand the currently selected range
%%%%%
function S_field = widen_range(S_field)
if ~isempty(S_field.range_index);
    ri = S_field.range_index;
    sr = S_field.alt_ranges(ri,:);
    ylims = get(S_field.axes_h,'ylim');
    if sr(2) < max(ylims) % Don't let the user expand (too far) off screen
        sr(2) = sr(2) + 0.1; 
        S_field.alt_ranges(ri,:) = sr;
        S_field = redraw_ranges(S_field);
    end
end

%%%%%
% Shrink the currently selected range
%%%%%
function S_field = shrink_range(S_field)
if ~isempty(S_field.range_index);
    ri = S_field.range_index;
    sr = S_field.alt_ranges(ri,:);
    if sr(2) - sr(1) > 0.2 %Keep the user from shrinking the range too much
        sr(2) = sr(2) - 0.1;
        S_field.alt_ranges(ri,:) = sr;
        S_field = redraw_ranges(S_field);
    end
end

%%%%%
% Move the currently selected range up
%%%%%
function S_field = move_range_up(S_field)
if ~isempty(S_field.range_index);
    ri = S_field.range_index;
    sr = S_field.alt_ranges(ri,:);
    ylims = get(S_field.axes_h,'ylim');
    if sr(2) + 0.1 < max(ylims) % Don't let the user move the range off screen
        sr = sr + 0.1;
        S_field.alt_ranges(ri,:) = sr;
        S_field = redraw_ranges(S_field);
    end
end

%%%%%
% Move the currently selected range down
%%%%%
function S_field = move_range_down(S_field)
if ~isempty(S_field.range_index);
    ri = S_field.range_index;
    sr = S_field.alt_ranges(ri,:);
    ylims = get(S_field.axes_h,'ylim');
    if sr(1) - 0.1 > min(ylims) % Don't let the user move the range off screen
        sr = sr - 0.1;
        S_field.alt_ranges(ri,:) = sr;
        S_field = redraw_ranges(S_field);
    end
end

%%%%%
% Draw the map!
%%%%%
function draw_map(S)
lat = S.Merge.Data.LATITUDE.Values;
lon = S.Merge.Data.LONGITUDE.Values - 360;
alt = S.alt;
xx = S.utc_xx;

usa = shaperead('usastatehi.shp');
for b=1:numel(usa);
    if any(strcmpi(usa(b).Name,{'Alaska','Hawaii'}));
    else
        line(usa(b).X, usa(b).Y, 'color','k','parent',S.map_axes);
    end
end

hold(S.map_axes);
scatter(lon(xx), lat(xx), 16, alt(xx), 'parent', S.map_axes);
hold(S.map_axes);
lonlim = [floor(min(lon(:))), ceil(max(lon(:)))];
latlim = [floor(min(lat(:))), ceil(max(lat(:)))];
set(S.map_axes,'xlim',lonlim); set(S.map_axes,'ylim',latlim);
cb=colorbar(S.map_axes); ylabel(cb,'Altitude (km)');



