function varargout = aerosol_no2_profiles_gui(varargin)
% AEROSOL_NO2_PROFILES_GUI MATLAB code for aerosol_no2_profiles_gui.fig
%      AEROSOL_NO2_PROFILES_GUI, by itself, creates a new AEROSOL_NO2_PROFILES_GUI or raises the existing
%      singleton*.
%
%      H = AEROSOL_NO2_PROFILES_GUI returns the handle to a new AEROSOL_NO2_PROFILES_GUI or the handle to
%      the existing singleton*.
%
%      AEROSOL_NO2_PROFILES_GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in AEROSOL_NO2_PROFILES_GUI.M with the given input arguments.
%
%      AEROSOL_NO2_PROFILES_GUI('Property','Value',...) creates a new AEROSOL_NO2_PROFILES_GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before aerosol_no2_profiles_gui_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to aerosol_no2_profiles_gui_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help aerosol_no2_profiles_gui

% Last Modified by GUIDE v2.5 03-Apr-2015 11:19:46

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @aerosol_no2_profiles_gui_OpeningFcn, ...
                   'gui_OutputFcn',  @aerosol_no2_profiles_gui_OutputFcn, ...
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


% --- Executes just before aerosol_no2_profiles_gui is made visible.
function aerosol_no2_profiles_gui_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to aerosol_no2_profiles_gui (see VARARGIN)


%%%%% MY STUFF %%%%%

% Require two tables to be passed as arguments, otherwise throw an error
E = JLLErrors;
if numel(varargin) ~= 2 || any(~iscellcontents(varargin,'istable'))
    E.badinput('Exactly 2 tables must be passed to this function');
end
handles = load_tables(handles,varargin);

% Initialize the field which will describe which row is currently selected
% in the table. We'll initialize it to 0, and test later that this number
% is > 0. This will be set in cat_table_CellSelectionCallback
handles.curr_row = 0;

%%%% GUI INNARDS %%%%%

% Choose default command line output for aerosol_no2_profiles_gui
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes aerosol_no2_profiles_gui wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = aerosol_no2_profiles_gui_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes when selected cell(s) is changed in cat_table.
function cat_table_CellSelectionCallback(hObject, eventdata, handles)
% hObject    handle to cat_table (see GCBO)
% eventdata  structure with the following fields (see MATLAB.UI.CONTROL.TABLE)
%	Indices: row and column indices of the cell(s) currently selecteds
% handles    structure with handles and user data (see GUIDATA)

handles.curr_row = eventdata.Indices(1);
guidata(hObject,handles);


% --- Executes on selection change in campaign_menu.
function campaign_menu_Callback(hObject, eventdata, handles)
% hObject    handle to campaign_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns campaign_menu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from campaign_menu

handles.curr_campaign = get(hObject,'Value');
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function campaign_menu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to campaign_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% Also set the possible campaigns here. The strings here should be of the
% same format that merge_field_names is expecting.

handles.campaigns = {'DISCOVER-MD','DISCOVER-CA','DISCOVER-TX','DISCOVER-CO'};
handles.curr_campaign = 1;

% Initialize the popup menu
set(hObject,'String',handles.campaigns);
set(hObject,'Value',handles.curr_campaign);

guidata(hObject,handles);



% --- Executes on button press in make_plot_btn.
function make_plot_btn_Callback(hObject, eventdata, handles)
% hObject    handle to make_plot_btn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%fprintf('Will plot for campaign %s\n',handles.campaigns{handles.curr_campaign});
if handles.curr_row > 0
    plot_curr_profile(handles);
else
    fprintf('No row currently selected\n')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% ADDITIONAL FUNCTIONS %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function S = load_tables(S,tables)
set(S.cat_table,'data',table2cell(tables{1}));
set(S.cat_table,'columnname',tables{1}.Properties.VariableNames);
set(S.cat_table,'rowname',tables{1}.Properties.RowNames);

set(S.crit_table,'data',table2cell(tables{2}));
set(S.crit_table,'columnname',tables{2}.Properties.VariableNames);
set(S.crit_table,'rowname',tables{2}.Properties.RowNames);

function plot_curr_profile(S)
campaign_name = S.campaigns{S.curr_campaign};
[Names,dates,mdir] = merge_field_names(campaign_name);
Data = get(S.cat_table,'Data');
r = S.curr_row;
sel_profnum = Data{r,1};
sel_date = datenum(Data{r,2});

% Check that the date makes sense
if sel_date < min(datenum(dates)) || sel_date > max(datenum(dates))
    fprintf('Date selected outside dates of chosen campaign.\n')
    return
end

file_pat = strcat(upper(regexprep(campaign_name,'\W','_')),'*%s_%s_%s.mat');
F = wildcard_load(mdir,file_pat,sel_date,'Merge');

if numel(sel_profnum) == 1
    profnums = remove_merge_fills(F.Merge,Names.profile_numbers);
else
    utc = remove_merge_fills(F.Merge, 'UTC');
end
no2 = remove_merge_fills(F.Merge,Names.no2_lif);
aer = remove_merge_fills(F.Merge,Names.aerosol_extinction);
alt = remove_merge_fills(F.Merge,Names.gps_alt);

if numel(sel_profnum) == 1
    xx = profnums == sel_profnum;
else
    xx = utc >= sel_profnum(1) & utc < sel_profnum(2);
end
prof_no2 = no2(xx);
prof_no2_alt = alt(xx);
prof_aer = aer(xx);
prof_aer_alt = alt(xx);

[no2_bins, no2_bin_alt] = bin_vertical_profile(prof_no2_alt, prof_no2, 0.25);
[aer_bins, aer_bin_alt] = bin_vertical_profile(prof_aer_alt, prof_aer, 0.25);

figure;
[hax,h1,h2] = plotxx(no2_bins,no2_bin_alt,aer_bins,aer_bin_alt,{'NO2','Aerosol'},{'Alt/NO2','Alt/Aerosol'});
set(h1,'Linewidth',5);
set(h2,'Linewidth',5);
ylim1 = hax(1).YLim; ylim2 = hax(2).YLim;
newylim = [0, max(ylim1(2),ylim2(2))];
set(hax(1),'ylim',newylim);
set(hax(2),'ylim',newylim);

[no2_bin_alt, no2_bins] = fill_nans(no2_bin_alt, no2_bins);
[aer_bin_alt, aer_bins] = fill_nans(aer_bin_alt, aer_bins);

line(no2_bins,no2_bin_alt,'parent',hax(1),'color','b','linewidth',2);
line(aer_bins,aer_bin_alt,'parent',hax(2),'color',[0 0.5 0],'linewidth',2);

