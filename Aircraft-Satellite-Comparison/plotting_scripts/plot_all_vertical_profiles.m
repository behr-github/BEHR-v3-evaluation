function [  ] = plot_all_vertical_profiles( Merge, data_field, binwidth, varargin )
%plot_all_vertical_profiles(Merge, data_field, binwidth, [opt] groupmode): Plots all vertical profiles for "data_field" in "Merge" binned to "binwidth" (km) layers
%   This function will generate plots of the vertical profile for the
%   variable described in "data_field" for all unique non-zero identifiers in either
%   profile number or site flag. There are three options for this argument:
%
%       'profnum' = Group measurments by raw profile number. Active by
%          default.
%       'profnumgroup' = Group measurements by the leading digit of the
%          profile number. That is, 1001, 1002, and 1003 would all be
%          grouped together. ***CURRENTLY NOT IMPLEMENTED***
%       'siteflag' = Group measurements by their geographic site flag.
%
%   There is also a parameter value, 'opmode' which can be set to 'rolling'
%   to make this function use rolling bins instead of normal ones.

p = inputParser;
p.addRequired('Merge',@isstruct);
p.addRequired('data_field', @(x) ischar(x) || (iscell(x) && length(x)<=2));
p.addRequired('binwidth',@isscalar);
p.addOptional('groupmode','profnum',@(x) any(strcmpi(x,{'profnum','siteflag','profnumgroup'})))
p.addParamValue('opmode','std', @(x) any(strcmpi(x,{'rolling','standard','std','std.'})));

p.parse(Merge,data_field, binwidth, varargin{:});
pout = p.Results;

Merge = pout.Merge;
binwidth = pout.binwidth;
groupmode = pout.groupmode;
opmode = pout.opmode;

if iscell(pout.data_field)
    field = pout.data_field{1};
    field2 = pout.data_field{2};
else
    field = pout.data_field;
end

% Set the criteria by which we will group measurements based on the
% optional fourth argument
if strcmpi(groupmode,'siteflag')
    groupids = Merge.Data.discoveraqSiteFlag1sec.Values;
    groupmethod = 'siteflag';
% elseif strcmpi(groupmode,'profnumgroup');
%     groupids = Merge.Data.ProfileSequenceNum.Values;
%     groupids = floor(groupids/1000);
%     groupmethod = 'profnum';
else
    groupids = Merge.Data.ProfileSequenceNum.Values;
    groupmethod = 'profnum';
end

% Find all unique, non-zero entries in our grouping criteria
id_vals = unique(groupids);
id_vals = id_vals(id_vals>0);

% For each unique value, plot the vertical profile.

for a = 1:numel(id_vals)
    plot_vertical_profile(Merge,field,'binwidth',binwidth,groupmethod,id_vals(a),'opmode',opmode);
    
    
    if exist('field2','var') %If there was a second field to plot, do so now.
        pos = get(gca,'position');
        h2 = axes('position',pos,'YAxisLocation','right','color','none','XTick','none');
        hold on
        plot_vertical_profile(Merge,field2,'binwidth',binwidth2,groupmethod,id_vals(a),'axes',h2);
    end
end

end

