function [ ] = ncatts( ni, varname )
%NCATTS Prints out the value of each attribute for a variable in a netCDF file
%   NCATTS( NI, VARNAME ) prints out the name and value for each attribute
%   in the variable VARNAME within the structure NI returned from NCINFO.
%
%   NCATTS( NI ) or NCATTS( NI, '/' ) prints the root attributes.

E = JLLErrors;

if ~isstruct(ni) || any(~isfield(ni,{'Filename','Name','Dimensions','Variables','Attributes'}))
    E.badinput('NI does not appear to be a structure returned from NCINFO')
end

if ~exist('varname','var')
    varname = '/';
elseif ~ischar(varname)
    E.badinput('VARNAME must be a string');
end

if strcmp(varname, '/')
    print_atts(ni)
else
    for a = 1:numel(ni.Variables)
        if strcmp(ni.Variables(a).Name, varname)
            print_atts(ni.Variables(a));
            return
        end
    end
    fprintf('Variable "%s" not found.\n', varname);
end
end

function print_atts(ncdf_var)
for b=1:numel(ncdf_var.Attributes)
    name = ncdf_var.Attributes(b).Name;
    val = ncdf_var.Attributes(b).Value;
    if isnumeric(val)
        if mod(val,1) == 0
            fprintf('%s: %d\n', name, val);
        else
            fprintf('%s: %g\n', name, val);
        end
    else
        fprintf('%s: %s\n', name, val);
    end
end
end