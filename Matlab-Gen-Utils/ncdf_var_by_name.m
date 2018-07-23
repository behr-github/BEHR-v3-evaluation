function  var_struct  = ncdf_var_by_name(info_struct, var_name)
%   VAR_STRUCT = NCDF_VAR_BY_NAME(INFO_STRUCT, VAR_NAME) Finds the variable
%   named VAR_NAME in the netCDF schema INFO_STRUCT (i.e. a structure
%   returned from NCINFO) and returns the variable's structure as
%   VAR_STRUCT.
%
%   Currently this only searches top level variables; searching recursively
%   through netCDF4 groups is not implemented.

if ~isstruct(info_struct) || ~isfield(info_struct, 'Variables') || ~isstruct(info_struct.Variables) || ~isfield(info_struct.Variables, 'Name')
    E.badinput('INFO_STRUCT does not appear to be a netCDF schema struct (either does not have "Variables" field or the field is not a struct with the field "Name")');
end
if ~ischar(var_name)
    E.badinput('VAR_NAME must be a char array')
end

all_var_names = {info_struct.Variables.Name};
xx = strcmp(var_name, all_var_names);
if sum(xx) == 1
    var_struct = info_struct.Variables(xx);
elseif sum(xx) < 1
    E.callError('var_not_found', 'Could not find variable named "%s" in %s', var_name, info_struct.Filename);
else
    E.callError('multiple_var_matches', 'Multiple variables names matched "%s" in %s', var_name, info_struct.Filename);
end

end

