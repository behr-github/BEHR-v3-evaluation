function  geos_ncdf_monthly_avg(input_ncdfs, save_dir, varargin)
%GEOS_NCDF_MONTHLY_AVG Compute a monthly average of ND51 netCDF outputs
%   GEOS_NCDF_MONTHLY_AVG( INPUT_NCDFS, SAVE_DIR ) Takes a list of netCDF
%   files (INPUT_NCDFS) as either a cell array or structure returned by
%   DIRFF() and calculates monthly averages of any variables that have a
%   "time" dimension. This will be a simple average, no weighting applied.
file_names = files_input(input_ncdfs);

unique_months = get_all_months(file_names);
for i_month = 1:numel(unique_months)
    month_files = list_files_for_months(unique_months(i_month), file_names);
    average_files(month_files, save_dir, unique_months(i_month));
end

end

function average_files(files, save_dir, month_dnum)
% Load the first file, figure out which variables we need to average and
% which need to be checked that they are identical. Also generate the save
% name and copy the schema.
E = JLLErrors;
fprintf('Averaging %s\n', datestr(month_dnum, 'mmm yyyy'));
gci = ncinfo(files{1});

save_name = regexprep(files{1}, '\d{8}', sprintf('%savg', datestr(month_dnum, 'yyyymm')));
save_name = fullfile(save_dir, save_name);
ncwriteschema(save_name, gci);

var_names = {gci.Variables.Name};
[do_average_var, var_field_names] = make_empty_struct_from_cell(var_names, false);
var_running_avg = make_empty_struct_from_cell(var_field_names, []);
% If a variable has a "time" dimension, then it will be averaged. Mark it
% as true in the "do_average_var" struct. Variables marked as "false" will
% be required to be identical in all files.
for i_var = 1:numel(var_names)
    if isempty(gci.Variables(i_var).Dimensions)
        avg_this_var = false;
    else
        var_dims = {gci.Variables(i_var).Dimensions.Name};
        xx = strcmpi(var_dims, 'time');
        avg_this_var = sum(xx) == 1;
        if sum(xx) > 1
            E.callError('multiple_time_dims', 'Found multiple "time" dimensions in variable %s', var_names{i_var});
        end
    end
    do_average_var.(var_field_names{i_var}) = avg_this_var;
    % Initialize an array for the running average for averaged variables;
    % read in the non averaged variables now. That will make the loop over
    % files cleaner by avoiding a check if we're in the first iteration.
    if avg_this_var
        var_running_avg.(var_field_names{i_var}) = zeros(gci.Variables(i_var).Size);
    else
        var_running_avg.(var_field_names{i_var}) = ncread(gci.Filename, var_names{i_var});
    end
end

for i_file = 1:numel(files)
    for i_var = 1:numel(var_names)
        var_value = ncread(files{i_file}, var_names{i_var});
        if do_average_var.(var_field_names{i_var})
            var_running_avg.(var_field_names{i_var}) = var_running_avg.(var_field_names{i_var}) + var_value;
        elseif ~isequal(var_value, var_running_avg.(var_field_names{i_var}))
            E.callError('Variable "%s" is not the same in %s as in the first file (%s)', var_names{i_var}, files{i_file}, files{1});
        end
    end
end

for i_var = 1:numel(var_names)
    if do_average_var.(var_field_names{i_var})
        n = numel(files);
    else
        n = 1;
    end
    ncwrite(save_name, var_names{i_var}, var_running_avg.(var_field_names{i_var}) / n);
end

end

function month_files = list_files_for_months(month_dnum, file_names)
E = JLLErrors;
file_dates = get_file_datenums(file_names);
xx = year(file_dates) == year(month_dnum) & month(file_dates) == month(month_dnum);
n_days = eomday(year(month_dnum), month(month_dnum));
if sum(xx) ~= n_days
    E.callError('insufficient_files_for_month', 'Wrong number of files (%d instead of %d) identified for %s', sum(xx), n_days, datestr(month_dnum, 'mmm yyyy'));
end
month_files = file_names(xx);
end

function month_dnums = get_all_months(file_names)
% Get all unique months that exist in the input files.
file_dates = get_file_datenums(file_names);

year_month = [year(file_dates(:)), month(file_dates(:))];
year_month = unique(year_month, 'rows');
month_dnums = nan(size(year_month,1),1);
for a=1:numel(month_dnums)
    month_dnums(a) = datenum(year_month(a,1), year_month(a,2), 1);
end

end

function file_dates = get_file_datenums(file_names)
file_dates = nan(size(file_names));
for a=1:numel(file_names)
    [~, base_name] = fileparts(file_names{a});
    file_dates(a) = datenum(regexp(base_name, '\d{8}','match','once'),'yyyymmdd');
end
end