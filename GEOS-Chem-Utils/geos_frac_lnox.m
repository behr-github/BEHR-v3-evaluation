function [ frac_lnox, frac_surface, gc_datevec, Prod_Amt ] = geos_frac_lnox( prod_struct )
%geos_percent_lnox Calculate percent of NO production due to lightning.
%   When considering the impact of MPN and slower HNO3 production, it is
%   useful to find GEOS-Chem grid cells that have a certain percent of the
%   NOx produced due to lightning.  This function will calculate that for
%   each time entry in an input structure produced by read_geos_output.
%   This structure is the only input required to this function.
%
%   This returns 4 variables: 
%       frac_lnox - the fraction (in decimal) of NO production due to
%       lightning in each grid cell for each time entry.
%
%       frac_surface - the fraction of NO produced due to surface sources
%       (anthropogenic, biomass burning, biofuel, fertilizer, and soil).
%
%       gc_datevec - a vector a datenums corresponding to the time entries
%       in the input structure.
%
%       Prod_Amt - a structure containing the summed production for each
%       category, plus total NO production (this will be useful if you wish
%       to examine each source independently)
%
%   This function expects the input structure to have 8 fields.  The
%   structure must be of the form produces by read_geos_output.  The 8
%   fields must include the following fullCat names:
%       1) 'Aircraft NO'
%       2) 'Anthropogenic NO'
%       3) 'Biomass NO'
%       4) 'Biofuel NO'
%       5) 'Fertilizer NO'
%       6) 'Lightning NO'
%       7) 'Soil NO'
%       8) 'Stratospheric NO'
%   

E = JLLErrors;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% INPUT VALIDATION %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% First make sure that this is a structure output from read_geos_output
if ~check_struct_type(prod_struct)
    error(E.badinput('Input prod_struct is not an output from read_geos_output (one or more expected fields is missing)'));
end

% Second check that all the GEOS-Chem output fields we expect to be in the
% input structure are present. Add categories here if needed.  Be sure to
% update isLNOx and isSurfaceNOx too.

categories = {'Aircraft NO','Anthropogenic NO','Biomass NO','Biofuel NO','Fertilizer NO','Lightning NO','Soil NO','Stratopsheric NO'};
isLNOx = logical([0 0 0 0 0 1 0 0]);
isSurfaceNOx = logical([0 1 1 1 1 0 1 0]);

if numel(isLNOx) ~= numel(categories)
    error(E.numelMismatch('isLNOx','categories'));
elseif numel(isSurfaceNOx) ~= numel(categories)
    error(E.numelMismatch('isSurfaceNOx','categories'));
end

struct_fields = {prod_struct(:).fullCat};
fields_found = ismember(categories,struct_fields);

if ~all(fields_found)
    % Create a message that lists the fields not found
    spec = [repmat('%s, ',1,sum(fields_found))-1,'%s'];
    msg_spec = sprintf('The fields %s were not found in the input structure',spec);
    msg = sprintf(msg_spec,categories{~fields_found});
    error(E.badinput(msg));
end

% Second make sure that all the date vectors are the same - they should be,
% but let's be careful.

for a=2:numel(prod_struct)
    if all(size(prod_struct(a).tVec == prod_struct(1).tVec))
        if ~all(prod_struct(a).tVec == prod_struct(1).tVec)
            error(E.badinput('The tVec for tracer #%d has different values than the first one',a));
        end
    else
        error(E.dimMismatch(sprintf('prod_struct(%d).tVec',a),'prod_struct(1).tVec'));
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% VARIABLE PREP %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%

% Go ahead and just set the date vector - we've already confirmed
% everything is kosher about it
gc_datevec = prod_struct(1).tVec;

% Find each production category - this will both ensure that any extra
% fields are not counted and allows us to easily break it down into
% lightning/surface/other

all_inds = zeros(size(categories));
for a=1:numel(categories)
   all_inds(a) = find(strcmp(categories{a},struct_fields));
end

% Double check that all of the categories were found - they should be
% (because we already checked this once) - but we don't want one of these
% indices to be unset
if any(isempty(all_inds))
    error(E.unknownError('One of the category indicies is unset - this should not have happened'));
end

% The matrix for fertilizer production should only be a 3D matrix - lon x
% lat x time. Since these are the dimensions we want each of our output
% matrices to have, we'll grab it.
gc_sz = size(prod_struct(all_inds(5)).dataBlock);
total_no = nan(gc_sz);
frac_lnox = nan(gc_sz);
frac_surface = nan(gc_sz);

NOAbs = make_empty_struct_from_cell(cat(2,categories,{'Total_NO'}));
[NOFrac, cat_fieldnames] = make_empty_struct_from_cell(categories);
Prod_Amt = struct('NO_molec_cm2_sec',NOAbs,'NO_frac',NOFrac);

%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% MAIN FUNCTION %%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%

% Go through each category specified - calculate the total emission
% (summing over altitude if multiple levels are present). Keep adding to
% the total production count

for a=1:numel(cat_fieldnames)
    if ndims(prod_struct(all_inds(a)).dataBlock) > 3
        sum_mat = squeeze(nansum2(prod_struct(all_inds(a)).dataBlock,3));
    else
        sum_mat = prod_struct(all_inds(a)).dataBlock;
    end
    Prod_Amt.NO_molec_cm2_sec.(cat_fieldnames{a}) = sum_mat;
    total_no = nansum2(cat(4,total_no,sum_mat),4);
end
Prod_Amt.NO_molec_cm2_sec.Total_NO = total_no;

% Now calculate the fractions.  We'll use isLNOx and isSurfaceNOx to
% determine when to add the fraction to the frac_lnox and frac_surface
% outputs
for a=1:numel(cat_fieldnames)
    frac_mat = Prod_Amt.NO_molec_cm2_sec.(cat_fieldnames{a}) ./ total_no;
    Prod_Amt.NO_frac.(cat_fieldnames{a}) = frac_mat;
    if isLNOx(a)
        frac_lnox = nansum2(cat(4,frac_lnox,frac_mat),4);
    elseif isSurfaceNOx(a)
        frac_surface = nansum2(cat(4,frac_surface,frac_mat),4);
    end
end

end

function chk  = check_struct_type(chk_struct)
    fns = fieldnames(chk_struct);
    expected_fns = {'dataBlock','dataUnit','fullName','fullCat','tVec','modelName','modelRes','dataScale','molMass','tEdge'};
    chk = all(ismember(expected_fns,fns));
end
