function [ new_struct, safe_names ] = make_empty_struct_from_cell( names, default_val )
%make_empty_struct_from_cell Makes an empty structure with given field names 
%   Takes an input cell array of field names as strings and generates a
%   structure with those field names referring to empty matrices.  The
%   second input is the default value you wish each field to have.
%
%   This function will return the structure (of course) but will also
%   return the field names as a second structure - this is useful if it
%   changes some names to make them structure safe (and saves you a call to
%   fieldnames).

E = JLLErrors;

%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% INPUT CHECKING %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%

narginchk(1,2);
if ~iscell(names) || ~all(iscellcontents(names,'ischar'))
    error(E.badinput('"names" must be a cell array of strings'));
end

if nargin < 2
    default_val = [];
end

%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% MAIN FUNCTION %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%

% Make the names safe as field names, pass a warning if they are changed

safe_names = regexprep(names,'\W','_');
for a=1:numel(safe_names)
    if ~isempty(regexp(safe_names{a}(1),'[0-9]', 'once'))
        safe_names{a} = strcat('FN_',safe_names{a});
    elseif ~isempty(regexp(safe_names{a}(1),'_', 'once'))
        safe_names{a} = strcat('FN',safe_names{a});
    end
end
if any(~ismember(safe_names,names))
    warning('Some field names were made safe for field names and will not match the original name');
end

tmp_cell = cell(1,2*numel(names));
for a=1:2*numel(names)
    % Make the odd entries the field names and the even entries empty
    % matrices - this way tmp_cell{:} will generate the proper inputs to
    % the "struct" function (note that cell(m,n) produces a cell array with
    % empty matrices by default, so we only need to set the odd entries,
    % unless the user overrides the default value of the fields).
    if mod(a,2)==1
        ind = (a+1)/2;
        tmp_cell{a} = safe_names{ind};
    elseif mod(a,2) == 0
        tmp_cell{a} = default_val;
    end
end

new_struct = struct(tmp_cell{:});

end

