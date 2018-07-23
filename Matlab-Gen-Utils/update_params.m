function parameters = update_params(parameters,varargin)
%UPDATE_PARAMETERS Update parameter values in a cell array of name-value pairs
%   PARAMETERS = UPDATE_PARAMS(PARAMETERS, NAME1, VALUE1, NAME2, VALUE2,
%   ...) will look for a char array or string matching each NAME in
%   PARAMETERS and replace the following element with the respective VALUE.
%   This is intended for updating "varargin" input arrays that contain
%   name-value pairs that need to be modified before passing to subordinate
%   functions.
%
%   PARAMETERS = UPDATE_PARAMS('remove', PARAMETERS, NAME1, NAME2, ...)
%   will remove the parameters called NAME1, NAME2, etc.
%
%   PARAMETERS = UPDATE_PARAMS('remove', PARAMETERS, S) will take the field
%   names of the structure S as the parameters to remove. This is useful if
%   you have a structure returned by an input parser.
%
%   PARAMETERS = UPDATE_PARAMS('missing', PARAMETERS, NAME1, VALUE1, NAME2,
%   VALUE2, ...) will add the parameters NAME1, NAME2, etc. if they are
%   missing but will not replace them if they already exist.

E = JLLErrors;

add_parameters = true;
only_add = false;
if ischar(parameters)
    if strcmpi(parameters, 'remove')
        add_parameters = false;
        parameters = varargin{1};
        varargin = varargin(2:end);
    elseif strcmpi(parameters, 'missing')
        only_add = true;
        parameters = varargin{1};
        varargin = varargin(2:end);
    else
        E.badinput('The first input must be the char array "remove" or a cell array');
    end
end

if ~iscell(parameters)
    E.badinput('PARAMETERS must be a cell array')
end



% We're not going to check that parameters has an even number of elements
% or that it alternates name and value because it might include some
% optional but unnamed input values.

if add_parameters
    if numel(varargin) < 2
        E.badinput('Must provide at least one parameter name-value pair')
    elseif mod(numel(varargin),2) ~= 0
        E.badinput('There must be an even number of arguments after the original parameters cell array (alternating parameter name and value)');
    end
else
    if numel(varargin) < 1
        E.badinput('Must provide at least one parameter name to remove');
    end
    
    if isstruct(varargin{1})
        if numel(varargin) == 1
            varargin = fieldnames(varargin{1});
        else
            E.badinput('If providing a structure whose field names are the parameter names to remove, you can only give one such structure and no other parameter names'); 
        end
    end
end

if add_parameters
    incr = 2;
else
    incr = 1;
end

for i_param = 1:incr:numel(varargin)
    this_param = varargin{i_param};
    
    if ~ischar(this_param)
        if add_parameters
            E.badinput('The first, third, fifth, etc. inputs after the parameters cell must be char arrays or strings (these are the parameter names)')
        else
            E.badinput('In "remove" mode, all additional arguments after the parameter array must be char arrays')
        end
    end
    
    parameter_missing = false;
    param_idx = find(strcmp(parameters, this_param));
    if isempty(param_idx)
        %E.callError('unknown_parameter', 'The parameter named "%s" does not exist in the parameters cell array', this_param);
        parameter_missing = true;
    elseif isscalar(param_idx) && mod(param_idx, 2) ~= 1
        %E.callError('parameter_name_matched_value', 'The parameter name "%s" only appears as a parameter value', this_param);
        parameter_missing = true;
    elseif ~isscalar(param_idx)
        % Try removing the even parameter indices; it's possible that the
        % parameter name is supposed to show up as a value
        param_idx(mod(param_idx,2) ~= 1) = [];
        if isempty(param_idx)
            %E.callError('parameter_name_matched_value', 'The parameter name "%s" only appears as a parameter value', this_param);
            parameter_missing = true;
        elseif ~isscalar(param_idx)
            E.callError('multiple_matched_parameters', 'The parameter name "%s" exists multiple times in the parameters cell array', this_param);
        end
    end
    
    % All that error checking, now all we need to do is change the
    % parameter value that follows the name or remove the parameter and its
    % value.
    if add_parameters
        if parameter_missing
            parameters(end+1:end+2) = varargin(i_param:i_param+1);
        elseif ~only_add
            % If the "missing" flag was given, then we should only add the
            % parameter if it is missing, not overwrite an existing
            % parameter value.
            parameters{param_idx+1} = varargin{i_param+1};
        end
    else
        if ~parameter_missing
            parameters(param_idx:param_idx+1) = [];
        end
    end
end

end

