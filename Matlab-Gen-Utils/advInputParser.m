classdef advInputParser < handle
    %ADVINPUTPARSER Improved version of the standard input parser
    %   This class recreates most of the functionality of the standard
    %   input parser, but with a few improvements:
    %
    %       1. Ability to add "flags", which are true/false parameters set
    %       to true simply by the presence of the proper string in the
    %       input parameters.
    %
    %       2. Order of parsing is required, flags, parameters, optional;
    %       this circumvents an issue in the regular input parser where an
    %       optional character argument will consume parameter names.
    
    properties
        AllowFlagRepeats;
        CaseSensitive;
        FunctionName = '';
        KeepUnmatched = false;
    end
    
    properties(SetAccess = private)
        Unmatched;
        Required;
        Optional;
        Parameters;
        Flags;
        Results;
        AdvResults;
    end
    
    methods
        
        function obj = advInputParser()
            
            obj.Flags = struct;
            obj.Parameters = struct;
            obj.Optional = struct;
            obj.Required = struct;
            obj.AllowFlagRepeats = true;
        end
        
        function val = get.AdvResults(obj)
            % Provided for compatibility with the previous version of
            % advInputParser.
            val = obj.Results;
        end
        
        function addRequired(obj, name, varargin)
            % addRequired(name)
            % addRequired(name, validation_fxn)
            % addRequired(name, validation_fxn, error_message)
            %
            %   Adds a required argument to the parser. 
            
            obj.check_duplicate_fields(name);
            % This method requires a default value, but the require fields
            % don't use it, so just provide a dummy value
            obj.Required.(name) = advInputParser.make_substruct(name, 'required param - default not used', varargin{:});
        end
        
        function addOptional(obj, name, default, varargin)
            % addOptional(name, default)
            % addOptional(name, default, validation_fxn)
            % addOptional(name, default, validation_fxn, error_message)
            %
            %   Adds an optional positional argument to the parser. 
            obj.check_duplicate_fields(name);
            obj.Optional.(name) = advInputParser.make_substruct(name, default, varargin{:});
        end
        
        function addParameter(obj, name, default, varargin)
            % addParameter(name, default)
            % addParameter(name, default, validation_fxn)
            % addParameter(name, default, validation_fxn, error_message)
            %
            %   Adds an optional parameter argument to the parser. 
            obj.check_duplicate_fields(name);
            obj.Parameters.(name) = advInputParser.make_substruct(name, default, varargin{:});
        end
        
        function addFlag(obj, name)
            % addFlag(name)
            %
            %   Adds an optional flag argument to the parser. 
            
            obj.check_duplicate_fields(name);
            obj.Flags.(name) = advInputParser.make_substruct(name, false);
        end
        
        function parse(obj, varargin)
            % First check for the right number of required arguments, and
            % remove them. Then look for flags, parameters, and optional
            % arguments, in that order.
            
            % Of the properties, we will need to implement CaseSensitive
            % and StructExpand ourselves. We could implement partial
            % matching, but the documentation (for R2017b) indicates that
            % that is only for parameters anyway.
            
            
            obj.create_default_results();
            
            inputs = varargin;
            inputs = obj.parse_required(inputs);
            inputs = obj.parse_flags(inputs);
            inputs = obj.parse_parameters(inputs);
            inputs = obj.parse_optional(inputs);
            
            if ~isempty(inputs) && ~obj.KeepUnmatched
                obj.printError('%d arguments could not be parsed', numel(inputs));
            else
                obj.Unmatched = inputs;
            end
        end
    end
    
    methods(Access = private)
        function printError(obj, message, varargin)
            % I'm not sure what to call on inputParser (or even if it can
            % be called from a subclass) to do this, so we'll write our
            % own. This method will print the requested error message,
            % optionally prefaced by the "Error in {}" where {} is the
            % property FunctionName
            if ~isempty(obj.FunctionName)
                msg_string = sprintf(message, varargin{:});
                formatted_msg = sprintf('Error using <a href="matlab: matlab.internal.language.introspective.errorDocCallback(''%1$s'')">%1$s</a>:\n%2$s', obj.FunctionName, msg_string);
            else
                formatted_msg = sprintf(message, varargin{:});
            end
            
            err_struct = struct('message', formatted_msg,...
                'identifier', 'MATLAB:InputParser:ArgumentFailedValidation',...
                'stack', dbstack(3));
            % We omit the first levels in the stack because those will
            % point to the input parser, and we don't want that - we want
            % the calling function.
            error(err_struct);
        end
        
        function parser_structs = all_fields(obj)
            parser_structs = advInputParser.combine_structs(obj.Required, obj.Optional);
            parser_structs = advInputParser.combine_structs(parser_structs, obj.Parameters);
            parser_structs = advInputParser.combine_structs(parser_structs, obj.Flags);
        end
        
        function check_duplicate_fields(obj, new_name)
            if ~ischar(new_name)
                error('advInputParser:add:badInput','NAME must be a character array');
            end
            
            compare_fxn = obj.get_comparison_fxn();
            
            parser_fields = obj.all_fields();
            fns = fieldnames(parser_fields);
            if any(compare_fxn(fns, new_name))
                error('advInputParser:add:badInput', 'The name "%s" is already used in the parser', new_name);
            end
        end
        
        function compare_fxn = get_comparison_fxn(obj)
            if obj.CaseSensitive
                compare_fxn = @strcmp;
            else
                compare_fxn = @strcmpi;
            end
        end
        
        function create_default_results(obj)
            obj.Results = struct();
            parser_fields = obj.all_fields();
            fns = fieldnames(parser_fields);
            for i_fn = 1:numel(fns)
                p_name = fns{i_fn};
                obj.Results.(p_name) = parser_fields.(p_name).default;
            end
        end
        
        function inputs = parse_required(obj, inputs)
            req_fns = fieldnames(obj.Required);
            n_req = numel(req_fns);
            if numel(inputs) < n_req
                obj.printError('At least %d arguments required (only %d given)', n_req, numel(inputs));
            end
            
            for i_fn = 1:n_req
                obj.Results.(req_fns{i_fn}) = inputs{i_fn};
            end
            
            inputs = inputs(n_req+1:end);
        end
        
        function inputs = parse_flags(obj, inputs)
            compare_fxn = obj.get_comparison_fxn();
            flags = fieldnames(obj.Flags);
            inputs_to_remove = false(size(inputs));
            for i_in = 1:numel(inputs)
                for i_flags = 1:numel(flags)
                    % If the input doesn't match the current flag, no need
                    % to do anything else
                    if ~compare_fxn(inputs{i_in}, flags{i_flags})
                        continue
                    end
                    
                    % If it does match, then we need to check if it was
                    % already matched. If AllowFlagRepeats is true
                    % (default) we just mark it to be removed and move on.
                    % Otherwise, it is an error.
                    inputs_to_remove(i_in) = true;
                    if obj.Results.(flags{i_flags}) && ~obj.AllowFlagRepeats
                        obj.printError('Input parsing error: The flag "%s" was given twice', flags{i_flags});
                    end
                    obj.Results.(flags{i_flags}) = true;
                end
            end
            
            inputs(inputs_to_remove) = [];
        end
        
        function inputs = parse_parameters(obj, inputs)
            compare_fxn = obj.get_comparison_fxn();
            params = fieldnames(obj.Parameters);
            for i_param = 1:numel(params)
                % For each parameter, find the last instance of it in the
                % inputs and get the value. Check that it passes the
                % validation function before adding it to results.
                
                xx = find(compare_fxn(inputs, params{i_param}), 1, 'last');
                if isempty(xx)
                    continue
                elseif xx == numel(inputs)
                    obj.printError('No input value given after the parameter name "%s"', params{i_param});
                end
                
                val = inputs{xx+1};
                val_fxn = obj.Parameters.(params{i_param}).validation_fxn;
                chk = val_fxn(val);
                if ~isscalar(chk) || ~islogical(chk)
                    obj.printError('The value for parameter "%s" return a non-scalar or non-logical value from the validation function: %s', params{i_param}, func2str(val_fxn));
                elseif ~chk
                    obj.printError(obj.Parameters.(params{i_param}).error_msg);
                end
                
                obj.Results.(params{i_param}) = val;
                inputs(xx:xx+1) = [];
            end
        end
        
        function inputs = parse_optional(obj, inputs)
            opts = fieldnames(obj.Optional);
            for i_opt = 1:numel(opts)
                % Get as many optional arguments as we can before running
                % out of inputs. If we used all the inputs, then we can
                % collapse inputs to an empty cell array.
                if i_opt > numel(inputs)
                    inputs = {};
                    return
                end
                
                val = inputs{i_opt};
                val_fxn = obj.Optional.(opts{i_opt}).validation_fxn;
                chk = val_fxn(val);
                if ~isscalar(chk) || ~islogical(chk)
                    obj.printError('The value for optional input #%d return a non-scalar or non-logical value from the validation function: %s', i_opt, func2str(val_fxn));
                elseif ~chk
                    obj.printError(obj.Optional.(opts{i_opt}).error_msg);
                end
                
                obj.Results.(opts{i_opt}) = val;
            end
            
            % If we haven't used all the inputs, then we should return
            % what's left for the unmatched structure.
            inputs = inputs(numel(opts)+1:end);
        end
        
        
    end
    
    methods(Access = private, Static)
        function S = make_substruct(name, default, varargin)
            if numel(varargin) < 1
                val_fxn = @(x) true;
            else
                val_fxn = varargin{1};
            end
            
            if numel(varargin) < 2
                msg = sprintf('The value for "%s" failed validation using the function: %s', name, func2str(val_fxn));
            else
                msg = varargin{2};
            end

            % Because of how struct() interprets cell arrays, any default values that are cell arrays
            % must be wrapped in an outer cell array to prevent the structure from being non-scalar
            if iscell(default)
                default = {default};
            end

            S = struct('name', name, 'default', default, 'validation_fxn', val_fxn, 'error_msg', msg);
        end
        
        
        function S = make_struct_default_val(fields, val)
            s_cell = cell(numel(fields)*2,1);
            s_cell(1:2:end) = fields;
            s_cell(2:2:end) = {val};
            S = struct(s_cell{:});
        end
        
        function S_combo = combine_structs(S1, S2)
            fields1 = fieldnames(S1);
            fields2 = fieldnames(S2);
            if any(ismember(fields1, fields2))
                error('advInputParser:combine_structs:dup_field', 'One or more fields are common in S1 and S2')
            end
            
            S_combo = advInputParser.make_struct_default_val(cat(1, fields1, fields2), []);
            for i_fields = 1:numel(fields1)
                S_combo.(fields1{i_fields}) = S1.(fields1{i_fields});
            end
            for i_fields = 1:numel(fields2)
                S_combo.(fields2{i_fields}) = S2.(fields2{i_fields});
            end
        end
        
    end
    
    
    
end

