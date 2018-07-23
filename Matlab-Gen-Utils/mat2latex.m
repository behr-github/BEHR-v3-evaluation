function varargout = mat2latex(M, format_spec, uncertainty_dim, asymmetry)
% MAT2LATEX( M ) Prints a matrix, M in Latex-table format
%   Takes a numeric matrix or a cell array  and prints out to the command window that matrix
%   formatted such that it can be copied into a Latex table. 
%
%   MAT2LATEX( M, FORMAT_SPEC ) Optional second argument is a format spec
%   string understood by fprintf to control how the numbers will be
%   formatted. Defaults to %g.
%
%   MAT2LATEX( M, UNCERTAINTY_DIM ) Allows you to include uncertainty
%   values that will be printed as v \pm u (x 10^e). Uncertainties must
%   alternate with values along the dimension given. For example, if this
%   is 1, then the first row of M should contain values, the second row the
%   uncertainties in the first row, the third row the next set of values,
%   the fourth row the uncertainties for the third row, and so on.
%   
%   MAT2LATEX( M, FORMAT_SPEC, UNCERTAINTY_DIM ) combines the previous two
%   syntaxes.
%
%   MAT2LATEX( M, 'uncertainty', UNCERTAINTY_DIM ) will format any numbers
%   such that the last digit in the value is the first digit in the
%   uncertainty.
%
%   MAT2LATEX( M, 'u', UNCERTAINTY_DIM ) is a shorthand for the last
%   syntax.
%
%   MAT2LATEX( M, 'u10', UNCERTAINTY_DIM ) will alter the output slightly
%   so that any value that would be printed in exponential with '10^{1}' is
%   printed in regular notation instead (i.e. 9. \pm 4 \times 10^{1} would
%   be given as just 90 \pm 40). This does not alter values with negative
%   exponents. Note that you can use any power of 10 instead; 'u100',
%   'u1000', 'u1000000' are all valid as well.
%
%   MAT2LATEX( M, ___, UNCERTAINTY_DIM, 'asym' ) will allow there to be
%   asymmetric uncertainty along the specified dimension. The matrix, M,
%   must rotate [value, lower uncertainty, upper uncertainty] along the
%   UNCERTAINTY_DIM. You may omit a format string to use the default. If
%   one of the uncertainty values is a NaN, then that value will use
%   symmetrical uncertainty formatting. If both are NaNs, an error is
%   thrown.

default_format_spec = '%g';
if nargin < 2
    format_spec = default_format_spec;
    uncertainty_dim = 0; % do not include uncertainty. Deafult asym_bool set below
elseif nargin == 2
    % If only two arguments, determine if the format spec or uncertainty
    % dimension was passed
    if isnumeric(format_spec)
        uncertainty_dim = format_spec;
        format_spec = '%g';
    else
        uncertainty_dim = 0;
    end
elseif nargin == 3
    % Only two possibilities: either a format string and uncertainty dim
    % have been passed, or an uncertainty dim and the 'asym' string
    if isnumeric(format_spec)
        asymmetry = uncertainty_dim;
        uncertainty_dim = format_spec;
        format_spec = default_format_spec;
    end
end

asym_bool = false;
if exist('asymmetry','var') && ischar(asymmetry) && strcmpi(asymmetry,'asym')
    asym_bool = true;
end
if asym_bool && uncertainty_dim == 0
    E.badinput('A non-zero uncertainty dimension must be specified to use asymmetric uncertainty')
end


E = JLLErrors;
if ~iscell(M) && ~isnumeric(M)
    E.badinput('Expecting a numeric matrix or a cell array')
elseif ~ismatrix(M)
    E.badinput('Higher dimension arrays don''t make sense as a table')
end

if ~ischar(format_spec)
    E.badinput('FORMAT_SPEC must be a string')
end
if ~isnumeric(uncertainty_dim) || uncertainty_dim < 0 || uncertainty_dim > 2
    E.badinput('UNCERTAINTY_DIM must be 0 (no uncertainties given), 1, or 2')
elseif uncertainty_dim > 0 && mod(size(M,uncertainty_dim),2) ~= 0 && ~asym_bool
    E.badinput('M must have an even number of entries along the UNCERTAINTY_DIM when using symmetrical uncertainties')
elseif uncertainty_dim > 0 && mod(size(M,uncertainty_dim),3) ~= 0 && asym_bool
    E.badinput('M must have a multiple of 3 entries along the UNCERTAINTY_DIM when using asymmetrical uncertainties')
end

if ~iscell(M)
    M = num2cell(M);
end

warn_neg_uncert = false;

% The regular print string
fstr1 = sprintf('%s & ', format_spec);
fstr1b = '%s & ';
% The end-of-line print string
fstr2 = sprintf(tex_in_printf('%s \ \n',3), format_spec);
fstr2b = tex_in_printf('%s \ \n',3);

sout = '';
print_to_screen = nargout == 0;

sz = size(M);
uncert_step = asym_bool + 2; % a little trickery to give 2 if using asymmetrical uncertainty, 1 if using symmetrical uncertainty, and 0 if using no uncertainty
if uncertainty_dim == 1
    astep = uncert_step;
else
    astep = 1;
end

if uncertainty_dim == 2
    bstep = uncert_step;
else
    bstep = 1;
end
for a=1:astep:sz(1)
    for b=1:bstep:sz(2)
        if isnumeric(M{a,b})
            if isnan(M{a,b})
                this_fstr = format_nan;
            elseif strcmpi(format_spec,'uncertainty') || (~isempty(regexp(format_spec, 'u\d*', 'once')) && isempty(strfind(format_spec, '%')))
                maxplace = str2double(regexp(format_spec,'\d*','match','once'));
                
                if asym_bool
                    this_fstr = format_asym_uncert(maxplace);
                else
                    this_fstr = format_uncert(maxplace);
                end
            else
                this_fstr = format_normal;
            end
            
            lprintf('$%s',this_fstr);
        else
            if b < (sz(2) - bstep+1)
                lprintf(fstr1b, M{a,b});
            else
                lprintf(fstr2b, M{a,b});
            end
        end
    end
end

lprintf('\n')

if warn_neg_uncert
    warning('Negative values of uncertainty have been replaced with positive ones');
end

if nargout > 0
    varargout{1} = sout;
end

% Nested functions for output
    function lprintf(s, varargin)
        if print_to_screen
            fprintf(s, varargin{:});
        else
            stemp = sprintf(s, varargin{:});
            sout = [sout, stemp];
        end
    end

% Nested functions to parse numbers
    function fstr = format_asym_uncert(maxplace)
        % This is similar to format_uncert, except that we need to
        % accomodate two different uncertainties. This means that the
        % string itself will be different (the upper uncertainty as
        % superscript and lower as subscript) but we also need to choose
        % the place to round to from the smaller of the two uncertainties.
        v = M{a,b};
        if uncertainty_dim == 1
            u = [M{a+1,b}, M{a+2,b}]; % lower then upper uncertainty
        elseif uncertainty_dim == 2
            u = [M{a,b+1}, M{a,b+2}];
        else
            E.badinput('Cannot format the values by rounding to the first place of the uncertainty if no uncertainty given (uncertainty_dim must be >0).');
        end
        
        % If there is a NaN as one of the uncertainties, then assume we
        % want a symmetrical uncertainty mixed in with the asymmetrical
        % ones.
        u_nans = isnan(u);
        if sum(u_nans) == 1
            u = u(~u_nans);
        elseif sum(u_nans) > 1
            E.badinput('Both uncertainties cannot be NaNs');
        end
        
        fstr = format_general_uncert(v, u, maxplace);
    end

    function fstr = format_uncert(maxplace)
        % Get the value and uncertainty following the usual rules, but
        % round them to the first non-zero place in the uncertainty first.
        v = M{a,b};
        if uncertainty_dim == 1
            u = M{a+1,b};
        elseif uncertainty_dim == 2
            u = M{a,b+1};
        else
            E.badinput('Cannot format the values by rounding to the first place of the uncertainty if no uncertainty given (uncertainty_dim must be >0).');
        end
        
        fstr = format_general_uncert(v, u, maxplace);
    end

    function fstr = format_general_uncert(v, u, maxplace)
        if any(u<0)
            warn_neg_uncert = true;
        end
        u = abs(u);
        place = min(10 .^ floor(log10(u)));
        v = round(v / place)*place;
        u = round(u ./ place) .* place;
        % The effect is that the larger uncertainty may have extra places,
        % but the smaller one will always have one.
        
        % Calculate the number of significant figures to be left in v
        % Testing if log10 is infinite rather than just u == 0 or v == 0
        % because it is if log10 = -Inf that directly causes the problem,
        % so I want to avoid any little floating point errors
        
        if isinf(log10(v)) && all(isinf(log10(u)))
            if asym_bool
                fstr = '0^{+0}_{-0}$';
            else
                fstr = '0 \pm 0$';
            end
        elseif ~isinf(log10(v))
            % Once again we will round off to the smaller of the two
            % uncertainties, if using asymmetric uncertainty.
            
            nsig = floor(log10(v)) - min(floor(log10(u))) + 1;
            fstr = sprintf('%#.*g',nsig,v);
        elseif any(~isinf(log10(u)))
            % If the value has become 0, then we will need to use u to set
            % the format, as long as it isn't 0! We cheat by formatting
            % u to have one significant digit, then replacing that
            % digit with a zero to become v, e.g. if u is 0.1, this
            % produces 0.0 for v (not u) which has the right number of
            % places. Again we use the smaller of u,
            ustr = sprintf('%.1g',min(u(~isinf(u))));
            fstr = regexprep(ustr,'\d','0');
        else
            E.unknownError('How did we get here?')
        end
        
        if b < (sz(2) - bstep+1)
            fstr=sprintf(fstr1b, fstr);
        else
            fstr=sprintf(fstr2b, fstr);
        end
        fstr = format_exponent(fstr);
        if ~isnan(maxplace)
            fstr = remove_exponent(fstr, maxplace);
        end
        fstr = insert_uncertainty(u, fstr);
        fstr = strrep(fstr,'\\\\','\\');
    end

    function fstr = format_normal
        if b < (sz(2) - bstep+1)
            fstr=sprintf(fstr1, M{a,b});
        else
            fstr=sprintf(fstr2, M{a,b});
        end
        
        fstr = format_exponent(fstr);
        
        % Add in uncertainty, if desired
        if uncertainty_dim > 0
            if uncertainty_dim == 1
                if ~asym_bool
                    u = M{a+1,b};
                else
                    u = [M{a+1,b}, M{a+2,b}];
                end
            elseif uncertainty_dim == 2
                if ~asym_bool
                    u = M{a,b+1};
                else
                    u = [M{a,b+1}, M{a,b+2}];
                end
            end
            if any(u < 0)
                warn_neg_uncert = true;
                u = abs(u);
            end
            
            % If there is a NaN as one of the uncertainties, then assume we
            % want a symmetrical uncertainty mixed in with the asymmetrical
            % ones.
            u_nans = isnan(u);
            if sum(u_nans) == 1
                u = u(~u_nans);
            elseif sum(u_nans) > 1
                E.badinput('Both uncertainties cannot be NaNs');
            end
            
            fstr = insert_uncertainty(u, fstr);
        end
    end

    function fstr = format_nan
        if b < (sz(2) - bstep+1)
            fstr = sprintf(fstr1b, '-$');
        else
            fstr = sprintf(fstr2b, '-$');
            fstr = strrep(fstr, '\\\\','\\');
        end
    end

end

function fstr = format_exponent(fstr)
% Replace e notation with x10^, remove + sign and leading 0s
fstr = regexprep(fstr,'e',' \\times 10^{');
fstr = regexprep(fstr,'(?<={-?)+?0*(?=\d+)','');
%fstr = regexprep(fstr,'(?<={)\+?0*(?=\d+)','');
fstr = regexprep(fstr,'(?<={-?\d+) (?=&|\\)','}$ ');
if ~isempty(strfind(fstr,'{')) && isempty(strfind(fstr,'}'))
    fstr = strcat(fstr,'}$'); % handles the final line
end
% Numbers without an exponent will not have the ending $ added
% yet
if isempty(strfind(fstr,'$'))
    fstr = regexprep(fstr,' (?=&|\\)','$ ');
end
end

function fstr = remove_exponent(fstr, maxplace)
% This will convert numbers from exponential notation to normal up to the
% power specified, so remove_exponent(fstr, 10) will make any powers of
% 10^1 or 10^-1 into just regular notation, (fstr, 100) up to 10^2, 10^-2,
% and so on.

% First, see if there is an exponent, if not we can return right now.
s = regexp(fstr, '10\^{','once');
if isempty(s)
    return
end

% Otherwise find the exponent. If it is outside our range, we can return.
[es, ee] = regexp(fstr,'(?<=10\^{)-?\d*(?=})');

e = str2double(fstr(es:ee));
if e < 0 || e > log10(maxplace)
    return
end

% Get the value, scale it as necessary, then convert back into string
origval = regexp(fstr,'\d*\.?\d*(?=\o{40})','once','match');
val = str2double(origval)*10^e;
valstr = num2str(val);

% Remove the x 10^{} part
fstr = regexprep(fstr, '\o{40}\\times.*}', '');
% and replace the original value
fstr = regexprep(fstr, origval, valstr);

end

function fstr = insert_uncertainty(u, fstr)
% convert u to the same power of 10 as the value. first
% find the exponent in the string and convert it to a
% number.
[es, ee] = regexp(fstr,'(?<=10\^{)-?\d*(?=})');
if isempty(es)
    val_exp = 0;
else
    val_exp = str2double(fstr(es:ee));
end
u = u .* 10^(-val_exp);
% figure out how many figures after the decimal point there
% are in the value
[ds,de] = regexp(fstr,'\.\d*');
if ~isempty(ds)
    ndec = de - ds;
else
    ndec = 0;
    de = regexp(fstr, '(\$|\o{40}\\times)')-1;
end

if isscalar(u)
    % This handles symmetrical uncertainty pretty easily. Convert u to have
    % the right number of decimal places, then insert
    ustr = sprintf('%.*f', ndec, u);
    % insert the uncertainty
    fstr = sprintf('%s \\pm %s%s',fstr(1:de),ustr,fstr(de+1:end));
else
    % If using asymmetrical uncertainty, it's very similar because we've
    % deliberately rounded u so that it's last digit is the first digit of
    % the smaller uncertainty. So if we turn both uncertainties into
    % numbers with the same amount of post-decimal digits, the larger one
    % may have several non-zero digits, but the smaller one should only
    % have 1.
    ustr_l = sprintf('%.*f',ndec,u(1));
    ustr_u = sprintf('%.*f',ndec,u(2));
    fstr = sprintf('%s_{-%s}^{+%s} %s',fstr(1:de),ustr_l,ustr_u,fstr(de+1:end));
end
end

function [str, dec_ind] = round_str(str, index)
% Rounds the number in the string to the given index
dec_ind = regexp(str,'\.');
if isempty(dec_ind)
    dec_ind = length(str)+1;
elseif numel(dec_ind) > 1
    error('Multiple decimal points found');
end

% Handle rounding. Round numerically to the requested precision, then we
% will reformat the string afterwards.
num = str2double(str);
if index < dec_ind
    place = 10^(index - dec_in + 1);
else
    place = 10^(index - dec_in);
end
num = round(num * place)/place;

% Now either make the following numbers 0s or remove them, as appropriate
if index < dec_ind
    str(index+1:dec_ind-1)='0';
    str(dec_ind:end)=[];
else
    str(index+1:end)=[];
end
end
