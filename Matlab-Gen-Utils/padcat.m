function [ M ] = padcat( dim, varargin )
%M = PADCAT( DIM, ... ) Concatenates different sizes matrices along DIM
%   Similar to the built-in CAT, this function takes any number of input
%   arrays and concatenates them along the dimension DIM. However, this
%   function will allow the concatenation of differently sized arrays by
%   padding them with NaNs in the non-concatenated dimensions to match the
%   size of the largest array.
%
%   M = PADCAT( DIM, ... , 'padval', PADVAL ) allows you to change the
%   value used for padding to PADVAL.

%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% INPUT CHECKING %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%

E = JLLErrors;

if numel(varargin) < 2
    % Yes only one array input isn't really being concatenated, but let's
    % not error in a case we don't need to.
    E.badinput('PADCAT requires at least two inputs')
end

% Parse the input arrays, reading and removing the PADVAL parameter
xx = strcmpi('padval',varargin);
if any(xx)
    i = find(xx);
    try
        padval = varargin{i+1};
    catch err
        if strcmpi(err.identifier, 'MATLAB:badsubscript')
            E.badinput('The parameter name ''padval'' must be followed by a padding value');
        else
            rethrow(err)
        end
    end
    varargin(i:i+1) = [];
else
    padval = NaN;
end

ndimmax = 0;
for i=1:numel(varargin)
    ndimmax = max(ndimmax, ndims(varargin{i}));
end

if ~isnumeric(dim) || ~isscalar(dim) || dim > ndimmax
    E.badinput('dim must be a valid dimension to concatenate along')
end

%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% MAIN FUNCTION %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%

% Find the find the maximum dimensions of the arrays passed in
sz = zeros(numel(varargin), ndimmax);
for i=1:numel(varargin)
    sz_i = size(varargin{i});
    sz(i, 1:numel(sz_i)) = sz_i;
end
szmax = max(sz,[],1);

% For each array input, pad it with the padding value along the
% non-concatenation dimensions and concatenate

M = [];
for i=1:numel(varargin)
    this_array = varargin{i};
    padvec = szmax - sz(i,:);
    padvec(dim) = 0;
    this_array = padarray(this_array,padvec,padval,'post');
    M = cat(dim, M, this_array);
end

end

