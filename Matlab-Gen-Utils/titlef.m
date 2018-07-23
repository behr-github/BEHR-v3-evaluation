function [  ] = titlef( varargin )
%titlef Uses formatted strings for title
%   To use special characters, such as a line break (\n) in a title, the
%   string must be returned from an sprintf function.  This will do that
%   automatically, to save typing.

if numel(varargin) == 1;
    title(sprintf(varargin{1}));
else
    title(sprintf(varargin{1}),varargin{2:end});
end


end

