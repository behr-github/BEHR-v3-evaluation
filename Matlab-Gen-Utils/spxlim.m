function [ ] = spxlim( xlimvec, varargin )
%SPXLIM( XLIM ) Set the x limits of all subplots in current fig to XLIM
%   Iterates through all subplot child axes of the current figure and sets
%   their x limits to the two vector element XLIM. Will skip any with the
%   tag 'suptitle' to avoid breaking them.
%
%   SPXLIM( FH, XLIM ) will instead set the x limits for the child subplots
%   of the figure specified by the handle FH.

E = JLLErrors;
%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% INPUT PARSING %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%

if isa(xlimvec, 'figure_handle')
    fh = xlimvec;
    if numel(varargin) < 1
        E.badinput('Must pass a two element vector as well as the figure handle')
    else
        xlimvec = varargin{1};
    end
else
    fh = gcf;
end

    

%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% MAIN FUNCTION %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%

ch = get(fh,'children');
for i=1:numel(ch)
    if strcmpi(ch(i).Tag,'suptitle')
        continue
    end
    xlim(ch(i), xlimvec);
end
end

