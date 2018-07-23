function [ xx ] = inranges( Data, ranges )
%INRANGE Returns true for any value found in the ranges defined.
%   In several of my functions, I define or return a list of ranges as an
%   n-by-2 matrix where each row is a pair of bounding values. This
%   function takes such a matrix and applies it to numerical data,
%   returning a logical matrix true for any value that can be found in any
%   of the ranges given.
%
%   Josh Laughner <joshlaugh5@gmail.com> 19 June 2015

E = JLLErrors;

%%%%%%%%%%%%%%%%%%%%%%%%
%%%% INPUT PARSING %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%

if ~isnumeric(Data)
    E.badinput('Data must be numeric');
end

if size(ranges,2) ~= 2
    E.badinput('ranges must be n-by-2')
elseif ~isnumeric(ranges)
    E.badinput('ranges must be numeric');
elseif any(ranges(:,1) >= ranges(:,2))
    E.badinput('A range is defined incorrectly: the first column must contain the lower bound')
elseif any(isnan(ranges(:)))
    E.badinput('NaNs are not allowed as a range boundary')
end

%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% MAIN FUNCTION %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%

xx = false(size(Data));
for a=1:size(ranges,1)
    xx = xx | (Data >= ranges(a,1) & Data < ranges(a,2));
end

end

