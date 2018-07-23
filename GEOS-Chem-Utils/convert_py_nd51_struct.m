function [ OutStruct ] = convert_py_nd51_struct( SatNO2, varargin )
%cleanup_py_nd51_struct Tidies up the structure output from read_gc_nd51.py
%   Python's bpch module makes it easier to access fields in the binary
%   punch (bpch) files that GEOS-Chem outputs, which in made it easier to
%   write something that could loop over all the files output when dealing
%   with satellite column output. While SciPy makes it possible to save
%   Matlab structures, it's not very easy to make structure arrays. This
%   function will take at least one structure (typically output from the
%   read_gc_nd51.py script, although structures output from
%   read_geos_output.m will also work).  It will ensure that all the
%   expected fields are present, as well as ensuring the expected field
%   order, and will concatenate all given structures into one.
%
%   Josh Laughner <joshlaugh5@gmail.com> 31 Mar 2015

E = JLLErrors;
narginchk(1,Inf)

%%%%% INPUT CHECKING %%%%%

% Convert SatNO2 to a structure if given as a cell array (which is how
% Python saves it)
if iscell(SatNO2)
    for a=1:numel(SatNO2)
        T(a) = SatNO2{a};
    end
    SatNO2 = T;
    clear('T');
elseif ~isstruct(SatNO2)
    E.badinput('SatNO2 must be a structure or a cell array')
end

if (~isempty(varargin) && ~all(iscellcontents(varargin,'isstruct')))
    E.badinput('All optional inputs must be structures')
end

expected_fields = {'dataBlock', 'dataUnit', 'fullName', 'fullCat', 'tVec', 'modelName', 'modelRes', 'dataScale', 'molMass', 'tEdge'};

flds = fieldnames(SatNO2);

if ~all(ismember(expected_fields,flds))
    E.badinput('All input structures are expected to have the fields %s',strjoin(expected_fields,', '));
end

for a=1:numel(varargin)
    flds = fieldnames(varargin{a});
    if ~all(ismember(expected_fields,flds))
        E.badinput('All input structures are expected to have the fields %s',strjoin(expected_fields,', '));
    end
end


%%%%% MAIN PROGRAM %%%%%

n=numel(SatNO2);

for a=1:numel(varargin)
    n=n+numel(varargin{a});
end

OutStruct = repmat(make_empty_struct_from_cell(expected_fields),1,n);
i=1;

[OutStruct,i] = append_struct(SatNO2,OutStruct,i);
for a=1:numel(varargin)
    [OutStruct,i] = append_struct(varargin{a},OutStruct,i);
end

end

function [OutStruct,i] = append_struct(struct_to_append, OutStruct, i)
    % Check that the structure to append has the same time vector as the
    % original structure. For now, we'll only worry about the day being the
    % same (hence the floor operator). 
    
    E = JLLErrors;
    if ~isempty(OutStruct(1).tVec) && any(floor(OutStruct(1).tVec) ~= floor(struct_to_append(1).tVec))
        E.callError('time_mismatch','The structure to append does not have the same time vector as the current structure');
    end

    expected_fields = {'dataBlock', 'dataUnit', 'fullName', 'fullCat', 'tVec', 'modelName', 'modelRes', 'dataScale', 'molMass', 'tEdge'};
    S = orderfields(struct_to_append,expected_fields); 
    j = i+numel(struct_to_append)-1;

    OutStruct(i:j) = S(:);
    i = j+1;
end

