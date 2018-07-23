function [ GC ] = subset_gc_timeper( GC, start_date, end_date )
% GC = SUBSET_GC_TIMEPER( GC, START_DATE, END_DATE )
%   Cuts down the structure GC to the time period specified by START_DATE
%   and END_DATE, both of which must be in the time period defined by
%   GC.tVec. GC must be a structure output from either READ_GEOS_OUTPUT or
%   CONVERT_PY_ND51_STRUCT.

E = JLLErrors;

if ~isstruct(GC) || ~any(isfield(GC,{'dataBlock','tVec'}))
    E.badinput('GC must be a structure output from either READ_GEOS_OUTPUT or CONVERT_PY_ND51_STRUCT.')
end

try
    start_datenum = datenum(start_date);
    end_datenum = datenum(end_date);
catch err
    if strcmp(err.identifier, 'MATLAB:datenum:ConvertDateString')
        E.badinput('start_date and end_date must be valid date strings or date numbers')
    end
end

if start_datenum < min(GC.tVec) || start_datenum > max(GC.tVec) || end_datenum < min(GC.tVec) || end_datenum > max(GC.tVec)
    E.badinput('start_date and end_date must be within the range of dates in the GC structure');
end


xx = GC.tVec >= start_datenum & GC.tVec <= end_datenum;

n = ndims(GC.dataBlock);
permvec = [n, 1:n-1];
db = permute(GC.dataBlock,permvec);
sz = size(db);
db = db(xx,:);
db = reshape(db,[sum(xx), sz(2:end)]);
GC.dataBlock = ipermute(db,permvec);

GC.tVec = GC.tVec(xx);
xx(find(xx,1,'last')+1) = true;
GC.tEdge = GC.tEdge(xx);



end

