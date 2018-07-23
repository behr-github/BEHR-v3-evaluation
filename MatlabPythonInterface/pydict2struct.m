function [ S ] = pydict2struct( pdict )
%PYDICT2STRUCT Convert a Python dictionary or list of dicts to a Matlab structure
%   S = PYDICT2STRUCT ( PDICT ) Converts PDICT to a Matlab structure, S. If
%   PDICT is a list of dictionaries, then S will be a multi-element
%   structure. Converts fields using python2matlab.m.


if isa(pdict, 'py.list')
    pcell = cell(pdict);
    dict_test = true(1, numel(pcell));
    for a = 1:numel(pcell)
        dict_test(a) = isa(pcell{a}, 'py.dict');
    end
    if ~all(dict_test)
        error('pyinterface:badinput','PDICT must be of type PY.LIST or PY.DICT. If is a PY.LIST, it must be a list of dictionaries.')
    end
elseif isa(pdict, 'py.dict')
    pcell = {pdict};
else
    error('pyinterface:badinput','PDICT must be of type PY.LIST or PY.DICT');
end

S = repmat(struct, 1, numel(pcell));

for a = 1:numel(pcell)
    stemp = struct(pcell{a});
    fns = fieldnames(stemp);
    for b = 1:numel(fns)
        S(a).(fns{b}) = python2matlab(stemp.(fns{b}));
    end
end

end
