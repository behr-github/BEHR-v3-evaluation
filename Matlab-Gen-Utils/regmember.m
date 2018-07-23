function [ matches ] = regmember( test_cell, expr_cell )
%REGMEMBER Regular expression test similar to ISMEMBER()
%   MATCHES = REGMEMBER( TEST_CELL, EXPR_CELL ) returns a logical array the
%   same size as TEST_CELL that is true for any member of TEST_CELL that
%   match any of the regular expressions in EXPR_CELL. Both inputs must be
%   cell arrays of strings.

E = JLLErrors;
if ~iscellstr(test_cell)
    E.badinput('TEST_CELL must be a cell array of strings');
end
if ~iscellstr(expr_cell)
    E.badinput('EXPR_CELL must be a cell array of strings');
end

matches = false(size(test_cell));
for a=1:numel(expr_cell)
    reg_result = regcmp(test_cell, expr_cell{a});
    matches = matches | reg_result;
end

end

