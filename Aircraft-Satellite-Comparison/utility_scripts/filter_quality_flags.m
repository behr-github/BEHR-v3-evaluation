function [ xx ] = filter_quality_flags( q_flags, bits, bitvalues )
%filter_quality_flags Find profiles with the desired quality flags
%   sprial_verification_avg_pix2prof outputs quality flags as binary
%   numbers, indicating various problems of note in the processing. This
%   function takes:
%       1) A cell array or vector of quality flags (should be some uint
%       type). Currently this will only allow cell arrays to have a single
%       quality flag per cell.
%       2) The bit numbers to examine, as a vector
%       3) The values those bits should have.
%   This will return a logical matrix the same size as q_flags indicating
%   if each flag met the specified criteria.
%
%   Josh Laughner <joshlaugh5@gmail.com> 26 Mar 2015

E = JLLErrors;

%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% INPUT CHECKING %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%

if iscell(q_flags)
    if ~all(iscellcontents(q_flags,'isscalar')) || ~all(iscellcontents(q_flags,'isnumeric'))
        E.badinput('"q_flags" must contain only scalar numbers');
    end
    q_flags = cell2mat(q_flags);
end

if numel(bits) ~= numel(bitvalues)
    E.numelMismatch('bits','bitvalues')
end

%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% MAIN FUNCTION %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%

xx = true(size(q_flags));

for a=1:numel(q_flags)
    for b=1:numel(bits)
        if xx(a)
            this_bit = bitget(q_flags(a),bits(b));
            xx(a) = this_bit == bitvalues(b);
        end
    end
end

end

