function [ is_overlap ] = test_box_overlap( Ax, Ay, Bx, By )
%TEST_BOX_OVERLAP Tests if boxes A and B overlap at all.
%   Ax, Ay, Bx, and By should be 4-element vectors defining the x and y
%   coordinates of the corners of two boxes. This function will return true
%   if the boxes overlap at all. Note that the inputs must go around the
%   rectangle and not crisscross or this algorithm will error. Note that if
%   any input point is a NaN this will return false.
%
%   4 cases must be considered:
%       - A and B do not overlap at all
%       - A and B partially overlap
%       - A is totally contained within B
%       - B is totally contained within A
%

E = JLLErrors;

%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% INPUT CHECKING %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%

if any([numel(Ax), numel(Ay), numel(Bx), numel(By)] ~= 4) || ~isvector(Ax) || ~isvector(Ay) || ~isvector(Bx) || ~isvector(By)
    E.badinput('All inputs must be four element vectors')
elseif ~isnumeric(Ax) || ~isnumeric(Ay) || ~isnumeric(Bx) || ~isnumeric(By)
    E.badinput('All inputs must be numeric')
end

if any(isnan(Ax)) || any(isnan(Ay)) || any(isnan(Bx)) || any(isnan(By))
    is_overlap = false;
    return
end

% Test if the points are ordered correctly by checking that the cross
% products of successive edge vector all have the same sign.
% (http://debian.fmi.uni-sofia.bg/~sergei/cgsr/docs/clockwise.htm)
Av = nan(4,3);
Bv = nan(4,3);
for i=1:4
    i2=mod(i,4)+1; % following point with wrapping
    Av(i,:) = [Ax(i2), Ay(i2), 0] - [Ax(i), Ay(i), 0];
    Bv(i,:) = [Bx(i2), By(i2), 0] - [Bx(i), By(i), 0];
end
Across = nan(4,3);
Bcross = nan(4,3);
for i=1:4
    i2=mod(i,4)+1; % following point with wrapping
    Across(i,:) = cross(Av(i,:),Av(i2,:));
    Bcross(i,:) = cross(Bv(i,:),Bv(i2,:));
end
if ~all(Across(:,3) > 0) && ~all(Across(:,3) < 0)
    E.badinput('The quadrilateral A does not have its points defined correctly; successive points should not criss-cross')
elseif ~all(Bcross(:,3) > 0) && ~all(Bcross(:,3) < 0)
    E.badinput('The quadrilateral B does not have its points defined correctly; successive points should not criss-cross')
end

%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% MBIN FUNCTION %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%

is_overlap = false;

% First we check if we can immediately say the two boxed do not overlap. We
% will do so with simple logical tests. Consider two rectangle aligned with
% the x and y axes which each of the quadrilaterals given are inscribed in.
if min(Ax) > max(Bx) || min(Bx) > max(Ax) || min(Ay) > max(By) || min(By) > max(Ay)
    return
end

% Second we can check if any corner of A is contained in B and vice versa.
% This will handle both cases where A is totally contained in B and vice
% versa as well as most cases of partial overlap. By testing each rectangle
% first and returning if it has a point within the other, we can gain a
% slight performance boost by not carrying out additional tests.

aa = inpolygon(Ax,Ay,Bx,By);
if any(aa)
    is_overlap = true;
    return
end
bb = inpolygon(Bx,By,Ax,Ay);
if any(bb)
    is_overlap = true;
    return
end

% Finally we will check for the rare case where neither has a point within
% the other but nevertheless do overlap. We will check this by testing if
% the intersection between any combination of sides in A and B intersects.
Am = nan(4,1);
Bm = nan(4,1);
Ab = nan(4,1);
Bb = nan(4,1);
for i=1:4
    i2=mod(i,4)+1; % following point with wrapping
    Am(i) = (Ay(i2) - Ay(i)) / (Ax(i2) - Ax(i));
    Bm(i) = (By(i2) - By(i)) / (Bx(i2) - Bx(i));
    Ab(i) = Ay(i) - Am(i) * Ax(i);
    Bb(i) = By(i) - Bm(i) * Bx(i);
end
for i=1:4
    i2 = mod(i,4)+1;
    for j=1:4
        j2 = mod(j,4)+1;
        % same slope --> parallel, do not intersect
        if Am(i) == Bm(j)
            continue
        end
        % otherwise check intersection lies on the line segment for both
        % edges
        x_int = (Bb(j) - Ab(i)) / (Am(i) - Bm(j));
        if x_int >= min(Ax([i,i2])) && x_int <= max(Ax([i,i2])) && x_int >= min(Bx([j,j2])) && x_int <= max(Bx([j,j2]))
            is_overlap = true;
            return
        end
    end
end


end

