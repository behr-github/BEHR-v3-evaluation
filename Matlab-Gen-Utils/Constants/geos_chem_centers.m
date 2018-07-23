function [ lons, lats ] = geos_chem_centers( grid )
%geos_chem_centers Generates center points for GEOS-Chem cells.
%   This function returns the center longitude/latitude of GEOS-Chem model
%   cells for several resolutions.  If no argument is passed, the user is
%   asked interactively which resolution/nest is desired. Alternatively,
%   the resolution can be specified as a string (e.g. '2x25' or
%   '025x03125NA') or by passing the size of a data block imported from the
%   model. This replicates the center points given at
%   http://acmg.seas.harvard.edu/geos/doc/man/ in Appendix A2
%
%   Josh Laughner <joshlaugh5@gmail.com> 17 Jul 2015

E = JLLErrors;

lons = [];
lats = [];

if nargin < 1
    question = 'Select a grid\n  1 - 2 x 2.5\n  2 - 0.25 x 0.3125 NA\n ?: ';
    while true
        answer = input(question,'s');
        switch answer
            case '1'
                grid = '2x25';
                break
            case '2'
                grid = '025x03125NA';
                break
            otherwise
                fprintf('Selection not recognized\n');
                return
        end
    end
elseif isnumeric(grid)
    % This handles if the user passes the size vector of a matrix. This is
    % really intended to be used automatically where the user can call it
    % with size(M) where M is some data matrix. The first two sizes should
    % correspond to one of the resolutions/nests.
    if numel(grid) < 2
        E.badinput('If passing grid as a size vector, it must have at least 2 elements')
    elseif all(grid(1:2) == [144, 91])
        grid = '2x25';
    elseif all(grid(1:2) == [225, 202])
        grid = '025x03125NA';
    else
        E.badinput('Size not recognized');
    end
elseif ischar(grid)
    % if a case is passed as a string, remove any decimal points,
    % underscores, or dashes that the user might accidentally put in.
    grid = regexprep(grid,'[._-]','');
end

switch grid
    case '2x25'
        lons = -180:2.5:177.5;
        lats = [-89.5,-88:2:88,89.5];
    case '025x03125NA'
        lons = -130:0.3125:-60;
        lats = 9.75:0.25:60;
end

end

