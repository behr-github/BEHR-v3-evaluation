function [ loncorn, latcorn ] = geos_chem_corners(  )
%geos_chem_corners Generates corner points for the 2 x 2.5 deg GEOS-Chem cells.
%   The sole purpose of this function is to return the 2 x 2.5 resolution
%   GEOS-Chem pixel corners. This replicates the edge points given at
%   http://acmg.seas.harvard.edu/geos/doc/man/ in Appendix A2.3
%
%   Josh Laughner <joshlaugh5@gmail.com> 21 Aug 2014

lonvec = -181.25:2.5:178.75;
latvec = [-90, -89:2:89, 90]';

[loncorn, latcorn] = meshgrid(lonvec, latvec);

end

