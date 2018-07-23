function [ bin_centers, bin_edges ] = geos5_reduced_Pgrid( surfPres )
%geos5_reduced_Pgrid Returns the pressure bin centers and edges for the
%reduced (47 level) GEOS-5 pressure scheme.
%   The GEOS-5 pressure coordinate is a hybrid system that uses a sigma
%   (terrain following) coordinate for the lowest levels and fixed pressure
%   levels above that.  This takes a surface pressure in hPa and returns
%   the bin centers and bin edges in hPa.

% AP is in hPa, BP is unitless. 

AP = [0.000000e+00, 4.804826e-02, 6.593752e+00, 1.313480e+01,...
1.961311e+01, 2.609201e+01, 3.257081e+01, 3.898201e+01,...
4.533901e+01, 5.169611e+01, 5.805321e+01, 6.436264e+01,...
7.062198e+01, 7.883422e+01, 8.909992e+01, 9.936521e+01,...
1.091817e+02, 1.189586e+02, 1.286959e+02, 1.429100e+02,...
1.562600e+02, 1.696090e+02, 1.816190e+02, 1.930970e+02,...
2.032590e+02, 2.121500e+02, 2.187760e+02, 2.238980e+02,...
2.243630e+02, 2.168650e+02, 2.011920e+02, 1.769300e+02,...
1.503930e+02, 1.278370e+02, 1.086630e+02, 9.236572e+01,...
7.851231e+01, 5.638791e+01, 4.017541e+01, 2.836781e+01,...
1.979160e+01, 9.292942e+00, 4.076571e+00, 1.650790e+00,...
6.167791e-01, 2.113490e-01, 6.600001e-02, 1.000000e-02];

BP = [1.000000e+00, 9.849520e-01, 9.634060e-01, 9.418650e-01,...
9.203870e-01, 8.989080e-01, 8.774290e-01, 8.560180e-01,...
8.346609e-01, 8.133039e-01, 7.919469e-01, 7.706375e-01,...
7.493782e-01, 7.211660e-01, 6.858999e-01, 6.506349e-01,...
6.158184e-01, 5.810415e-01, 5.463042e-01, 4.945902e-01,...
4.437402e-01, 3.928911e-01, 3.433811e-01, 2.944031e-01,...
2.467411e-01, 2.003501e-01, 1.562241e-01, 1.136021e-01,...
6.372006e-02, 2.801004e-02, 6.960025e-03, 8.175413e-09,...
0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00,...
0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00,...
0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00,...
0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00];

bin_edges = AP + BP*surfPres;
bin_centers = (bin_edges(1:end-1) + bin_edges(2:end))/2;

end

