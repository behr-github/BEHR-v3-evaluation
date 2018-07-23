function [ SatDataStruct ] = calc_sat_ndens( SatDataStruct )
%CALC_SAT_NDENS Calculate the number density of boxes in GC ND51 output
%   The GEOS-Chem ND51 diagnostic does not offer the ability to include
%   number density in the output, however it can be either copied from the
%   regular 24-hr ctm.bpch file, or it could be calculated from the
%   temperature and pressure of each box using the ideal gas law:
%       
%       N / V = (n * A_v) / V = (P * A_v) / (R * T)
%
%   This function will take a structure of the output from read_gc_nd51.py
%   (run the resulting SatData through convert_py_nd51_struct.m using
%   cleanup_py_nd51(SatData{:}) to get this). It will see if the necessary
%   variables are present and, if so, calculate the number density,
%   appending it as an additional entry in the structure.
%
%   Josh Laughner <joshlaugh5@gmail.com> 2 Jul 2015

E = JLLErrors;

%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% INPUT CHECKING %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~isstruct(SatDataStruct)
    E.badinput('SatDataStruct must be a structure. Use convert_py_nd51 if necessary to convert the SatData cell array');
end

% Check that the proper variables are present in the structure.
req_vars = {'PEDGE-$_PSURF', 'DAO-3D-$_TMPU'};
present_vars = {SatDataStruct(:).fullName};
vars_bool = ismember(req_vars, present_vars);
if any(~vars_bool)
    E.badinput('The variables %s must be present in SatDataStruct, but are not', strjoin(req_vars(~vars_bool),', '));
end

%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% MAIN FUNCTION %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%

% Define constants
Av = 6.02e23; % molec./mol
R = 8.314; % J/(mol*K)

% Find which index refers to temperature and which to pressure
pres_ind = ~iscellcontents(strfind(present_vars, 'PEDGE-$_PSURF'), 'isempty');
temp_ind = ~iscellcontents(strfind(present_vars, 'DAO-3D-$_TMPU'), 'isempty');

% Calculate the number density from temperature and pressure.
T = SatDataStruct(temp_ind).dataBlock;
P = SatDataStruct(pres_ind).dataBlock;

ndens = (P *1e2 * Av) ./ (R * T);   % the factor of 1e2 accounts for pressure in hPa (need Pa)

i = numel(SatDataStruct);

SatDataStruct(i+1).dataBlock = ndens;
SatDataStruct(i+1).dataUnit = 'molec./m3'; % keep in same units as GC number density to be consistent
SatDataStruct(i+1).fullName = 'Number density of air (calculated)';
SatDataStruct(i+1).fullCat = [];
SatDataStruct(i+1).tVec = SatDataStruct(i).tVec;
SatDataStruct(i+1).modelName = SatDataStruct(i).modelName;
SatDataStruct(i+1).modelRes = SatDataStruct(i).modelRes;
SatDataStruct(i+1).dataScale = SatDataStruct(i).dataScale;
SatDataStruct(i+1).molMass = nan;
SatDataStruct(i+1).tEdge = SatDataStruct(i).tEdge;

end

