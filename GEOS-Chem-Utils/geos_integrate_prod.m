function [ total_prod, db ] = geos_integrate_prod( varargin )
%GEOS_INTEGRATE_PROD Calculate total production of a category
%   This will calculate the total number of MOLES (not molecules) of a
%   species produced in GEOS-Chem. It can take arguments as:
%
%   GEOS_INTEGRATE_PROD(Area, TimePer1, TimePer2,...) requires that all
%   inputs be scalar structures output from read_geos_output. The first one
%   must contain the 'Grid box surface area' tracer, and all the following
%   ones must have the same production tracer, currently it is required
%   that the units be molec/cm2/s. Multiple time periods for production can
%   be given so that it is convinient if your GEOS-Chem output is split up
%   into several .bpch files - this way, you don't need to try to
%   concatenate multiple structures before using this function.
%
%
%   Josh Laughner <joshlaugh5@gmail.com> 10 Aug 2015


E = JLLErrors;

%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% INPUT CHECKING %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%

req_fields = {'dataBlock','dataUnit','fullCat','fullName','tEdge'};


if numel(varargin) < 2
    E.badinput('At least two input structures (area and emissions) are required')
end

% Check that all inputs are scalar structures
for a=1:numel(varargin)
    if ~isstruct(varargin{a}) || ~isscalar(varargin{a}) || any(~ismember(req_fields, fieldnames(varargin{a})))
        E.badinput('All inputs must be scalar structures with the fields %s',strjoin(req_fields,', '));
    end
end

% Make sure that the first input is the area of the GEOS-Chem grid cells
AreaStruct = varargin{1};

if isempty(regexp(AreaStruct.fullName,'surface area', 'once'))
    E.badinput('The first input must contain a surface area variable')
elseif size(AreaStruct.dataBlock, 3) > 1
    warning('Only the first 2D slice of the Area dataBlock will be used, since it is assumed that the areas of the GEOS-Chem grid cells are static')
    AreaStruct.dataBlock(:,:,2:end) = [];
elseif isempty(regexp(AreaStruct.dataUnit, 'm\^?2','once'))
    E.badinput('Area is not in square meters.')
end

area = AreaStruct.dataBlock * 1e4; % convert to cm^2 (the emissions is assumed to be in molec/cm^2/s)

% Check the rest of the input structures: they should all have the same
% fullName and fullCat and the first two dimensions of the dataBlock are
% the same as the area matrix.  Also the units should be in molec/cm2/s

Emis = [varargin{2:end}]; % convert to single structure
for a=1:numel(Emis)
    if a>1 && (~strcmp(Emis(a).fullCat, Emis(1).fullCat) || ~strcmp(Emis(a).fullName, Emis(1).fullName))
        E.badinput('Input structure %d seems to refer to a different tracer than the second one (first emissions structure', a+1)
    end
    
    if size(Emis(a).dataBlock,1) ~= size(area,1) || size(Emis(a).dataBlock,2) ~= size(area,2)
        E.badinput('Input structure %d seems to have a different lat/lon size of dataBlock than the area input')
    elseif isempty(regexp(Emis(a).dataUnit, 'molec\.?/cm\^?2/s','once'))
        E.badinput('Emissions are required to have units of molec/cm2/s')
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% MAIN FUNCTION %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%

area = reshape(area,[],1);

% Total emissions for each cell = emission rate * area * time step. Sum up
% each structure and add it to the running total.

total_prod = 0;
area_list = [];
db.prod = [];
db.n_sec = [];
db.tvec = [];
db.prod_grid = [];


for a=1:numel(Emis)
    % Reshape so that all the spatial coordinates are strung out along the
    % first dimension. This saves having to check if there's multiple
    % levels or just one (which would only have 3 dims instead of 4).
    
    % tEdge has an extra entry over the number of time periods, since each
    % needs a start and end.
    n_times = length(Emis(a).tEdge) - 1;
    sz = size(Emis(a).dataBlock);
    this_emis = reshape(Emis(a).dataBlock,[],n_times);
    
    n_levels = size(this_emis,1) / numel(area);
    if mod(n_levels,1) ~= 0
        E.unknownError('Somehow there is not a full set of grid cells for each level in input number %d', a+1)
    end
    
    this_area = repmat(area, n_levels, 1);
    area_list = cat(2, area_list, this_area);
    
    for b=1:n_times
        t_vec = datevec(Emis(a).tEdge(b:b+1));
        n_sec = etime(t_vec(2,:), t_vec(1,:));
        
        prod_slice = this_emis(:,b) .* this_area * n_sec;
        total_prod = total_prod + nansum(prod_slice, 1);
        db.prod = cat(1,db.prod,total_prod);
        db.n_sec = cat(1,db.n_sec,n_sec);
        db.tvec = cat(1, db.tvec, t_vec');
        
        prod_slice = nansum(reshape(prod_slice, sz(1:3)),3);
        if isempty(db.prod_grid)
            db.prod_grid = prod_slice;
        else
            db.prod_grid = nansum(cat(3, db.prod_grid, prod_slice),3);
        end
    end
    
end

% convert molecules to moles
total_prod = total_prod / 6.022e23;
db.prod = db.prod / 6.022e23;
db.prod_grid = db.prod_grid / 6.022e23;

end

