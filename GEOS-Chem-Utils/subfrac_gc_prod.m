function [ subfrac ] = subfrac_gc_prod( Prod, frac_cat, denom_cats, verbose )
%SUBFRAC_GC_PROD Calculates a fraction of NO production from GC output
%   GEOS-Chem v. 9-02b has 8 categories for the production of NO: aircraft,
%   anthropogenic, biomass burning, biofuels, fertilizer, lightning, soil,
%   and stratospheric.  This function takes 3 arguments:
%       Prod - a structure output from read_geos_output that contains
%       production information.
%
%       frac_cat - the index of which category to want to calculate the
%       total fraction of emissions for.
%
%       denom_cats - which categories you want to use as the sum to
%       calculate the fraction of. Should be a vector.
%
%   This function will calculate the total production over all vertical
%   levels and so return a lon x lat x time matrix with the fraction of
%   frac_cat that makes up NO production.  For example, if you wished to
%   calculate the subfraction of NO production that's due to lightning out
%   of lightning, anthropogenic, and biomass burning, then:
%
%       subfrac_gc_prod(Prod, 6, [2,3]); 
%
%   and
%
%       subfrac_gc_prod(Prod, 6, [2,3,6]); 
%
%   will both produce the same result, as this function will always include
%   frac_cat in the denominator.
%
%   One optional argument: setting the final argument to false will prevent
%   printing out the category names as a double check.
%
%   Josh Laughner <joshlaugh5@gmail.com> 18 Nov 2015

E = JLLErrors;

%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% INPUT CHECKING %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~isstruct(Prod) || any(~isfield(Prod,{'dataBlock','fullCat','tVec'}))
    E.badinput('Prod should be a structure containing the fields (at least) dataBlock, tVec, and fullCat, i.e. it should be output from read_geos_output');
end

if ~isscalar(frac_cat) || ~isnumeric(frac_cat) || frac_cat > numel(Prod)
    E.badinput('frac_cat must be a valid numeric index for the Prod structure')
end

if ~isnumeric(denom_cats) || any(denom_cats > numel(Prod))
    E.badinput('denom_cats must contain valid numeric indices for the Prod structure')
end

if ~exist('verbose','var')
    verbose = true;
end

%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% MAIN FUNCTION %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%

% Make sure that the frac_cat is not repeated in the denom_cats as it will
% be added in automatically.

xx = frac_cat == denom_cats;
if sum(xx) > 0
    denom_cats(xx) = [];
end

% Print out the category names to double check
if verbose
    fprintf('Calculating fraction of %s with the other categories:\n',Prod(frac_cat).fullCat);
    for a=1:numel(denom_cats)
        fprintf('\t %s\n',Prod(denom_cats(a)).fullCat);
    end
end

% For each category, read it in and sum over the vertical levels.
frac_cat_prod = sum_vert_levels(Prod(frac_cat));
denom_cat_prod = frac_cat_prod; % go ahead and start with that as the total
for a=1:numel(denom_cats)
    tmp_prod = sum_vert_levels(Prod(denom_cats(a)));
    denom_cat_prod = nansum2(cat(4,denom_cat_prod,tmp_prod),4);
end

subfrac = frac_cat_prod ./ denom_cat_prod;

end

function prod_mat = sum_vert_levels(Prod_i)
% Subfunction to handle summing over vertical levels if necessary. Only
% pass one element of the Prod structure.
E = JLLErrors;
if ~isscalar(Prod_i) || ~isstruct(Prod_i)
    E.badinput('Only input one element of the Prod structure')
end

% Figure out if there are vertical levels or not.  If not, the third
% dimension will be time, and that will cause issues because we don't want
% to sum over time.
ndim_data = ndims(Prod_i.dataBlock);
ntimes = numel(Prod_i.tVec);
if ndim_data > 3 || (ntimes == 1 && ndim_data > 2)
    prod_mat = squeeze(nansum2(Prod_i.dataBlock,3));
else
    prod_mat = Prod_i.dataBlock;
end

end

