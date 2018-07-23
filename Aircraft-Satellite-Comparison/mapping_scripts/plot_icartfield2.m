function [  ] = plot_icartfield2( Merge_str, data_field, states, varargin )
%UNTITLED8 Summary of this function goes here
%   Detailed explanation goes here

p = inputParser;
p.addRequired('Merge_str');
p.addRequired('data_field',@isstr);
p.addRequired('states');
p.addOptional('skip',1,@isscalar);
p.addParamValue('cbrange',[],@ismatrix);

p.parse(Merge_str, data_field, states, varargin{:})
pout = p.Results;
Merge= pout.Merge_str;
field = pout.data_field;
states = pout.states;
skip = pout.skip;
cbrange = pout.cbrange;

if ~iscell(states); states = {states}; end

state_abbrev = {'al','ak','az','ar','ca','co','ct','de','fl','ga','hi','id','il','in','ia','ks','ky','la','me','md','ma','mi','mn','ms','mo','mt','ne','nv','nh','nj','nm','ny','nc','nd','oh','ok','or','pa','ri','sc','sd','tn','tx','ut','vt','va','wa','wv','wi','wy'};
state_names = {'alabama','alaska','arizona','arkansas','california','colorado','connecticut',...
    'delaware','florida','georgia','hawaii','idaho','illinois','indiana','iowa','kansas',...
    'kentucky','louisiana','maine','maryland','massachusetts','michigan','minnesota','mississippi',...
    'missouri','montana','nebraska','nevada','new_hampshire','new_jersey','new_mexico','new_york','north_carolina','north_dakota',...
    'ohio','oklahoma','oregon','pennsylvania','rhode_island','south_carolina','south_dakota','tennessee','texas','utah','vermont','virginia','washington','west_virginia',...
    'wisconsin','wyoming'};

% Read the states shapefile and plot
usa = shaperead('usastatehi.shp');
figure;
if strcmpi(states,'all');
    for a=1:numel(usa)
        plot(usa(a).X, usa(a).Y,'color','k');
        hold on
    end
else
    for b = numel(states);
        if length(states{b}) == 2
            xx = find(strcmpi(states{b},state_abbrev));
        else
            xx = find(strcmpi(states{b},state_names));
        end
        plot(usa(xx).X, usa(xx).Y,'color','k');
        hold on
    end
end

% Read the icart field, removing fills
data = eval(sprintf('Merge.Data.%s.Values',field));
lon = Merge.Data.LONGITUDE.Values - 360;
lat = Merge.Data.LATITUDE.Values;

if strcmp(field,'UTC');
    fills = NaN;
else
    fills = eval(sprintf('Merge.Data.%s.Fill',field));
end
ulod = Merge.metadata.upper_lod_flag; llod = Merge.metadata.lower_lod_flag;

yy = data == fills | data == ulod | data == llod;
data = data(~yy); lon = lon(~yy); lat = lat(~yy);

% Reduce the size of the plotted vectors
tmp = true(size(lat));
tmp(1:skip:end) = 0;

lat(tmp) = [];
lon(tmp) = [];
data(tmp) = [];

% Plot the data
scatter(lon,lat,20,data);
cb = colorbar;
if strcmp(field,'UTC');
    ylabel(cb,'s');
else
    ylabel(cb,eval(sprintf('Merge.Data.%s.Unit',field)));
end
title(sprintf('%s for %s',field,Merge.metadata.date),'fontsize',16)

end

