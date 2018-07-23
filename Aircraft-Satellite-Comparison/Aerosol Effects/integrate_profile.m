function [ column_struct ] = integrate_profile(  )
%integrate_profile Numerically integrate an aircraft profile
%   Carries out a trapezoidal integration of the specified field in a Merge
%   structure for each profile.  Returns a structure with the integrated
%   column, the profile number, and the profile bins (for comparison).
%
%   Josh Laughner <joshlaugh5@gmail.com> 23 Sept 2014

start_date = '01/16/2013';
end_date = '02/16/2013';

mat_dir = '/Volumes/share/GROUP/DISCOVER-AQ/Matlab Files/Aerosol_Merges/Properties/';
mat_pattern = 'Cal*%s_%s_%s.mat';
int_fieldname = 'Ext_green_ambient_TSIPSAP_LARGE';
alt_fieldname = 'Pressure_Altitude';
profnum_fieldname = 'ProfileNumber';

column_struct = struct('Dates',sprintf('%s to %s',start_date,end_date),'IntegratedField',int_fieldname,'IntegratedFieldUnits','','AltFieldUnits','','ProfNums',[],'IntegratedColumns',[],'Bins',{{}},'BinAlts',{{}});

dates = datenum(start_date):datenum(end_date);
units_unset = 1;
N=0;
for d=1:numel(dates)
    S=wildcard_load(mat_dir,mat_pattern,dates(d));
    if isempty(S); continue; end
    Merge = S.Merge;
    
    profnums = remove_merge_fills(Merge,profnum_fieldname);
    unique_profnums = unique(profnums(profnums>0));
    int_field = remove_merge_fills(Merge,int_fieldname);
    alt = remove_merge_fills(Merge,alt_fieldname);
    
    % The first time through, copy the units
    if units_unset
        column_struct.IntegratedFieldUnits = Merge.Data.(int_fieldname).Unit;
        column_struct.AltFieldUnits = Merge.Data.(alt_fieldname).Unit;
        units_unset = 0;
    end
    
    % Loop through all the profiles for a given day
    for p=1:numel(unique_profnums)

        
        % Find the data for this profile
        xx = profnums == unique_profnums(p);
        
        % Bin the column
        [bins, bin_alts] = bin_vertical_profile(alt(xx),int_field(xx),0.25);
        
        % Fill in any internal nans, and clip off those at the edge. If
        % only one value is not a nan, skip this profile
        if sum(~isnan(bins))<2; continue; end
        [bin_alts, bins] = fill_nans(bin_alts, bins);
        
        % Integrate the column using the trapezoidal rule
        int_col = trapz(bin_alts, bins);
        
        N=N+1;
        
        % Save everything to the structure
        column_struct.ProfNums(N) = unique_profnums(p);
        column_struct.IntegratedColumns(N) = int_col;
        column_struct.Bins{N} = bins;
        column_struct.BinAlts{N} = bin_alts;
    end
    
end

end

