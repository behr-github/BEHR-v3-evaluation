function [ Data ] = behrhdf2struct( hdfi )
%Parses a BEHR HDF file (actually an HDF5 or he5 file) into a Matlab data structure
%   Useful for comparing old values in BEHR HDF files to current Matlab
%   data structures. Especially useful in conjunction with
%   "Compare_Data_Fields", which allows you to compare an aribitrary field
%   in two data structures
%
%   DATA = BEHRHDF2STRUCT( HDFI ) returns a structure, DATA with all of the datasets
%   from the HDF5 file represented by HDFI as fields. Each orbit will be one element 
%   of DATA. HDFI is an info structure returned by H5INFO() on the HDF5 file to be
%   read.
%
%   Josh Laughner <joshlaugh5@gmail.com> 16 Apr 2014

Data=struct([]);
n = length(hdfi.Groups(1).Groups);
for i_orbit=1:n
    
    k = strfind(hdfi.Groups(1).Groups(i_orbit).Name,'Swath');
    swath = str2double(hdfi.Groups(1).Groups(i_orbit).Name(k+5:end));
    Data(i_orbit).Swath = swath;
    
    for i_dat=1:length(hdfi.Groups(1).Groups(i_orbit).Datasets)
        dset = h5read(hdfi.Filename, h5dsetname(hdfi,1,i_orbit,i_dat));
        fill_val = hdfi.Groups(1).Groups(i_orbit).Datasets(i_dat).FillValue;
        field_name = hdfi.Groups(1).Groups(i_orbit).Datasets(i_dat).Name;
        if fill_val < 0
            dset(dset < fill_val * 0.98) = nan; % allow for a little bit of floating point error in the fill value
        else
            dset(dset == fill_val) = nan; % if the fill value is positive, it's probably a flag field, so there shouldn't be floating point error
        end
        Data(i_orbit).(field_name) = dset;
    end
   
    % Also handle swath-level attributes. Try to convert each to a number,
    % if that returns a NaN, assume it should stay as a string (which is
    % most of them)
    for i_att = 1:length(hdfi.Groups(1).Groups(i_orbit).Attributes)
        attribute = hdfi.Groups(1).Groups(i_orbit).Attributes(i_att);
        att_value = str2double(attribute.Value);
        if isnan(att_value)
            att_value = attribute.Value;
        end
        % Some of the attributes get added directly to the .hdf files
        % during publication, and so aren't actually valid structure field
        % names. Fix that here
        att_name = strrep(attribute.Name, '-', '_');
        Data(i_orbit).(att_name) = att_value;
    end
end

end

