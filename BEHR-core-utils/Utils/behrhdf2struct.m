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
for a=1:n
    
    k = strfind(hdfi.Groups(1).Groups(a).Name,'Swath');
    swath = str2double(hdfi.Groups(1).Groups(a).Name(k+5:end));
    Data(a).Swath = swath;
    
    for b=1:length(hdfi.Groups(1).Groups(a).Datasets)
        dset = h5read(hdfi.Filename, h5dsetname(hdfi,1,a,b));
        fill_val = hdfi.Groups(1).Groups(a).Datasets(b).FillValue;
        field_name = hdfi.Groups(1).Groups(a).Datasets(b).Name;
        if fill_val < 0
            dset(dset < fill_val * 0.98) = nan; % allow for a little bit of floating point error in the fill value
        else
            dset(dset == fill_val) = nan; % if the fill value is positive, it's probably a flag field, so there shouldn't be floating point error
        end
        Data(a).(field_name) = dset;
    end
    
end

end

