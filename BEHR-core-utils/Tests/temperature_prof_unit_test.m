function [ success ] = temperature_prof_unit_test( test_data_input_dir )
%TEMPERATURE_PROF_UNIT_TEST Tests rNmcTmp2 against previous output
%   SUCCESS = TEMPERATURE_PROF_UNIT_TEST() Reads the file
%   "temperature_prof_test_data.mat" from the same directory as this
%   function and runs rNmcTmp2 with the latitude, longitude, and month data
%   stored in that file. It compares the resulting temperature profiles
%   against those stored in that mat file, they must be equal according to
%   ISEQUALN() in order for the test to pass. SUCCESS will be true if all
%   temperature profiles are equal, false otherwise.
%
%   TEMPERATURE_PROF_UNIT_TEST( TEST_DATA_INPUT_DIR ) will generate
%   temperature_prof_test_data.mat from the OMI_SP*.mat files found in
%   TEST_DATA_INPUT_DIR. I generally recommend that the directory be a unit
%   tests data directory in the BEHR-core repository, since those have a
%   good selection of days.

save_dir = fileparts(mfilename('fullpath'));
save_file = fullfile(save_dir, 'temperature_prof_test_data.mat');

if exist('test_data_input_dir', 'var')
    make_test_data(test_data_input_dir, save_file);
else
    success = run_test(save_file);
end

end

function make_test_data(data_dir, save_file)
F = dirff(fullfile(data_dir, 'OMI_SP*.mat'));
n = min(5, numel(F)); % limit to 5 files to keep the test data to something we can store on GitHub (< 100 MB)

lon = cell(n,1);
lat = cell(n,1);
mon = cell(n,1);
temp = cell(n,1);

id_str = sprintf('Produced from %d files:\n\t%s', numel(F), strjoin({F.name},'\n\t')); %#ok<NASGU>

fileTmp = fullfile(behr_paths.amf_tools_dir,'nmcTmpYr.txt');

for a=1:n
    fprintf('Calculating temperatures for file %d of %d\n', a, n);
    D = load(F(a).name);
    Data = D.Data;
    lon{a} = cat_sat_data(Data, 'Longitude');
    lat{a} = cat_sat_data(Data, 'Latitude');
    mon{a} = repmat(month(Data(1).Date), size(lon{a})); % each file is one day, so the month had better be the same for each orbit
    temp{a} = rNmcTmp2(fileTmp, behr_pres_levels, lon{a}, lat{a}, mon{a});
end

save(save_file, 'id_str', 'lon', 'lat', 'mon', 'temp');

end

function success = run_test(data_file)
D = load(data_file);
fileTmp = fullfile(behr_paths.amf_tools_dir,'nmcTmpYr.txt');
success = true;
for a=1:numel(D.lon)
    fprintf('Testing day %d of %d\n', a, numel(D.lon));
    test_temp = rNmcTmp2(fileTmp, behr_pres_levels, D.lon{a}, D.lat{a}, D.mon{a});
    if ~isequaln(test_temp, D.temp{a})
        success = false;
        return
    end
end
end