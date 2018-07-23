# BEHR-v3-evaluation
This repository contains the code used in Laughner, Zhu, and Cohen, AMT 2018 evaluating v3.0B of the BEHR OMI NO2 product. If you
just want to see how the analysis is done, look at `behr_validation_plots.m`, this is the driver script that produces the figures
and tables from the paper. It also contains a driver function that does the comparison between BEHR columns and aircraft/Pandora 
columns. This does not start from the ICARTT files for the aircraft and Pandora data, it relies on them being preproccessed into
Matlab `.mat` files by `read_merge_data.m` and `read_pandora_files.m` in `Aircraft-Satellite-Comparison/read_icart_files`. These
`.mat` files are stored in the supporting repository at https://doi.org/10.6078/D1JT28.

This code also uses `.mat` versions of BEHR files rather than the publicly available `.hdf` files. In these `.mat` files, the native 
BEHR data is represented as a structured named `Data` and the gridded data as one named `OMI`. You can generate these structures 
from the `.hdf` files using the function `behrhdf2struct` in the `BEHR-core-utils/Utils` subdirectory.

## Setup
To run this code, you will need to do several things:

1. Add all subdirectories here to your Matlab path.
1. Run `BEHR_initial_setup.m` in the `BEHR-core-utils` directory. This will create `BEHR-core-utils/Utils/Constants/behr_paths.m` 
which will contain some paths that other code relies on. (You can run `behr_paths.ValidatePaths` to check if any are invalid. Some, 
like the MODIS and OMI directories, aren't necessary, just set those to some existing directory.'). It will also build the gridding 
code, which should only be needed if you're going to redo the uncertainty analysis yourself.
1. Update the directories at the top of `misc_behr_v3_validation.m` to where you downloaded the BEHR v2.1C files and converted to Matlab 
files and the version 2 WRF profiles.
1. Update the main path in `merge_field_names.m` (in `Aircraft-Satellite-Comparison/utility_scripts/`) to where you downloaded my preprocessed
aircraft data.
1. Update the path in the subfunction `verify_pandora_dir` of `Aircraft-Satellite-Comparison/satellite_verification_scripts/run_pandora_verification.m`
to where you downloaded my preprocessed Pandora data.
