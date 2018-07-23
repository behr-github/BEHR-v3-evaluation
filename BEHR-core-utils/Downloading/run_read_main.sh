#!/bin/bash -i
# This script will figure out what SP .mat files we need to 
# keep BEHR up to date. It assumes that the record is continuous up to the last
# file is continuous, and that it only needs to run forward.  For simplicity, it
# will run up to the last day that MODIS albedo data is available for, since
# that release should be the limiting factor.
#
# Needs the environmental variables "SPDIR" (directory where the OMI_SP .mat
# files are kept), "MODDIR" (root directory for MODIS files - it will assume
# that the MCD43D* products are in the MCD43D subdirectory within there), and
# MATRUNDIR (directory where it will generate files related to running MATLAB)
# Also relies on alias startmatlab being defined to open a terminal window
# Matlab session (no GUI).
#
# Should be set to run in crontab about a day after all the downloading is
# handled, that way we can be pretty sure all the data is in.

source ~/.bashrc

# Debugging level, set higher to print more information
DEBUG=2
# How many days back in time to look for OMI_SP .mat files. Must be < 0
# as it is tested against offset back in time.
stopoffset=-90

SPDIR="$(python get_behr_path.py sp_mat_dir)"
MODDIR="$(python get_behr_path.py mcd43d_dir)"

if [[ -z $MATRUNDIR ]]
then
    echo "ERROR run_read_main.sh: env. var. MATRUNDIR unset"
    automessage.sh "ERROR run_read_main" "Env. variable MATRUNDIR unset"
    exit 1
elif [[ $DEBUG -gt 0 ]]
then
    echo -e "SPDIR=${SPDIR}\nMODDIR=${MODDIR}\nMATRUNDIR=${MATRUNDIR}"
fi

# Check that the SPDIR and MODDIR folders exist, if not, the file server
# likely isn't mounted
if [[ ! -d $SPDIR ]];
then
    echo "ERROR run_read_main.sh: $SPDIR is not a directory. Is the file server mounted?"
    automessage.sh "ERROR run_read_main.sh" "$SPDIR does not exist"
    exit 1
fi
if [[ ! -d $MODDIR ]];
then
    echo "ERROR run_read_main.sh: $MODDIR is not a directory. Is the file server mounted?"
    automessage.sh "ERROR run_read_main.sh" "$MODDIR does not exist"
    exit 1
fi


# Find the last existing OMI_SP_YYYYMMDD.mat file
offset=0
foundit=false
while true
do
    startdate=$(date -d "${offset} days" +'%Y%m%d')
    testfile="${SPDIR}/OMI_SP_*${startdate}.mat"

    if [[ $DEBUG -gt 1 ]]; then echo "Checking for $testfile"; fi

    for f in $testfile
    do
        if [[ -f $f ]]
        then
            if [[ $DEBUG -gt 0 ]]; then echo "Found $f"; fi
            offset=$((offset+1))
            startdate=$(date -d "${offset} days" +'%Y-%m-%d')
            foundit=true
            break
        fi
    done
    if $foundit
    then
        break
    elif [[ $offset -lt $stopoffset ]]
    then
        automessage.sh "run_read_main.m failed" "No OMI_SP files found within $((-stopoffset)) days."
        exit 1
    else
        offset=$((offset - 1))
    fi
done

# Find the last MCD43D31 (BRDF quality) file
offset=0
while true
do
    y=$(date -d "${offset} days" +'%Y')
    doy=$(date -d "${offset} days" +'%j')
    fpat="MCD43D31.A${y}${doy}.*.hdf"
    testfile="${MODDIR}/MCD43D31/${y}/${fpat}"

    if [[ $DEBUG -gt 1 ]]; then echo "Checking for $testfile"; fi

    ls $testfile >& /dev/null
    if [[ $? -eq 0 ]]
    then
        if [[ $DEBUG -gt 0 ]]; then echo "Found $testfile"; fi
        enddate=$(date -d "${offset} days" +'%Y-%m-%d')
        break
    elif [[ $offset -lt $stopoffset ]]
    then
        automessage.sh "run_read_main.m failed" "No MCD43D31 files found within $((-stopoffset)) days."
        exit 1
    else
        offset=$((offset-1))
    fi
done

echo "Start: $startdate End: $enddate"

# With the new behr_paths class, it is much easier to set up over SSH connections, so
# we will use the default paths in that class for inputs and outputs. We just need to 
# run the add paths method at the beginning
echo "addpath('${HOME}/Documents/MATLAB/BEHR/BEHR-core-utils/Utils/Constants');" >> ${MATRUNDIR}/runscript.m
echo "behr_paths.AddCodePaths('nosave');" >> ${MATRUNDIR}/runscript.m
echo "warning('off', 'all');" >${MATRUNDIR}/runscript.m 
echo "read_main('${startdate}', '${enddate}', 'overwrite', false, 'DEBUG_LEVEL', 1); exit(0)" >> ${MATRUNDIR}/runscript.m

startmatlab -r "run('${MATRUNDIR}/runscript.m')" > "${MATRUNDIR}/mat.log"

matexit=$?
if [[ $matexit -ne 0 ]]
then
    automessage.sh "MATLAB: read_main.m failed" -f "$MATRUNDIR/mat.log"
else
    automessage.sh "MATLAB: read_main.m succeeded" -f "$MATRUNDIR/mat.log"
fi
exit 0
