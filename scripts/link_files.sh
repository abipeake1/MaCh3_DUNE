#!/bin/bash
#
# A script to link required files to the proper place i.e. where the sample config files will look for them

MACH3DIR=`pwd`
#FILESDIR=/vols/t2k/users/ljw20/data/DUNE_2021/DUNE_2021_splines_tdr_v8
FILESDIR=/cvmfs/dune.osgstorage.org/pnfs/fnal.gov/usr/dune/persistent/stash/MaCh3/inputs/TDR/v2/


if [ ! -d "$MACH3DIR/inputs/DUNE_CAF_files" ]
then
  mkdir $MACH3DIR/inputs/DUNE_CAF_files
fi
ln -sf ${FILESDIR}/DUNE_2023_FD_CAFs/*root inputs/DUNE_CAF_files



if [ ! -d "$MACH3DIR/inputs/DUNE_spline_files" ]
then
  mkdir $MACH3DIR/inputs/DUNE_spline_files
fi

ln -sf ${FILESDIR}/DUNE_2023_FD_splines/*root inputs/DUNE_spline_files
