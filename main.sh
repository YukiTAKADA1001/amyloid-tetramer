#!/bin/bash -eu

readonly HOME_DIR=/home/users/bip
readonly SAVE_DIR=/save/users/bip
readonly PROJ_DIR=${SAVE_DIR}/simulation4
readonly PROGRAMS_DIRNAME=programs


function usage() {
cat <<_EOT_
Usage:
    $0 acd|nat loop

Description:
    This script does following things
    - [tleap]   Make extended structure of Amyloid beta
    - [OpenMM]  Equil the peptide without water or metal ions
    - [VMD]     Make system by integrate pdb files
    - [tleap]   Add water and metal ions to the system
    - [Python]  Add chain ID to distinguish peptides
    - [OpenMM]  Execute equilibration and production run

Options:
    This script has no option.

_EOT_
exit 1
}



# Check usage
if [ $# -ne 2 ]; then
    usage
    exit 1
fi
# Exit if parameter condition was not acd or nat.
CONDITION_DIRNAME=$1
if [ $CONDITION_DIRNAME != "acd" ] && [ $CONDITION_DIRNAME != "nat" ]; then
    echo "Error:"
    echo "Specify condition directory as acd or nat."
    exit 1
fi
# Exit if paremeter loop was not between 0 and 100.
LOOP=$2
if [ $LOOP -lt 1 ] || [ $LOOP -gt 100 ]; then
    echo "Error:"
    echo "Set loop between 0 and 100."
    exit 1
fi



# Main operation
# PROJ_DIR以下は完全に閉じたディレクトリ構造（いかなる絶対パスも参照しない）
# のため、大元のPROJ_DIRなどの絶対パスを変更しても問題なし。
cd $PROJ_DIR
echo "Preparing for production run..."

for ((i=0; i < $LOOP; i++))
do
    cd $CONDITION_DIRNAME


    # Make new directory for new simulation. (Home and Save area)
    sim_dirname=simulation_`date +%Y%m%d_%H%M%S`
    mkdir $sim_dirname && cd $sim_dirname
    mkdir prepare && cd prepare


    # Make extended structure (tleap script)
    # stdout is tleap log. Can get them in leap.log, so throw stdout away to /dev/null
    tleap -f ${PROJ_DIR}/${PROGRAMS_DIRNAME}/make_extended_${CONDITION_DIRNAME}.in > /dev/null


    # Execute long time process as background
    {
        # Equil peptide without water and ions
        jsub_stdout=`jsub -q PN ${PROJ_DIR}/${PROGRAMS_DIRNAME}/equil_only_pep.sh`
        jobid=`echo $jsub_stdout | python -c "print( input().split('.')[0] )"`


        # Wait until equil simulation finishes.
        # stdout file is the finish flag.
        while [ ! -e equil_only_pep.sh.o${jobid} ]
        do
            sleep 10s
        done
        # Modify atom index (N-terminus H -> H1)
        # Without this operation, some errors occur with tleap.
        # This is attributed to OpenMM PDB output format.
        for f in amyloid*.pdb
        do
            sed -i -e 's/ATOM      2  H   ASP A   1/ATOM      2  H1  ASP A   1/g' $f
            if [ $CONDITION_DIRNAME = "acd" ]; then
                sed -i -e 's/HIS/HIP/g' $f
            fi
        done


        while :
        do
            ####
            # Configure system using VMD
            # Move peptide
            vmd -dispdev none -e ${PROJ_DIR}/${PROGRAMS_DIRNAME}/config_system.tcl > vmd.log
            # Integrate sigle peptide files into one pdb.
            if [ -e initial.pdb ]; then
                rm initial.pdb
            fi
            for f in amyloid_transed*.pdb
            do
                sed -e 1d -e 's/END/TER/g' $f >> initial.pdb
            done
            echo "END" >> initial.pdb


            # Add water and ions using Amber tleap
            tleap -f ${PROJ_DIR}/${PROGRAMS_DIRNAME}/add_solvent_ion_${CONDITION_DIRNAME}.in > /dev/null


            # Add chain ID to distinguish peptides
            python ${PROJ_DIR}/${PROGRAMS_DIRNAME}/add_chnid.py


            # Check if distances between peptides are not too close.
            # If OK, this python script returns status 0, and exit loop.
            # Temporary invalid -e option to use error code 1 as loop continue
            # flag ($dist_check).
            set +e
            python ${PROJ_DIR}/${PROGRAMS_DIRNAME}/dist_check.py
            dist_check=$?
            if test $dist_check -eq 0; then
                break
            fi
            set -e
            echo ${sim_dirname}"   Some peptide pairs are too close. Configuring them agein..."
        done


        # Equil the system and execute production run
        mv initial_wat_ion.prmtop initial_wat_ion.inpcrd initial_wat_ion.pdb ../
        cd ../
        jsub -q PN ${PROJ_DIR}/${PROGRAMS_DIRNAME}/production.sh


        echo ${sim_dirname}"   done!"
    } &

    # sleep 1s to avoid making files have same names
    sleep 1s
    cd $PROJ_DIR
done


wait
echo "Now, all preparation operations are successfully finished and main simulations are running!"
exit 0
