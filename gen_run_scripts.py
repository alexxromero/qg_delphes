"""
gen_event_scripts.py generates the scripts used for unpacking the MG jet
gridpack, and subsequent Pythia and Delphes implementations.
"""

import sys
import os
import re
import subprocess as sp
from argparse import ArgumentParser
from textwrap import dedent

def gen_event_specs_scr():
    script = dedent(
        """\
        #!/usr/bin/env bash

        gridpack=$1
        pythia_card=$2
        delphes_card=$3
        initial_seed=$4
        event_tag=$5

        echo "----- run specs -----"
        echo "start: $(date)"
        echo "input arguments:"
        echo "gridpack  ::  ${gridpack}"
        echo "pythia card ::  ${pythia_card}"
        echo "delphes card ::  ${delphes_card}"
        echo "initial_seed ::  ${initial_seed}"
        echo "tag ::  ${tag}"

        echo "Setting up ROOT..."
        export ATLAS_LOCAL_ROOT_BASE=/cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase
        source ${ATLAS_LOCAL_ROOT_BASE}/user/atlasLocalSetup.sh
        lsetup "root 6.14.04-x86_64-slc6-gcc62-opt"

        output_dir="/DFS-L/DATA/atlas/${USER}/qg_eflow/${event_tag}"
        mkdir -p $output_dir
        mkdir -p ${output_dir}/${tag}
        mkdir -p ${output_dir}/${tag}/LHE      # store LHE events
        mkdir -p ${output_dir}/${tag}/Delphes  # store Delphes ROOT files
        mkdir -p ${output_dir}/${tag}/ProcessedJets  # store jets after
                                                     # pT- and eta-filtering

        temp_dir="/DFS-L/SCRATCH/atlas/${USER}/${event_tag}_${initial_seed}"
        mkdir -p $temp_dir
        cd $temp_dir

        echo "Setting up Pythia8..."
        wget http://home.thep.lu.se/~torbjorn/pythia8/pythia8235.tgz
        tar -xzvf pythia8235.tgz
        cd pythia8235
        ./configure --prefix=$(pwd)
        make install
        export PYTHIA8=$(pwd)  # variable needed if using Delphes
        cd ..

        echo "Setting up Delphes..."
        wget http://cp3.irmp.ucl.ac.be/downloads/Delphes-3.4.2.tar.gz
        tar -xzvf Delphes-3.4.2.tar.gz
        cd Delphes-3.4.2
        export DELPHESDIR=$(pwd)
        make HAS_PYTHIA8=true
        cd ..
        cp ${HOME}/qg_delphes/delphes_scripts/analyze_eflowJets.cpp ./Delphes-3.4.2
        cp ${HOME}/qg_delphes/delphes_scripts/analyze_eflowJets.h ./Delphes-3.4.2

        echo "Unpacking the gridpack..."
        temp_gridpack="${temp_dir}/gridpack_${event_tag}_${initial_seed}.tar.gz"
        cp $gridpack $temp_gridpack
        tar -zxvf $temp_gridpack

        echo "Compiling MG5..."
        cd madevent
        ./bin/compile
        cd ..
        chmod +x run.sh

        echo "Running MG5..."
        nEvents_run=10000  # it's recommended to run batches of no more than 10k
        seed=$initial_seed
        for i in {0..4}
        do
            ./run.sh $nEvents_run $seed  # generates events.lhe.gz
            gunzip ./events.lhe.gz
            echo "Running Pythia8+Delphes for seed $seed..."
            ${DELPHESDIR}/DelphesPythia8 $delphes_card $pythia_card delphes.root

            # process the jets
            cd Delphes-3.4.2
            root -l -b -q analyze_eflowJets.cpp'("../delphes.root", "processed_jets.root")'
            cd ..

            mv ./events.lhe ${output_dir}/${tag}/LHE/events_${seed}.lhe
            mv ./delphes.root ${output_dir}/${tag}/Delphes/delphes_${seed}.root
            mv ./Delphes-3.4.2/processed_jets.root ${output_dir}/${tag}/ProcessedJets/processed_jets_${seed}.root
            echo "done with seed $seed"
            ((seed++))
        done

        echo "DONE. :)"

        cd $HOME
        rm -rf $temp_dir
        """
    )
    return script


def gen_submit_jobs_scr(event_specs_scr, gridpack, pythia_card, delphes_card,
                        nbatch_50k, init_seed, event_tag):
    script = dedent(
        """\
        #!/usr/bin/env bash
        set -euo pipefail
        mkdir -p logs
        """
    )

    # each batch calls DelphesPythia8 10 times with 10k events per run
    seed = init_seed
    for i in range(nbatch_50k):
        script += dedent(
            """\

            FLGS='-t 300 -p atlas -c 2 -o logs/out-%j.txt -e logs/error-%j.txt'
            sbatch ${{FLGS}} {} {} {} {} {} {}
            """.format(event_specs_scr,
                       gridpack, pythia_card, delphes_card, seed, event_tag)
        )
        seed += 5
    return script


def generate_scripts(gridpack, pythia_card, delphes_card,
                     nbatch_50k, init_seed, event_tag,
                     event_dir):
    event_specs_scr = os.path.join(event_dir, "MG5PythiaDelphes_specs.sh")
    submit_jobs_scr = os.path.join(event_dir, "submit_quarkgluon.sh")

    with open(event_specs_scr, "w") as f:
        f.write(gen_event_specs_scr())

    with open(submit_jobs_scr, "w") as f:
        f.write(gen_submit_jobs_scr(event_specs_scr, gridpack, pythia_card,
                                    delphes_card, nbatch_50k, init_seed,
                                    event_tag))

if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument("-g", "--gridpack", type=str,  required=True)
    parser.add_argument("--pythia_card", type=str,  required=True)
    parser.add_argument("--delphes_card", type=str,  required=True)
    parser.add_argument("--nbatch_50k", type=int,  default=1)
    parser.add_argument("--init_seed", type=int,  default=1)
    parser.add_argument("--event_tag", type=str,  default="dummy_event")
    args = parser.parse_args()

    start_dir = os.getcwd()
    event_dir = os.path.join(start_dir, args.event_tag)
    if not os.path.exists(event_dir):
        os.makedirs(event_dir)
    os.chdir(event_dir)
    generate_scripts(args.gridpack, args.pythia_card, args.delphes_card,
                     args.nbatch_50k, args.init_seed, args.event_tag,
                     event_dir)
    os.chdir(start_dir)
