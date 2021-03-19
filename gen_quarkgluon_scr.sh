#!/usr/bin/env bash

starting_dir=$(pwd)
cards_dir="$(pwd)/cards"
gridpack_dir="$(pwd)/gridpacks"


# ----- gluon jets ----- #
jet_type="gluon"
gridpack="${gridpack_dir}/${jet_type}_500GeV_gridpack.tar.gz"
pythia_card="${cards_dir}/configLHE.cmnd"
delphes_card="${cards_dir}/delphes_card_ATLAS_uniformGranularity.tcl"
nbatch_50k=80
init_seed=1000
event_tag="gluons_500GeV_march162021"

python ./gen_run_scripts.py \
  -g $gridpack \
  --pythia_card $pythia_card \
  --delphes_card $delphes_card \
  --nbatch_50k $nbatch_50k \
  --init_seed $init_seed \
  --event_tag $event_tag


# ----- quark jets ----- #
jet_type="quark"
gridpack="${gridpack_dir}/${jet_type}_500GeV_gridpack.tar.gz"
pythia_card="${cards_dir}/configLHE.cmnd"
delphes_card="${cards_dir}/delphes_card_ATLAS_uniformGranularity.tcl"
nbatch_50k=80
init_seed=2000
event_tag="quarks_500GeV_march162021"

python ./gen_run_scripts.py \
  -g $gridpack \
  --pythia_card $pythia_card \
  --delphes_card $delphes_card \
  --nbatch_50k $nbatch_50k \
  --init_seed $init_seed \
  --event_tag $event_tag
