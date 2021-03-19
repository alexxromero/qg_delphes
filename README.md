# MG5+Pythia8+Delphes Jet Generation

This package is used to generate light-quark and gluon jets in UCI's green planet cluster. The jets are generated via MadGraph5, followed by shower and hadronization on Pythia8, and Delphes detector response simulation.

Run parameters -- specified in gen_quarkgluon_scr.sh:

gridpack :: MG gridpack -- see gridpacks for how to generate the grispacks
delphes_card :: Delphes card
pythia_card :: Pythia card
nbatch_50k :: No. of 50k jet event batches to generate
init_seed :: Initial seed passed to the MG gridpack. Increases per batch
event_tag :: User-set tag of the event


Usage:
From within this directory, event:

    bash gen_quarkgluon_scr.sh  # generates event scripts
    bash [event_tag]/submit_quarkgluon.sh   # submit jobs
