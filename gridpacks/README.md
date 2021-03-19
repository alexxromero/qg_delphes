# MG5 GRIDPACK generation

A gridpack is a tarball file that contains the calculations to run a specified MG5 event. To instructions to generate a gridpack for light-quark and gluon
jets are shown  below.

No need to setup Pythia8/Delphes at this stage.    
The parameters below modify the run_card.dat to generate pp collisions at 14 TeV, and the jets to have pT in the range [495, 444] GeV.  

From the MG5 directory:

### Gluon Jets

    ./bin/mg5_aMC   
    MG5_aMC>generate p p > g g   
    MG5_aMC>output gluon_jets_500GeV  
    MG5_aMC>launch gluon_jets_500GeV/  
    >  
    >set lpp1 1  
    >set lpp2 1  
    >set ebeam1 7000.0  
    >set ebeam2 7000.0  
    >set gridpack True  
    >set ptj 495  
    >set ptjmax 555  
    >  

### Quark Jets

    ./bin/mg5_aMC   
    MG5_aMC>define q = u d s
    MG5_aMC>define q~ = u~ d~ s~
    MG5_aMC>generate p p > q q
    MG5_aMC>add process p p > q q~   
    MG5_aMC>output quark_jets_500GeV  
    MG5_aMC>launch quark_jets_500GeV/  
    >  
    >set lpp1 1  
    >set lpp2 1  
    >set ebeam1 7000.0  
    >set ebeam2 7000.0  
    >set gridpack True  
    >set ptj 495  
    >set ptjmax 555  
    >  

This will generate the files run_01_gridpack.tar.gz in the gluon_jets_500GeV/ and quark_jets_500GeV/ directories. 
