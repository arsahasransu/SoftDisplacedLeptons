The files here create varibale plots from hepmc file. They have been used to compare the kinematic distributions between madgraph
samples generated at different locations. 

Pre-requisites:

The HepMC reader package needs to be installed from http://lcgapp.cern.ch/project/simu/HepMC/.
The following variables need to be set.

```
export HEPMC_DIR=<installation-dir>
export LD_LIBRARY_PATH=$HEPMC_DIR/lib:$LD_LIBRARY_PATH
```

To Run:
```
make example_CompareHepmc.exe
```
