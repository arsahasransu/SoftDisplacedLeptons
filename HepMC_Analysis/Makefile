#LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$HEPMC_DIR/lib

HepMCdir = $HEPMC_DIR
HepMCinc = $(HepMCdir)/include
HepMClib = $(HepMCdir)/lib -lHepMC
ROOT_Config = `root-config --cflags --libs`

CXX = g++

all: example_EventSelection.exe example_StudyHepMC.exe example_CompareHepmc.exe

example_EventSelection.exe:
	$(CXX) -o out_EventSelection.exe example_EventSelection.cc $(ROOT_Config) -I$(HEPMC_DIR)/include -L/$(HEPMC_DIR)/lib -lHepMC

example_StudyHepMC.exe:
	$(CXX) -o out_StudyHepMC.exe example_StudyHepMC.cc $(ROOT_Config) -I$(HEPMC_DIR)/include -L/$(HEPMC_DIR)/lib -lHepMC

example_CompareHepmc.exe:
	$(CXX) -o out_CompareHepmc.exe example_CompareHepmc.cc $(ROOT_Config) -I$(HEPMC_DIR)/include -L/$(HEPMC_DIR)/lib -lHepMC -ldl
