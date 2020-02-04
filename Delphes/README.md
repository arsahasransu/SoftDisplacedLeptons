# myCMS.tcl 
Customized detector card to run Delphes on.

# To Run Delphes
./DelphesHepMC ./cards/my_CMS.tcl Delp3_4_1_ppTobb_Cuts2_BatchCTR.root MG5_PY8_ppTobb_NLO_Cuts2_BatchCTR.hepmc

# delphesObjectSelector.C 

Run this on delphes produced root files to get ntuples with selected objects. 
Makes sure to give appropriate file input and output filenames.

Syntax:
root delphesObjectSelector.C
