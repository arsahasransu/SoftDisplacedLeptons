# plotVarDisb_Objects_*.C
Plots the variable histogram in form of a pdf.
Every file of this form plots variables of different category. The grouping is similar to that in the overleaf report.

# VariableMaker.C
Writes the interesting variables to a ROOT file. This ROOT file can be used by the classifier network to train upon.

How to compile:
```
root 
.L VariableMaker.C+ // compile, after that you can execute any of the functions from the command line
```

How to create one ROOT file (takes wild cards):

```
root 
.L VariableMaker.C+ // compile, after that you can execute any of the functions from the command line
runononesample("../../<pathtofiles>_*.root","outputfile.root")
```
