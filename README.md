# Data Card Creator (minimal example from PUAnalysis Framework)

## Setup instructions (do only once)

```
scram project CMSSW CMSSW_10_6_10
cmsenv
readlink -f $(which root)
```
The last line from the output should show that the environment
for ROOT 6.14.09 has been set up. Support for [RDataFrame](https://root.cern/doc/master/classROOT_1_1RDataFrame.html)
was only added in ROOT 6.14. Notably, ROOT 6.14 does not have
the libROOTDataFrame.so library automatically linked in the
default libraries (it is in ROOT 6.18, at least), so the
BuildFiles in the StatTools folder link the correct library.

## Usage
Every time the Data Card Creator is changed, re-compile from the top directory (or from StatTools if
something in StatTools was changed):
```
scram b USER_CXXFLAGS="-Wno-error=unused-but-set-variable"
```

### Standard TTree->Draw() version
To run the minimal example of the original Data Card Creator:
```
cd StatTools/data
bash makePlot 
```

### RDataFrame example
```
cd StatTools/data
bash makePlot_RDF
```

If the verbose option is on/ something is written to std::cout, it's a good
idea to pipe the std::cout to a text file to check if the histogram yields and number of entries remain unchanged as one iterates through edits:
```  
sort out.txt > out.sorted
sort outRDF.txt > outRDF.sorted
diff -y out.sorted outRDF.sorted
```

### Boost example (in-development)
If using a virtual environment on lxplus, just `pip install boost-histogram`. If you are not using environments,
you can do a user install (`--user`), which is not as good as an environmental in general, but works. 
```
cd StatTools/data
bash makePlot_Boost
```