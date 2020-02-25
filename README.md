# Data Card Creator (minimal example from PUAnalysis Framework)

## Setup instructions

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
To run the minimal example of the original Data Card Creator:
```
cd StatTools/data
bash makePlot
```

**(In-development!)** To run the minimal example of the RDataFrame creator:
```
cd StatTools/data
bash makePlot_RDF
```