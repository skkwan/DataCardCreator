/******************************************************/
/* DataCardCreatorHThTh_2016_RDF.h                    */
/******************************************************/

#include "TCanvas.h"
#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TH1.h"
#include "TF1.h"
#include "PhysicsTools/FWLite/interface/CommandLineParser.h" 
#include <math.h>
#include <boost/algorithm/string.hpp>
#include "DataCardCreator/DataFormats/interface/TH1Keys.h"
#include <TEfficiency.h>
#include <TMath.h>
#include <iostream>
#include <string>
#include <boost/lexical_cast.hpp>
//#include "FineBins.h"

/******************************************************/

pair<float,float> createHistogramAndShifts(string file,string name, string cut,float scaleFactor = 1, string postfix = "",bool normUC  = true, bool keys=false) {
    
  TFile *f  = new TFile(file.c_str());
  if(f==0) printf("Not file Found\n");
    
  TTree *t= (TTree*)f->Get((channel_+"EventTree/eventTree").c_str());
    
  if(t==0) printf("Not Tree Found in file %s\n",file.c_str());
  pair<float,float> yield;
    
  //if(!keys)
  yield =makeHistogram(t,filelabel_+postfix,name,cut,scaleFactor);
  //else
  //yield =makeTHKeys(t,filelabel_+postfix,name,cut,scaleFactor);
    
  f->Close();
  return yield;
}


/******************************************************/

