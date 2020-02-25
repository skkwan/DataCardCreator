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

//#include "ROOT/TDataFrame.hxx"
//using namespace ROOT::Experimental;

#include <chrono> 
using namespace std::chrono; 

using std::cout;
using std::string;
using std::vector;
using std::endl;
using std::pair;
using std::make_pair;

class DataCardCreatorHThTh_2016_RDF {
	public:
  
  DataCardCreatorHThTh_2016_RDF(optutl::CommandLineParser parser) {
    channel_  = parser.stringValue("channel");
    
    string name_=channel_;
    filelabel_="";
    
    //define the histogram binning
    bins_ = parser.integerValue("bins");
    min_  = parser.doubleValue("min");
    max_  = parser.doubleValue("max");
    variable_      = parser.stringValue("variable");    
    
    preSelection_ = parser.stringValue("preselection");
    dir_ = parser.stringValue("dir");
    fout_ = new TFile(parser.stringValue("outputfile").c_str(),"RECREATE");
    
  }
  
  
  void makeHiggsShape(string preselection, string categoryselection, string prefix) {
    
    pair<float,float> tmp;
    
    cout<<"Make Higgs Shape"<<endl;
    auto start = high_resolution_clock::now(); 

    std::string fullselection = "("+preselection+")";//*"+weight_;
    tmp= createHistogramAndShifts(dir_+"ggH120.root","ggH120",(fullselection),luminosity_,prefix);

    tmp= createHistogramAndShifts(dir_+"ggH125.root","ggH125",(fullselection),luminosity_,prefix);
    tmp= createHistogramAndShifts(dir_+"ggH130.root","ggH130",(fullselection),luminosity_,prefix);
    //if(doSys_>0)
    //createShiftsTES("ggH125",dir_+"ggH125.root",categoryselection+"&&"+trigSelection_+"&&"+osSignalSelection_,weight_,luminosity_*legCorr,prefix,tmp);
    
    tmp= createHistogramAndShifts(dir_+"vbfH120.root","qqH120",(fullselection),luminosity_,prefix);
    tmp= createHistogramAndShifts(dir_+"vbfH125.root","qqH125",(fullselection),luminosity_,prefix);
    tmp= createHistogramAndShifts(dir_+"vbfH130.root","qqH130",(fullselection),luminosity_,prefix);
    
    //if(doSys_>0)
    //createShiftsTES("qqH125",dir_+"vbfH125.root",categoryselection+"&&"+trigSelection_+"&&"+osSignalSelection_,weight_,luminosity_*legCorr,prefix,tmp);
    
    tmp= createHistogramAndShifts(dir_+"ZH120.root","ZH120",(fullselection),luminosity_,prefix);
    tmp= createHistogramAndShifts(dir_+"ZH125.root","ZH125",(fullselection),luminosity_,prefix);
    tmp= createHistogramAndShifts(dir_+"ZH130.root","ZH130",(fullselection),luminosity_,prefix);
    
    //if(doSys_>0)
    //createShiftsTES("ZH125",dir_+"ZH125.root",categoryselection+"&&"+trigSelection_+"&&"+osSignalSelection_,weight_,luminosity_,prefix,tmp);
    
    tmp= createHistogramAndShifts(dir_+"WpH120.root","WpH120",(fullselection),luminosity_,prefix);
    tmp= createHistogramAndShifts(dir_+"WpH125.root","WpH125",(fullselection),luminosity_,prefix);
    tmp= createHistogramAndShifts(dir_+"WpH130.root","WpH130",(fullselection),luminosity_,prefix);
    
    tmp= createHistogramAndShifts(dir_+"WmH120.root","WmH120",(fullselection),luminosity_,prefix);
    tmp= createHistogramAndShifts(dir_+"WmH125.root","WmH125",(fullselection),luminosity_,prefix);
    tmp= createHistogramAndShifts(dir_+"WmH130.root","WmH130",(fullselection),luminosity_,prefix);
    
    //createShiftsTES("WH125",dir_+"WH125.root",categoryselection+"&&"+trigSelection_+"&&"+osSignalSelection_,weight_,luminosity_,prefix,tmp);
    
    tmp= createHistogramAndShifts(dir_+"ttH120.root","ttH120",(fullselection),luminosity_,prefix);
    tmp= createHistogramAndShifts(dir_+"ttH125.root","ttH125",(fullselection),luminosity_,prefix);
    tmp= createHistogramAndShifts(dir_+"ttH130.root","ttH130",(fullselection),luminosity_,prefix);
    //std::cout<<"creating met systematics Higgs"<<std::endl;
    //createMETSystematicsHiggs(fullselection, luminosity_, prefix);
    //std::cout<<"creating jet systematics Higgs"<<std::endl;
    //createJETSystematicsHiggs(fullselection, luminosity_, prefix);

    // Get ending timepoint 
    auto stop = high_resolution_clock::now(); 
  
    // Get duration. Substart timepoints to  
    // get durarion. To cast it to proper unit 
    // use duration cast method 
    auto duration = duration_cast<microseconds>(stop - start); 
  
    cout << "Time taken by function: "
         << duration.count() << " microseconds" << endl; 
  }
  
  void close() {
    fout_->Close();
  }
  
  pair<float,float> createHistogramAndShifts(string file,string name, string cut,float scaleFactor = 1, string postfix = "",bool normUC  = true, bool keys=false) {
    
    TFile *f  = new TFile(file.c_str());
    if(f==0) printf("Not file Found\n");
    
    TTree *t= (TTree*)f->Get((channel_+"EventTree/eventTree").c_str());

    ROOT::Experimental::TDataFrame d(file, (channel_+"EventTree/eventTree").c_str(), {"m_sv"});
    // Example of a simple cut
    auto cut = d.Filter("m_sv > 0");
    auto nEntries = cut.Count();
    std::cout << *nEntries << " passed all filters" << std::endl;
    
    if(t==0) printf("Not Tree Found in file %s\n",file.c_str());
    pair<float,float> yield;
    
    //if(!keys)
    yield =makeHistogram(t,filelabel_+postfix,name,cut,scaleFactor);
      //else
      //yield =makeTHKeys(t,filelabel_+postfix,name,cut,scaleFactor);
    
    f->Close();
    return yield;
  }
  
  
  pair<float,float> makeHistogram(TTree* tree,string folder,TString name,string cut,float scaleFactor = 1.) {

    if(fout_->Get(folder.c_str())==0)
      fout_->mkdir(folder.c_str());
    
    TH1F *h=0;
    h= new TH1F(name,name,bins_,min_,max_);
    tree->Draw(variable_+">>"+name,cut.c_str());
    h->Sumw2();
    //h->Scale(scaleFactor);
    fout_->cd(folder.c_str());

    Double_t error=0.0;
    
    float yield = h->IntegralAndError(1,h->GetNbinsX(),error,"");
    if(yield == 0){
      h->SetBinContent(1,0.00001);
    }
    
    h->Write(h->GetName(),TObject::kOverwrite);

    return make_pair(yield,error);
  }
  
  
 private:		
  string channel_;
  string filelabel_;
  
  //files
  TFile *fout_;
  int verbose_;
  string preSelection_;
  
  //Luminosity and efficiency corrections
  float luminosity_;
  float luminosityErr_;
  
  //histogram options
  TString variable_;
  int bins_;
  float min_;
  float max_;
  
  string weight_;
  string dir_;
  
  
  
};