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

#include <ROOT/RDataFrame.hxx>
#include <ROOT/RResultPtr.hxx>

#include <chrono> 
using namespace std::chrono; 

using std::cout;
using std::string;
using std::vector;
using std::endl;
using std::pair;
using std::make_pair;
using std::map;

typedef std::map<string, std::pair<string, string>> cuts_t;
typedef std::vector<ROOT::RDF::RResultPtr<TH1D>> vecPtr_t;


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
    
    verbose_ = parser.integerValue("verbose");

    preSelection_ = parser.stringValue("preselection");
    dir_ = parser.stringValue("dir");
    fout_ = new TFile(parser.stringValue("outputfile").c_str(),"RECREATE");
    
  }
  
  
  void makeHiggsShape(string preselection, string categoryselection, string prefix) {
    
    
    pair<float,float> tmp;
    ROOT::EnableImplicitMT();

    cout<<"Make Higgs Shape"<<endl;
    auto start = high_resolution_clock::now();

    std::string fullselection = "("+preselection+")";//*"+weight_;

    /*
    tmp= createHistogramAndShifts(dir_+"ggH120.root","ggH120",(fullselection),luminosity_,prefix);

    tmp= createHistogramAndShifts(dir_+"ggH125.root","ggH125",(fullselection),luminosity_,prefix);
    tmp= createHistogramAndShifts(dir_+"ggH130.root","ggH130",(fullselection),luminosity_,prefix);
    //if(doSys_>0)
    //createShiftsTES("ggH125",dir_+"ggH125.root",categoryselection+"&&"+trigSelection_+"&&"+osSignalSelection_,weight_,luminosity_*legCorr,prefix,tmp);
    
    tmp= createHistogramAndShifts(dir_+"vbfH120.root","qqH120",(fullselection),luminosity_,prefix);
    tmp= createHistogramAndShifts(dir_+"vbfH125.root","qqH125",(fullselection),luminosity_,prefix);
    tmp= createHistogramAndShifts(dir_+"vbfH130.root","qqH130",(fullselection),luminosity_,prefix);
    */
    //if(doSys_>0)
    //createShiftsTES("qqH125",dir_+"vbfH125.root",categoryselection+"&&"+trigSelection_+"&&"+osSignalSelection_,weight_,luminosity_*legCorr,prefix,tmp);
    /*
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
    */
    //std::cout<<"creating met systematics Higgs"<<std::endl;
    //createMETSystematicsHiggs(fullselection, luminosity_, prefix);
    std::cout<<"creating jet systematics Higgs"<<std::endl;
    createJETSystematicsHiggs(fullselection, luminosity_, prefix);

    // Get ending timepoint 
    auto stop = high_resolution_clock::now(); 
  
    // Get duration. Substart timepoints to  
    // get durarion. To cast it to proper unit 
    // use duration cast method 
    auto duration = duration_cast<seconds>(stop - start); 
  
    cout << "Time taken by function: "
	 << duration.count() << " seconds" << endl; 
  }
  
  void close() {
    fout_->Close();
  }

  /**********************************************************************/

  void createJETSystematicsHiggs(string inputSelections, float scale, string prefix) {
    
    createJETSystematicsHiggsForAFile(inputSelections, scale, prefix, dir_+"ggH125.root",  "ggH125_CMS_scale_j_");
    createJETSystematicsHiggsForAFile(inputSelections, scale, prefix, dir_+"vbfH125.root", "qqH125_CMS_scale_j_");
    createJETSystematicsHiggsForAFile(inputSelections, scale, prefix, dir_+"ZH125.root",   "ZH125_CMS_scale_j_");
    createJETSystematicsHiggsForAFile(inputSelections, scale, prefix, dir_+"ttH125.root" , "ttH125_CMS_scale_j_"); 
  }
  
  /**********************************************************************/

  void createJETSystematicsHiggsForAFile(string inputSelections, float scale, string prefix, string filename, string histNamePrefix){
    std::vector<std::string> jetSysVec = {"Closure", "AbsoluteFlavMap", "AbsoluteMPFBias", "AbsoluteScale", "AbsoluteStat", "FlavorQCD", "Fragmentation", "PileUpDataMC", "PileUpPtBB", "PileUpPtEC1", "PileUpPtEC2", "PileUpPtHF", "PileUpPtRef", "RelativeBal", "RelativeFSR", "RelativeJEREC1", "RelativeJEREC2", "RelativeJERHF", "RelativePtBB", "RelativePtEC1", "RelativePtEC2", "RelativePtHF", "RelativeStatEC", "RelativeStatFSR", "RelativeStatHF", "SinglePionECAL", "SinglePionHCAL", "TimePtEta", "Total"};
    
    // Book the cuts in a map
    cuts_t mapCuts;

    for(auto jetSys : jetSysVec){
      
      //      cout << "Making cuts for" << jetSys << endl;

      std::string newSelectionUp=inputSelections;
      std::string newSelectionDown=inputSelections;
      ReplaceStringInPlace(newSelectionUp,   "njets", "njet_"   +jetSys+"Up");
      ReplaceStringInPlace(newSelectionUp,   "mjj"  , "vbfMass_"+jetSys+"Up");
      ReplaceStringInPlace(newSelectionDown, "njets", "njet_"   +jetSys+"Down");
      ReplaceStringInPlace(newSelectionDown, "mjj"  , "vbfMass_"+jetSys+"Down");

      mapCuts.insert(pair<string, pair<string, string>>(filename+"_"+jetSys+"_Up",   make_pair(newSelectionUp,   histNamePrefix+jetSys+"_13TeVUp")));
      mapCuts.insert(pair<string, pair<string, string>>(filename+"_"+jetSys+"_Down", make_pair(newSelectionDown, histNamePrefix+jetSys+"_13TeVDown")));
    }

    ROOT::RDataFrame d((channel_+"EventTree/eventTree").c_str(), filename);

    vecPtr_t results = MakeResultPtrVector(d, mapCuts);

    FinalizeAndWriteHists(results, scale, "", true, false);

  }

  /**********************************************************************/

  void ReplaceStringInPlace(std::string& subject, const std::string& search,
			    const std::string& replace) {
    size_t pos = 0;
    while ((pos = subject.find(search, pos)) != std::string::npos) {
      subject.replace(pos, search.length(), replace);
      pos += replace.length();
    }
  }

  /**********************************************************************/

  /* Applies the cuts in mapCuts to the RDataFrame d, returning a vector
     of the resulting smart pointers. */
  vecPtr_t MakeResultPtrVector(ROOT::RDataFrame d,
							       cuts_t mapCuts){
    
    vecPtr_t vHist;

    cuts_t::iterator itMap = mapCuts.begin();
    
    while (itMap != mapCuts.end())
      {
	// Accessing KEY:
	string key_ = itMap->first;
	  
	// Accessing CUT+KEY pair:
	std::pair<string, string> pair_ = itMap->second;
	  
	// Accessing CUT:
	string cut_ = pair_.first;

	// Accessing NAME:
	string name_ = pair_.second;

	//	cout << key_ << " :: " << cut_ << endl;

	vHist.push_back(d.Filter(cut_).Histo1D({(name_).data(), (name_).data(), bins_, min_, max_}, variable_));

	itMap++;
      }
    return vHist;

  }

  /**********************************************************************/

  /* Loop through a vector of TH1D RResultPtrs and process/write them. */
  void FinalizeAndWriteHists(vecPtr_t vHist, 
			     float scaleFactor = 1, string postfix = "", bool normUC  = true, bool keys=false){

    ROOT::EnableImplicitMT();

    string folder = filelabel_+postfix;
    if(fout_->Get(folder.c_str())==0)
      fout_->mkdir(folder.c_str());
    fout_->cd(folder.c_str());

    for (auto it = vHist.begin(); it != vHist.end(); ++it) {
      TH1D* h = it->GetPtr();
      h->Sumw2();
      h->Scale(scaleFactor);
      h->Write(h->GetName(), TObject::kOverwrite);

      cout << h->GetName() << ": " << h->Integral() << ", entries: " << h->GetEntries() << endl;
    }
  }

  /**********************************************************************/

  /* string name: name of the variable to draw. */

  pair<float,float> createHistogramAndShifts(string file, string name, string nominalCut, 
					     float scaleFactor = 1, string postfix = "",bool normUC  = true, bool keys=false) {
    ROOT::EnableImplicitMT();

    string folder = filelabel_+postfix;
    if(fout_->Get(folder.c_str())==0)
      fout_->mkdir(folder.c_str());
    fout_->cd(folder.c_str());

    //    std::cout << file << std::endl;

    ROOT::RDataFrame d((channel_+"EventTree/eventTree").c_str(), file);

    ////////////////////////////////////
    // Nominal histogram
    ////////////////////////////////////
    
    auto hNominal = d.Filter(nominalCut).Histo1D({name.data(), name.data(), bins_, min_, max_}, "m_sv");
    
    
    ////////////////////////////////////
    // Shifts TES
    ////////////////////////////////////
    // Do TES shifts if it's a H125.root file
    bool doTES = (file.find("125") != std::string::npos);
    
    if (doTES){
      std::cout << "Found a 125 ROOT file, doing TES shifts" << std::endl;

      cuts_t mapCuts;

      mapCuts.insert(pair<string, pair<string, string>>("ptSelectionDM0Up_",   make_pair("(((pt_2*0.988)>40&&decayMode_2==0)||((pt_2*1.010)>40&&decayMode_2==1)||((pt_2*1.004)>40&&decayMode_2==10))&&(((pt_1*0.988)>50&&decayMode_1==0)||((pt_1*1.010)>50&&decayMode_1==1)||((pt_1*1.004)>50&&decayMode_1==10))", "_CMS_scale_t_1prong_13TeVUp")));
      mapCuts.insert(pair<string, pair<string, string>>("ptSelectionDM0Down_", make_pair("(((pt_2*0.976)>40&&decayMode_2==0)||((pt_2*1.010)>40&&decayMode_2==1)||((pt_2*1.004)>40&&decayMode_2==10))&&(((pt_1*0.976)>50&&decayMode_1==0)||((pt_1*1.010)>50&&decayMode_1==1)||((pt_1*1.004)>50&&decayMode_1==10))", "_CMS_scale_t_1prong_13TeVDown")));  

      mapCuts.insert(pair<string, pair<string, string>>("ptSelectionDM1Up_",   make_pair("(((pt_2*0.982)>40&&decayMode_2==0)||((pt_2*1.016)>40&&decayMode_2==1)||((pt_2*1.004)>40&&decayMode_2==10))&&(((pt_1*0.982)>50&&decayMode_1==0)||((pt_1*1.016)>50&&decayMode_1==1)||((pt_1*1.004)>50&&decayMode_1==10))", "_CMS_scale_t_1prong1pizero_13TeVUp")));
      mapCuts.insert(pair<string, pair<string, string>>("ptSelectionDM1Down_", make_pair("(((pt_2*0.982)>40&&decayMode_2==0)||((pt_2*1.004)>40&&decayMode_2==1)||((pt_2*1.004)>40&&decayMode_2==10))&&(((pt_1*0.982)>50&&decayMode_1==0)||((pt_1*1.004)>50&&decayMode_1==1)||((pt_1*1.004)>50&&decayMode_1==10))", "_CMS_scale_t_1prong1pizero_13TeVDown")));
      
      mapCuts.insert(pair<string, pair<string, string>>("ptSelectionDM10Up_",   make_pair("(((pt_2*0.982)>40&&decayMode_2==0)||((pt_2*1.010)>40&&decayMode_2==1)||((pt_2*1.010)>40&&decayMode_2==10))&&(((pt_1*0.982)>50&&decayMode_1==0)||((pt_1*1.010)>50&&decayMode_1==1)||((pt_1*1.010)>50&&decayMode_1==10))", "_CMS_scale_t_3prong_13TeVUp")));
      mapCuts.insert(pair<string, pair<string, string>>("ptSelectionDM10Down_", make_pair("(((pt_2*0.982)>40&&decayMode_2==0)||((pt_2*1.010)>40&&decayMode_2==1)||((pt_2*0.998)>40&&decayMode_2==10))&&(((pt_1*0.982)>50&&decayMode_1==0)||((pt_1*1.010)>50&&decayMode_1==1)||((pt_1*0.998)>50&&decayMode_1==10))", "_CMS_scale_t_3prong_13TeVDown")));

      vecPtr_t vHist;
      
      // Fill the vector of RResultPtr by looping through the map
      cuts_t::iterator itMap = mapCuts.begin();

      while (itMap != mapCuts.end())
	{
	  // Accessing CUTNAME:
	  string cutname_ = itMap->first;
	  
	  // Accessing CUT+KEY pair:
	  std::pair<string, string> pair_ = itMap->second;
	  
	  // Accessing CUT:
	  string cut_ = pair_.first;

	  // Accessing KEY:
	  string key_ = pair_.second;

	  //	  cout << key_ << " :: " << cut_ << endl;

	  vHist.push_back(d.Filter(nominalCut+"&&"+cut_).Histo1D({(name+key_).data(), (name+key_).data(), bins_, min_, max_}, variable_));

	  itMap++;
	}
      
      for (auto it = vHist.begin(); it != vHist.end(); ++it) {
	TH1D* h = it->GetPtr();
	h->Sumw2();
	h->Scale(scaleFactor);
	h->Write(h->GetName(), TObject::kOverwrite);

      }
    }
    // Read & write nominal histogram
    hNominal->Sumw2();
    hNominal->Scale(scaleFactor);
    hNominal->Write(name.data(), TObject::kOverwrite);
    
    Double_t error=0.0;
    
    float yield = hNominal->IntegralAndError(1,hNominal->GetNbinsX(),error,"");

    if(yield == 0){
      hNominal->SetBinContent(1,0.00001);
    }
			  
    return make_pair(yield,error);


  }

  /**********************************************************************/

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
