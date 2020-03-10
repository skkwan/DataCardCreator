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
typedef std::map<string, cuts_t>                     dictionary_t;

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
    trigSelection_        = parser.stringValue("trigSelection");
    trigSelectionData_    = parser.stringValue("trigSelectionData");

    dir_ = parser.stringValue("dir");
    fout_ = new TFile(parser.stringValue("outputfile").c_str(),"RECREATE");
    
    //read systematic uncertainties 
    luminosity_    = parser.doubleValue("luminosity");
    
  }
  
  
  void makeHiggsShape(string preselection, string categoryselection, string prefix) {
    
    
    pair<float,float> tmp;
    ROOT::EnableImplicitMT();

    cout<<"Make Higgs Shape"<<endl;
    auto start = high_resolution_clock::now();

    std::string fullselection = "("+preselection+")";//*"+weight_;

    // SYNTAX:      BookCut(string filename, string cutName, string cut, string histName)
    //BookCut(dir_+"ggH120.root", "ggH120_nominal", fullselection, "ggH120");

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
    
    std::cout<<"creating jet systematics Higgs with luminosity " << luminosity_ << endl;
    createJETSystematicsHiggs(fullselection, luminosity_, prefix);

    FinishUp(luminosity_, prefix, false, false);

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

    // TEMPORARY: Testing making 108 histograms per file
    std::string trigselection = "("+trigSelection_+")";
    
    createJETSystematicsHiggsForAFile(trigselection, scale, prefix, dir_+"ggH125.root",  "dummy_ggH125_CMS_scale_j_");

    createJETSystematicsHiggsForAFile(trigselection, scale, prefix, dir_+"vbfH125.root", "dummy_qqH125_CMS_scale_j_");

    createJETSystematicsHiggsForAFile(trigselection, scale, prefix, dir_+"ZH125.root",   "dummy_ZH125_CMS_scale_j_");
    createJETSystematicsHiggsForAFile(trigselection, scale, prefix, dir_+"ttH125.root" , "dummy_ttH125_CMS_scale_j_");
  

  // TEMPORARY: Testing making 162 histograms per file
  std::string bothselection = "("+trigSelection_+"&&"+inputSelections+")";
  createJETSystematicsHiggsForAFile(bothselection,  scale, prefix, dir_+"ggH125.root",  "dummy2_ggH125_CMS_scale_j_");

  createJETSystematicsHiggsForAFile(bothselection, scale, prefix, dir_+"vbfH125.root", "dummy2_qqH125_CMS_scale_j_");

  createJETSystematicsHiggsForAFile(bothselection, scale, prefix, dir_+"ZH125.root",   "dummy2_ZH125_CMS_scale_j_");
  createJETSystematicsHiggsForAFile(bothselection, scale, prefix, dir_+"ttH125.root" , "dummy2_ttH125_CMS_scale_j_");
  }
  /**********************************************************************/

  void createJETSystematicsHiggsForAFile(string inputSelections, float scale, string prefix, string filename, string histNamePrefix){
    std::vector<std::string> jetSysVec = {"Closure","AbsoluteFlavMap", "AbsoluteMPFBias", "AbsoluteScale", "AbsoluteStat", "FlavorQCD", "Fragmentation", "PileUpDataMC", "PileUpPtBB", "PileUpPtEC1", "PileUpPtEC2", "PileUpPtHF", "PileUpPtRef", "RelativeBal", "RelativeFSR", "RelativeJEREC1", "RelativeJEREC2", "RelativeJERHF", "RelativePtBB", "RelativePtEC1", "RelativePtEC2", "RelativePtHF", "RelativeStatEC", "RelativeStatFSR", "RelativeStatHF", "SinglePionECAL", "SinglePionHCAL", "TimePtEta", "Total"};
    
    // Book the cuts in a map
    cuts_t mapCuts;

    for(auto jetSys : jetSysVec){
      std::string newSelectionUp=inputSelections;
      std::string newSelectionDown=inputSelections;
      ReplaceStringInPlace(newSelectionUp,   "njets", "njet_"   +jetSys+"Up");
      ReplaceStringInPlace(newSelectionUp,   "mjj"  , "vbfMass_"+jetSys+"Up");
      ReplaceStringInPlace(newSelectionDown, "njets", "njet_"   +jetSys+"Down");
      ReplaceStringInPlace(newSelectionDown, "mjj"  , "vbfMass_"+jetSys+"Down");

      // SYNTAX:      BookCut(string filename, string cutName, string cut, string histName)
      BookCut(filename, prefix+filename+"_"+jetSys+"_Up",  newSelectionUp,   histNamePrefix+jetSys+"_13TeVUp");
      BookCut(filename, prefix+filename+"_"+jetSys+"_Down",newSelectionDown, histNamePrefix+jetSys+"_13TeVDown");
 
    }
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

  /* "filename" : name of ROOT file. "cutName": must be unique for each 
     cut. "cut": the actual selections. "histName": the TH1D name. */

  void BookCut(string filename, string cutName, string cut, string histName){

    // Check if an entry for that filename already exists in dict_
    dictionary_t::iterator itDict = dict_.find(filename);
    if (itDict != dict_.end())
      {
	// Element filename already exists
	cuts_t::iterator itCut;
	itCut = (itDict->second).find(cutName);
	if (itCut != (itDict->second).end())
	  cout<< "[WARNING: trying to add cut with existing name -- skipped]" << endl;
	else
	  (itDict->second).insert({cutName+histName, make_pair(cut, histName)});
	//	cout << filename << " is already in dict_, inserting " << histName << endl;
      }
    else
      {
	// Element filename does not exist yet. Make a new cuts_t object.
	cuts_t newCuts;
	newCuts.insert({cutName+histName, make_pair(cut, histName)});
	
	dict_.insert({filename, newCuts});
	//	cout << "Adding " << filename << " to dict, inserting " << histName << endl;
      }

  }


  /**********************************************************************/

  /* Applies the cuts in mapCuts to the RDataFrame d, returning a vector
     of the resulting smart pointers. */
  vecPtr_t MakeResultPtrVector(ROOT::RDataFrame d, cuts_t mapCuts){
    
    vecPtr_t vHist;

    cout << "Making result pointer vector: ... " << endl;

    for (cuts_t::iterator itMap = mapCuts.begin(); itMap != mapCuts.end(); ++itMap)
      {
	// Accessing the key (unique name of the cut)
	// string key_ = itMap->first;
	  
	// Accessing the value (a pair consisting of the cutstring and the desired hist name)
	std::pair<string, string> pair_ = itMap->second;
	string cut_ = pair_.first;
	string name_ = pair_.second;

	// cout << key_ << " :: " << name_ << endl;
	cout << "name :: " << name_ << endl;

	vHist.push_back(d.Filter(cut_).Histo1D({(name_).data(), (name_).data(), bins_, min_, max_}, variable_));
      }

    cout << "-------------- Done with making result pointer ------- " << endl;
    return vHist;

  }

  /**********************************************************************/

  /* Loop through a vector of TH1D RResultPtrs and process/write them. */

  void FinalizeAndWriteHists(vecPtr_t vHist, 
			     float scaleFactor = 1, string postfix = "", bool normUC  = true, bool keys=false){

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

  /*  Accesses the private class member dict_ and executes all the cuts
      booked in it. */

  void FinishUp(float scaleFactor = 1, string postfix = "", bool normUC  = true, bool keys=false){
   
    cout << "Now in FinishUp: ... " <<endl;

    for (dictionary_t::iterator itDict = dict_.begin(); itDict != dict_.end(); ++itDict)
      {
	// Get the filename
	string filename = itDict->first;
	
	// Get the map of cuts
	cuts_t allCuts = itDict->second;

	// Make RDataFrame
	ROOT::RDataFrame d((channel_+"EventTree/eventTree").c_str(), filename);
	
	// Make the vector of smart pointers 
	vecPtr_t vecHists = MakeResultPtrVector(d, allCuts);
	
	// Execute the smart pointers (slow step)
	FinalizeAndWriteHists(vecHists, scaleFactor, postfix, normUC, keys);
	
      }
  }

  /**********************************************************************/

 private:		
  string channel_;
  string filelabel_;
  
  //files
  TFile *fout_;
  int verbose_;

  string preSelection_;
  string trigSelection_;
  string trigSelectionData_;
  
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
  
  dictionary_t dict_;
  
};
