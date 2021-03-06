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


#include <chrono> 
using namespace std::chrono; 

using std::cout;
using std::string;
using std::vector;
using std::endl;
using std::pair;
using std::make_pair;
using std::endl;

class DataCardCreatorHThTh_2016 {
	public:
  
  DataCardCreatorHThTh_2016(optutl::CommandLineParser parser) {
    channel_  = parser.stringValue("channel");
    
    string name_=channel_;
    filelabel_="";
    
    //define the histogram binning
    bins_ = parser.integerValue("bins");
    min_  = parser.doubleValue("min");
    max_  = parser.doubleValue("max");
    variable_      = parser.stringValue("variable");    
    
    preSelection_ = parser.stringValue("preselection");
    signalSelection_      = parser.stringValue("signalSelection");
    trigSelection_        = parser.stringValue("trigSelection");
    trigSelectionData_    = parser.stringValue("trigSelectionData");

    verbose_ = parser.integerValue("verbose");

    dir_ = parser.stringValue("dir");
    fout_ = new TFile(parser.stringValue("outputfile").c_str(),"RECREATE");

    luminosity_ = parser.doubleValue("luminosity");
  }
  
  
  void makeHiggsShape(string preselection, string categoryselection, string prefix) {
    
    pair<float,float> tmp;
    
    cout<<"Make Higgs Shape"<<endl;
    auto start = high_resolution_clock::now(); 

    std::string fullselection = "("+preselection+"&&"+signalSelection_+")";//*"+weight_;

    /*
    tmp= createHistogramAndShifts(dir_+"ggH120.root","ggH120",(fullselection),luminosity_,prefix);

    tmp= createHistogramAndShifts(dir_+"ggH125.root","ggH125",(fullselection),luminosity_,prefix);
    tmp= createHistogramAndShifts(dir_+"ggH130.root","ggH130",(fullselection),luminosity_,prefix);

    // if(doSys_>0)
    //      createShiftsTES("ggH125",dir_+"ggH125.root",categoryselection+"&&"+trigSelection_+"&&"+osSignalSelection_,weight_,luminosity_*legCorr,prefix,tmp);

    tmp= createHistogramAndShifts(dir_+"vbfH120.root","qqH120",(fullselection),luminosity_,prefix);
    tmp= createHistogramAndShifts(dir_+"vbfH125.root","qqH125",(fullselection),luminosity_,prefix);
    tmp= createHistogramAndShifts(dir_+"vbfH130.root","qqH130",(fullselection),luminosity_,prefix);
    
    //if(doSys_>0)
    //createShiftsTES("qqH125",dir_+"vbfH125.root",categoryselection+"&&"+trigSelection_+"&&"+osSignalSelection_,weight_,luminosity_*legCorr,prefix,tmp);
    
    tmp= createHistogramAndShifts(dir_+"ZH120.root","ZH120",(fullselection),luminosity_,prefix);
    tmp= createHistogramAndShifts(dir_+"ZH125.root","ZH125",(fullselection),luminosity_,prefix);
    tmp= createHistogramAndShifts(dir_+"ZH130.root","ZH130",(fullselection),luminosity_,prefix);
    
    */

    //if(doSys_>0)
    //createShiftsTES("ZH125",dir_+"ZH125.root",categoryselection+"&&"+trigSelection_+"&&"+osSignalSelection_,weight_,luminosity_,prefix,tmp);
    //    createShiftsTES("ZH125", dir_+"ZH125.root", "", weight_, luminosity_,prefix,tmp);
    /*
    tmp= createHistogramAndShifts(dir_+"WpH120.root","WpH120",(fullselection),luminosity_,prefix);
    tmp= createHistogramAndShifts(dir_+"WpH125.root","WpH125",(fullselection),luminosity_,prefix);
    tmp= createHistogramAndShifts(dir_+"WpH130.root","WpH130",(fullselection),luminosity_,prefix);
    
    tmp= createHistogramAndShifts(dir_+"WmH120.root","WmH120",(fullselection),luminosity_,prefix);
    tmp= createHistogramAndShifts(dir_+"WmH125.root","WmH125",(fullselection),luminosity_,prefix);
    tmp= createHistogramAndShifts(dir_+"WmH130.root","WmH130",(fullselection),luminosity_,prefix);
    */
    //createShiftsTES("WH125",dir_+"WH125.root",categoryselection+"&&"+trigSelection_+"&&"+osSignalSelection_,weight_,luminosity_,prefix,tmp);
    /*
    tmp= createHistogramAndShifts(dir_+"ttH120.root","ttH120",(fullselection),luminosity_,prefix);
    tmp= createHistogramAndShifts(dir_+"ttH125.root","ttH125",(fullselection),luminosity_,prefix);
    tmp= createHistogramAndShifts(dir_+"ttH130.root","ttH130",(fullselection),luminosity_,prefix);
    */
    // std::cout<<"creating met systematics Higgs"<<std::endl;
    // createMETSystematicsHiggs(fullselection, luminosity_, prefix);
    std::cout<<"creating jet systematics Higgs"<<std::endl;

    std::string trigselection = "("+trigSelection_+")";
    createJETSystematicsHiggs(fullselection, luminosity_, prefix);

    // Uncomment the following after making sure that the "signalSelection" cuts work
    //    createJETSystematicsHiggs(trigselection, luminosity_, prefix+"_dummy");
    //    std::string bothselection = "("+trigSelection_+"&&"+fullselection+")";
    //    createJETSystematicsHiggs(bothselection, luminosity_, prefix+"_dummy2");
    
    // Get ending timepoint 
    auto stop = high_resolution_clock::now(); 
  
    // Get duration. Substart timepoints to  
    // get duration. To cast it to proper unit 
    // use duration cast method 
    auto duration = duration_cast<seconds>(stop - start); 
  
    cout << "Time taken by function: "
         << duration.count() << " seconds" << endl; 
  }
  
  void close() {
    fout_->Close();
  }

  void createJETSystematicsHiggs(string inputSelections, float scale, string prefix){
    std::vector<std::string> jetSysVec = {"Closure", "AbsoluteFlavMap", "AbsoluteMPFBias", "AbsoluteScale", "AbsoluteStat", "FlavorQCD", "Fragmentation", "PileUpDataMC", "PileUpPtBB", "PileUpPtEC1", "PileUpPtEC2", "PileUpPtHF", "PileUpPtRef", "RelativeBal", "RelativeFSR", "RelativeJEREC1", "RelativeJEREC2", "RelativeJERHF", "RelativePtBB", "RelativePtEC1", "RelativePtEC2", "RelativePtHF", "RelativeStatEC", "RelativeStatFSR", "RelativeStatHF", "SinglePionECAL", "SinglePionHCAL", "TimePtEta", "Total"};
    for(auto jetSys : jetSysVec){

      std::string newSelectionUp=inputSelections;
      std::string newSelectionDown=inputSelections;
      ReplaceStringInPlace(newSelectionUp,   "njets", "njet_"   +jetSys+"Up");
      ReplaceStringInPlace(newSelectionUp,   "mjj"  , "vbfMass_"+jetSys+"Up");
      ReplaceStringInPlace(newSelectionDown, "njets", "njet_"   +jetSys+"Down");
      ReplaceStringInPlace(newSelectionDown, "mjj"  , "vbfMass_"+jetSys+"Down");
      
      pair<float,float> ggH125JetUp   = createHistogramAndShifts(   dir_+"ggH125.root" ,"ggH125_CMS_scale_j_"+jetSys+"_13TeVUp"  ,newSelectionUp,scale,prefix);
      pair<float,float> ggH125JetDown = createHistogramAndShifts(   dir_+"ggH125.root" ,"ggH125_CMS_scale_j_"+jetSys+"_13TeVDown",newSelectionDown,scale,prefix);

      
      pair<float,float> qqH125JetUp   = createHistogramAndShifts(  dir_+"vbfH125.root","qqH125_CMS_scale_j_" +jetSys+"_13TeVUp"  ,newSelectionUp,scale,prefix);
      pair<float,float> qqH125JetDown = createHistogramAndShifts(  dir_+"vbfH125.root","qqH125_CMS_scale_j_" +jetSys+"_13TeVDown",newSelectionDown,scale,prefix);
        
      pair<float,float> ZH125JetUp    = createHistogramAndShifts(  dir_+"ZH125.root"  ,"ZH125_CMS_scale_j_"  +jetSys+"_13TeVUp"  ,newSelectionUp,scale,prefix);
      pair<float,float> ZH125JetDown  = createHistogramAndShifts(  dir_+"ZH125.root"  ,"ZH125_CMS_scale_j_"  +jetSys+"_13TeVDown",newSelectionDown,scale,prefix);

      // File not found
      //      pair<float,float> WH125JetUp    = createHistogramAndShifts(  dir_+"WH125.root"  ,"WH125_CMS_scale_j_"  +jetSys+"_13TeVUp"  ,newSelectionUp,scale,prefix);
      //      pair<float,float> WH125JetDown  = createHistogramAndShifts(  dir_+"WH125.root"  ,"WH125_CMS_scale_j_"  +jetSys+"_13TeVDown",newSelectionDown,scale,prefix);
          
      pair<float,float> ttH125JetUp   = createHistogramAndShifts(  dir_+"ttH125.root" ,"ttH125_CMS_scale_j_" +jetSys+"_13TeVUp"  ,newSelectionUp,scale,prefix);
      pair<float,float> ttH125JetDown = createHistogramAndShifts(  dir_+"ttH125.root" ,"ttH125_CMS_scale_j_" +jetSys+"_13TeVDown",newSelectionDown,scale,prefix);


    }
  }

  /////////////////////////////////////////////////////

  void ReplaceStringInPlace(std::string& subject, const std::string& search,
			    const std::string& replace) {
    size_t pos = 0;
    while ((pos = subject.find(search, pos)) != std::string::npos) {
      subject.replace(pos, search.length(), replace);
      pos += replace.length();
    }
  }

  
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
  
  ///custom var
  //addme
  pair<float,float> createHistogramAndShiftsCustomVar(string variableIn, string file,string name, string cut,float scaleFactor = 1, string postfix = "",bool normUC  = true, bool keys=false,bool ShapeUncertainty=true) {
    TFile *f  = new TFile(file.c_str());
    if(f==0) printf("Not file Found\n");
    //get the nominal tree first
    TTree *t= (TTree*)f->Get((channel_+"EventTree/eventTree").c_str());
    if(t==0) printf("Not Tree Found in file %s\n",file.c_str());
    pair<float,float> yield;

    //    if(!keys)
      yield =makeHistogramCustomVar(variableIn, t,filelabel_+postfix,name,cut,scaleFactor);
      //    else
      //      yield =makeTHKeys(t,filelabel_+postfix,name,cut,scaleFactor);

    f->Close();
    return yield;
  }
  
  pair<float,float> makeHistogram(TTree* tree,string folder,TString name,string cut,float scaleFactor = 1.) {

    if(fout_->Get(folder.c_str())==0)
      fout_->mkdir(folder.c_str());

    TH1F *h=0;
    h= new TH1F(name,name,bins_,min_,max_);

    h->Sumw2();
    tree->Draw(variable_+">>"+name,cut.c_str());
    h->Scale(scaleFactor);

    fout_->cd(folder.c_str());

    Double_t error=0.0;
    
    float yield = h->IntegralAndError(1,h->GetNbinsX(),error,"");
    if(yield == 0){
      h->SetBinContent(1,0.00001);
    }
    
    h->Write(h->GetName(),TObject::kOverwrite);

    if(verbose_>0)
      cout << name<<": "<<h->Integral() << ", entries: " << h->GetEntries() <<endl;

    return make_pair(yield,error);
  }
  

  

  void createShiftsTES(string name, string inputFile, string inputSelections, string inputWeight, float scale, string prefix,  pair<float,float> nominalYield){

    std::string ptSelectionDM0Up_    = "(((pt_2*0.988)>40&&decayMode_2==0)||((pt_2*1.010)>40&&decayMode_2==1)||((pt_2*1.004)>40&&decayMode_2==10))&&(((pt_1*0.988)>50&&decayMode_1==0)||((pt_1*1.010)>50&&decayMode_1==1)||((pt_1*1.004)>50&&decayMode_1==10))";
    std::string ptSelectionDM0Down_  = "(((pt_2*0.976)>40&&decayMode_2==0)||((pt_2*1.010)>40&&decayMode_2==1)||((pt_2*1.004)>40&&decayMode_2==10))&&(((pt_1*0.976)>50&&decayMode_1==0)||((pt_1*1.010)>50&&decayMode_1==1)||((pt_1*1.004)>50&&decayMode_1==10))";
    std::string ptSelectionDM1Up_    = "(((pt_2*0.982)>40&&decayMode_2==0)||((pt_2*1.016)>40&&decayMode_2==1)||((pt_2*1.004)>40&&decayMode_2==10))&&(((pt_1*0.982)>50&&decayMode_1==0)||((pt_1*1.016)>50&&decayMode_1==1)||((pt_1*1.004)>50&&decayMode_1==10))";
    std::string ptSelectionDM1Down_  = "(((pt_2*0.982)>40&&decayMode_2==0)||((pt_2*1.004)>40&&decayMode_2==1)||((pt_2*1.004)>40&&decayMode_2==10))&&(((pt_1*0.982)>50&&decayMode_1==0)||((pt_1*1.004)>50&&decayMode_1==1)||((pt_1*1.004)>50&&decayMode_1==10))";
    std::string ptSelectionDM10Up_   = "(((pt_2*0.982)>40&&decayMode_2==0)||((pt_2*1.010)>40&&decayMode_2==1)||((pt_2*1.010)>40&&decayMode_2==10))&&(((pt_1*0.982)>50&&decayMode_1==0)||((pt_1*1.010)>50&&decayMode_1==1)||((pt_1*1.010)>50&&decayMode_1==10))";
    std::string ptSelectionDM10Down_ = "(((pt_2*0.982)>40&&decayMode_2==0)||((pt_2*1.010)>40&&decayMode_2==1)||((pt_2*0.998)>40&&decayMode_2==10))&&(((pt_1*0.982)>50&&decayMode_1==0)||((pt_1*1.010)>50&&decayMode_1==1)||((pt_1*0.998)>50&&decayMode_1==10))";

    pair<float,float> ZTT_DM0_UP  ;
    pair<float,float> ZTT_DM0_DOWN  ;
    pair<float,float> ZTT_DM1_UP    ;
    pair<float,float> ZTT_DM1_DOWN  ;
    pair<float,float> ZTT_DM10_UP   ;
    pair<float,float> ZTT_DM10_DOWN ;

    //std::cout<<"DM0Up Selection "<<"("+ptSelectionDM0Up_    +"&&"+nominalSelection_+"&&"+inputSelections+")*"+inputWeight<<std::endl;
    //    ZTT_DM0_UP    = createHistogramAndShiftsCustomVar(variable_+"_DM0_UP"   , inputFile,name+"_CMS_scale_t_1prong_13TeVUp",  ("("+ptSelectionDM0Up_    +"&&"+nominalSelection_+"&&"+inputSelections+")*"+inputWeight),scale,prefix);
    ZTT_DM0_UP    = createHistogramAndShiftsCustomVar((variable_+"_DM0_UP").Data()   , inputFile,name+"_CMS_scale_t_1prong_13TeVUp",  "("+ptSelectionDM0Up_    +"&&"+preSelection_+")", scale, prefix);
    ZTT_DM0_DOWN  = createHistogramAndShiftsCustomVar((variable_+"_DM0_DOWN").Data() , inputFile,name+"_CMS_scale_t_1prong_13TeVDown",  "("+ptSelectionDM0Down_ +"&&"+preSelection_+")", scale, prefix);
    ZTT_DM1_UP    = createHistogramAndShiftsCustomVar((variable_+"_DM1_UP").Data(), inputFile,name+"_CMS_scale_t_1prong1pizero_13TeVUp", "("+ptSelectionDM1Up_   + "&&"+preSelection_+")", scale, prefix);
    ZTT_DM1_DOWN  = createHistogramAndShiftsCustomVar((variable_+"_DM1_DOWN").Data(), inputFile,name+"_CMS_scale_t_1prong1pizero_13TeVDown", "("+ptSelectionDM1Down_  +"&&"+preSelection_+")", scale, prefix);
    ZTT_DM10_UP   = createHistogramAndShiftsCustomVar((variable_+"_DM10_UP").Data(), inputFile,name+"_CMS_scale_t_3prong_13TeVUp", "("+ptSelectionDM10Up_ + "&&"+preSelection_+")", scale, prefix);
    ZTT_DM10_DOWN = createHistogramAndShiftsCustomVar((variable_+"_DM10_DOWN").Data(), inputFile,name+"_CMS_scale_t_3prong_13TeVDown", "("+ptSelectionDM10Down_ + "&&"+preSelection_+")", scale, prefix);

    /*
    ZTT_DM0_DOWN  = createHistogramAndShiftsCustomVar(variable_+"_DM0_DOWN" , inputFile,name+"_CMS_scale_t_1prong_13TeVDown",       ("("+ptSelectionDM0Down_  +"&&"+nominalSelection_+"&&"+inputSelections+")*"+inputWeight),scale,prefix);
    ZTT_DM1_UP    = createHistogramAndShiftsCustomVar(variable_+"_DM1_UP"   , inputFile,name+"_CMS_scale_t_1prong1pizero_13TeVUp",  ("("+ptSelectionDM1Up_    +"&&"+nominalSelection_+"&&"+inputSelections+")*"+inputWeight),scale,prefix);
    ZTT_DM1_DOWN  = createHistogramAndShiftsCustomVar(variable_+"_DM1_DOWN" , inputFile,name+"_CMS_scale_t_1prong1pizero_13TeVDown",("("+ptSelectionDM1Down_  +"&&"+nominalSelection_+"&&"+inputSelections+")*"+inputWeight),scale,prefix);
    ZTT_DM10_UP   = createHistogramAndShiftsCustomVar(variable_+"_DM10_UP"  ,inputFile,name+"_CMS_scale_t_3prong_13TeVUp",         ("("+ptSelectionDM10Up_   +"&&"+nominalSelection_+"&&"+inputSelections+")*"+inputWeight),scale,prefix);
    ZTT_DM10_DOWN = createHistogramAndShiftsCustomVar(variable_+"_DM10_DOWN",inputFile,name+"_CMS_scale_t_3prong_13TeVDown",       ("("+ptSelectionDM10Down_ +"&&"+nominalSelection_+"&&"+inputSelections+")*"+inputWeight),scale,prefix);
						        

    //a quick sanity check... 
    std::cout<<"ZTT_DM0 nominal "<< nominalYield.first <<" Up "<<ZTT_DM0_UP.first<<" Down "<<ZTT_DM0_DOWN.first<<std::endl;
    if((abs(nominalYield.first-ZTT_DM0_UP.first)/nominalYield.first)>0.2)
      std::cout<<"Difference between Nominal and Up is greater than 20% this might indicate a problem"<<std::endl;
    if((abs(nominalYield.first-ZTT_DM0_DOWN.first)/nominalYield.first)>0.2)
      std::cout<<"Difference between Nominal and Down is greater than 20% this might indicate a problem"<<std::endl;
    */
  }

  pair<float,float> makeHistogramCustomVar(string variableIn, TTree* tree,string folder,string name,string cut, float scaleFactor = 1.) {

    if(fout_->Get(folder.c_str())==0)
      fout_->mkdir(folder.c_str());
    TH1F *h=0;


    //    if(binning_.size()==0)
    h= new TH1F(name.c_str(),name.c_str(),bins_,min_,max_);
    //else 
    //h = new TH1F(name.c_str(),name.c_str(),binning_.size()-1,&binning_[0]);
    h->Sumw2();

    tree->Draw((variableIn+">>"+name).c_str(),cut.c_str());

    h->Scale(scaleFactor);

    if(verbose_>0)
      cout<< " " <<name<<": "<<h->Integral()<<endl;

    fout_->cd(folder.c_str());

    Double_t error=0.0;
    //LD
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
  
  string trigSelection_;
  string trigSelectionData_;

  string signalSelection_;
  
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
