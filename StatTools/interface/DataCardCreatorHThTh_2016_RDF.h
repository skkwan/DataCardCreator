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
    ROOT::EnableImplicitMT();

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
    /*
    tmp= createHistogramAndShifts(dir_+"ZH120.root","ZH120",(fullselection),luminosity_,prefix);
    tmp= createHistogramAndShifts(dir_+"ZH125.root","ZH125",(fullselection),luminosity_,prefix);
    tmp= createHistogramAndShifts(dir_+"ZH130.root","ZH130",(fullselection),luminosity_,prefix);
    */
    //if(doSys_>0)
    //createShiftsTES("ZH125",dir_+"ZH125.root",categoryselection+"&&"+trigSelection_+"&&"+osSignalSelection_,weight_,luminosity_,prefix,tmp);
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
    //std::cout<<"creating met systematics Higgs"<<std::endl;
    //createMETSystematicsHiggs(fullselection, luminosity_, prefix);
    //std::cout<<"creating jet systematics Higgs"<<std::endl;
    //createJETSystematicsHiggs(fullselection, luminosity_, prefix);

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


  pair<float,float> createHistogramAndShifts(string file,string name, string cut, float scaleFactor = 1, string postfix = "",bool normUC  = true, bool keys=false) {
    string folder = filelabel_+postfix;

    ROOT::EnableImplicitMT();

    auto startCreate = high_resolution_clock::now();

    auto start = high_resolution_clock::now();

    ROOT::RDataFrame d((channel_+"EventTree/eventTree").c_str(), file);


    
    auto end = high_resolution_clock::now();
    auto duration = duration_cast<microseconds>(end - start);
    cout << "Time taken to read dataframe: "
         << duration.count() << " microseconds" << endl;

    auto cbHisto = high_resolution_clock::now();
    auto h = d.Filter(cut).Histo1D({name.data(), name.data(), bins_, min_, max_}, "m_sv");

    auto startSumW2 = high_resolution_clock::now();
    h->Sumw2();
    auto endSumW2 = high_resolution_clock::now();
    duration = duration_cast<microseconds>(endSumW2 - startSumW2);
    cout << "Time for sumw2: " << duration.count() << "microseconds" << endl;

    // Alternatively:
    // auto dCut = d.Filter(cut);
    // auto h = dCut.Histo1D({name.data(), name.data(), bins_, min_, max_}, "m_sv");

    auto caHisto = high_resolution_clock::now();
    // END

    duration = duration_cast<microseconds>(caHisto - cbHisto);
    cout << "Time taken to Filter and call Histo1D: "
         << duration.count() << "  microseconds" << endl;
    
    if(fout_->Get(folder.c_str())==0)
      fout_->mkdir(folder.c_str());
    fout_->cd(folder.c_str());


    /*
    
    auto startSumW2 = high_resolution_clock::now();
    h->Sumw2();
    auto endSumW2 = high_resolution_clock::now();

    duration = duration_cast<microseconds>(endSumW2 - startSumW2);
    cout << "Time taken to SumW2: " 
	 << duration.count() << "  microseconds" << endl;
    */
    // START
    auto cbWrite = high_resolution_clock::now();
    h->Write(name.data(), TObject::kOverwrite);
    
    auto caWrite = high_resolution_clock::now();
    duration = duration_cast<microseconds>(caWrite - cbWrite);
    // END
    cout << "Time taken to write histo: "
         << duration.count() << "  microseconds" << endl;

    Double_t error=0.0;
    
    float yield = h->IntegralAndError(1,h->GetNbinsX(),error,"");

    if(yield == 0){
      h->SetBinContent(1,0.00001);
    }
			  
    auto endCreate = high_resolution_clock::now();
    duration = duration_cast<microseconds>(endCreate - startCreate);
    cout << "Time inside create: "
	 << duration.count() <<"  microseconds" << endl;

    return make_pair(yield,error);
  }

  /**********************************************************************/

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
    ZTT_DM0_UP    = createHistogramAndShiftsCustomVar(variable_+"_DM0_UP"   , inputFile,name+"_CMS_scale_t_1prong_13TeVUp",         ("("+ptSelectionDM0Up_    +"&&"+nominalSelection_+"&&"+inputSelections+")*"+inputWeight),scale,prefix);
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
