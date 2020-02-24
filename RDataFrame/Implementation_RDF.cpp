/******************************************************/
/* Implementation_RDF.cpp                             */
/* Usage: root -l Implementation_RDF.cpp              */
/******************************************************/

/******************************************************/

#include <ROOT/TDataFrame.hxx>
using namespace std;
using namespace ROOT::Experimental;

/******************************************************/

int Implementation_RDF()
{
  // Directory
  auto dir_ =  "/hdfs/store/user/ojalvo/ztt_weighted_Apr7_svFit/ggH125.root";
  auto treeName = "diTauEventTree/eventTree";

  // Create a RDataFrame, a class that allows us to interact with the data
  // contained in the tree.
  EnableImplicitMT();
  TDataFrame d(treeName, fileName);
  
  // Example of a simple cut
  auto preselection = "isoTight_1>0.5&&isoTight_2>0.5&& abs(eta_1)<2.1&&abs(eta_2)<2.1&&  \
                       againstMuonTight3_1>0&&againstMuonTight3_2>0&& \
                       againstElectronVLooseMVA6_1>0&&againstElectronVLooseMVA6_2>0&& \
                       npv>0&&diLeptons==0&&extraelec_veto==0&&extramuon_veto==0";
                         
  auto cut = d.Filter(preselection);
  auto nEntries = cut.Count();
  std::cout << *nEntries << " passed all filters" << std::endl;

  // Draw the results
  auto hCut = cut.Histo1D("m_sv");
  auto c = new TCanvas();
  hCut->DrawCopy();

  return 0;
}






/******************************************************/
