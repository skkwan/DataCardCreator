
#include "DataCardCreator/StatTools/interface/DataCardCreatorHThTh_2016_RDF.h"
#include "PhysicsTools/FWLite/interface/CommandLineParser.h" 


int main (int argc, char* argv[]) 
{

	optutl::CommandLineParser parser ("Background subtraction ");

	//Input Files
	parser.addOption("channel",optutl::CommandLineParser::kString,"Channel  ","mutau");
	parser.addOption("preselection",     optutl::CommandLineParser::kString,"preselection",    "");
	parser.addOption("folder",optutl::CommandLineParser::kString,"folder","_inclusive");

	//Other Options
	parser.addOption("luminosity",optutl::CommandLineParser::kDouble,"Luminosity",10.);
	parser.addOption("variable",optutl::CommandLineParser::kString,"Shape variable ","mass");
	parser.addOption("weight",optutl::CommandLineParser::kString,"Weight for MC (Multiply Weight Factors here for efficiencies)","__WEIGHT__");
	parser.addOption("min",optutl::CommandLineParser::kDouble,"Minimum value",0.);
	parser.addOption("max",optutl::CommandLineParser::kDouble,"Maximum Value ",500.);
	parser.addOption("bins",optutl::CommandLineParser::kInteger,"Number of Bins",50);
	parser.addOption("verbose",optutl::CommandLineParser::kInteger,"verbose",0);
	parser.addOption("dir",optutl::CommandLineParser::kString,"dir","../inputs/mutau");

	//category options

	parser.parseArguments (argc, argv);

	DataCardCreatorHThTh_2016_RDF creator(parser);
	std::string inclSel = parser.stringValue("preselection"); 
	std::string foldername = parser.stringValue("folder");
	creator.makeHiggsShape(inclSel,inclSel,foldername);
	creator.close();
}
