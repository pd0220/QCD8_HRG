// including used header
#include "AnalysisTools.hh"

// ------------------------------------------------------------------------------------------------------------

// PDF file name (hadron list)
std::string const PDG = "../PDG.txt";

// ------------------------------------------------------------------------------------------------------------

// main function
int main(int, char **)
{
    // hadrons from PDG
    std::vector<Hadron> hadronList = HadronList(PDG);
    // temperature (GeV)
    double T = .15;
    // kCut
    int kCut = 2;
    // values at mu = 0
    std::cout << "pressure: " << PressureHRG(hadronList, T, kCut) << std::endl;
    std::cout << "energy density: " << EnergyDensityHRG(hadronList, T, kCut) << std::endl;
}