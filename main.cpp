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
    // temperatures (in GeV)
    double const T_140MeV = 0.14, T_150MeV = 0.15;
    // kCuts
    int const kCut_Boltzmann = 1, kCut_Bessel = 30;

    // create hadrons
    // baryon octet
    Hadron const PROTON("proton", 0.9383, "fermion", 1, 1, 0, 2);
    Hadron const NEUTRON("neutron", 0.9396, "fermion", 1, 0, 0, 2);
    Hadron const SIGMA_MINUS("sigma-", 1.1975, "fermion", 1, -1, -1, 2);
    Hadron const SIGMA_ZERO("sigma0", 1.1926, "fermion", 1, 0, -1, 2);
    Hadron const SIGMA_PLUS("sigma+", 1.1894, "fermion", 1, 1, -1, 2);
    Hadron const LAMBDA("lambda", 1.11568, "fermion", 1, 0, -1, 2);
    Hadron const XI_MINUS("xi-", 1.32131, "fermion", 1, -1, -2, 2);
    Hadron const XI_ZERO("xi0", 1.31483, "fermion", 1, 0, -2, 2);
    // meson octet
    Hadron const KAON_ZERO("K0", 0.498, "boson", 0, 0, 1, 1);
    Hadron const KAON_PLUS("K+", 0.494, "boson", 0, 1, 1, 1);
    Hadron const PION_MINUS("pi-", 0.14, "boson", 0, -1, 0, 1);
    Hadron const PION_ZERO("pi0", 0.135, "boson", 0, 0, 0, 1);
    Hadron const PION_PLUS("pi+", 0.14, "boson", 0, 1, 0, 1);
    Hadron const ETA("eta", 0.547, "boson", 0, 0, 0, 1);
    Hadron const KAON_MINUS("K-", 0.494, "boson", 0, -1, -1, 1);
    Hadron const ANTI_KAON_ZERO("anti-K0", 0.498, "boson", 0, 0, -1, 1);

    // collect hadrons into a container
    std::vector<Hadron> hadronContainer{PROTON, NEUTRON, SIGMA_MINUS, SIGMA_ZERO, SIGMA_PLUS, LAMBDA, XI_MINUS, XI_ZERO,
                                        KAON_ZERO, KAON_PLUS, PION_MINUS, PION_ZERO, PION_PLUS, ETA, KAON_MINUS, ANTI_KAON_ZERO};

    // calculate pressures and energy densities
    std::vector<double> pressures_T140_Boltzmann{};
    std::vector<double> pressures_T140_Bessel{};
    std::vector<double> pressures_T150_Boltzmann{};
    std::vector<double> pressures_T150_Bessel{};
    std::vector<double> energyDensity_T140_Boltzmann{};
    std::vector<double> energyDensity_T140_Bessel{};
    std::vector<double> energyDensity_T150_Boltzmann{};
    std::vector<double> energyDensity_T150_Bessel{};
    for (int i = 0; i < static_cast<int>(hadronContainer.size()); i++)
    {
        pressures_T140_Boltzmann.push_back(iPartialPressure(T_140MeV, hadronContainer[i], kCut_Boltzmann));
        pressures_T140_Bessel.push_back(iPartialPressure(T_140MeV, hadronContainer[i], kCut_Bessel));
        pressures_T150_Boltzmann.push_back(iPartialPressure(T_150MeV, hadronContainer[i], kCut_Boltzmann));
        pressures_T150_Bessel.push_back(iPartialPressure(T_150MeV, hadronContainer[i], kCut_Bessel));
        energyDensity_T140_Boltzmann.push_back(iPartialEnergyDensity(T_140MeV, hadronContainer[i], kCut_Boltzmann));
        energyDensity_T140_Bessel.push_back(iPartialEnergyDensity(T_140MeV, hadronContainer[i], kCut_Bessel));
        energyDensity_T150_Boltzmann.push_back(iPartialEnergyDensity(T_150MeV, hadronContainer[i], kCut_Boltzmann));
        energyDensity_T150_Bessel.push_back(iPartialEnergyDensity(T_150MeV, hadronContainer[i], kCut_Bessel));
    }
    /*
    // write to screen
    std::cout << "NAME p_Boltzmann/p_Bessel(140MeV) p_Boltzmann/p_Bessel(150 MeV) eps_Boltzmann/eps_Bessel(140MeV) eps_Boltzmann/eps_Bessel(140MeV)" << std::endl;
    for (int i = 0; i < static_cast<int>(hadronContainer.size()); i++)
    {
        std::cout << hadronContainer[i].getName() << " " << std::fixed
                  << std::setprecision(6) << pressures_T140_Boltzmann[i] / pressures_T140_Bessel[i] << " "
                  << std::setprecision(6) << pressures_T150_Boltzmann[i] / pressures_T150_Bessel[i] << " "
                  << std::setprecision(6) << energyDensity_T140_Boltzmann[i] / energyDensity_T140_Bessel[i] << " "
                  << std::setprecision(6) << energyDensity_T150_Boltzmann[i] / energyDensity_T150_Bessel[i] << " "
                  << std::endl;
    }
    */
    // code testing
    // Ts
    std::vector<double> testTs{110., 120., 130., 140., 150., 160., 170., 180., 190., 200., 210., 220., 230., 240., 250.};
    // pressure test
    for (int i = 0; i < static_cast<int>(testTs.size()); i++)
    {
        std::cout << testTs[i] << " " << PressureHRG(hadronList, testTs[i] / 1000, 1000) / std::pow(testTs[i] / 1000, 4) << std::endl;
    }
}