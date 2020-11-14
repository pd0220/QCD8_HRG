// including used header(s)
#include "AnalysisTools.hh"

// ------------------------------------------------------------------------------------------------------------

// PDG file name (hadron list)
std::string const PDG = "../PDG.txt";

// ------------------------------------------------------------------------------------------------------------

// main function
int main(int, char **)
{
    // hadrons from PDG
    std::vector<Hadron> hadronList = HadronList(PDG);

    // where to cut the summation in the partial pressure
    int kCut = 1000;

    // temperature values (MeV)
    std::vector<double> Ts{110.0, 115.0, 120.0, 125.0, 130.0, 135.0, 140.0, 145.0, 150.0, 155.0, 160.0, 165.0, 170.0, 175.0, 180.0, 185.0,
                           190.0, 195.0, 200.0, 205.0, 210.0, 215.0, 220.0, 225.0, 230.0, 235.0, 240.0, 245.0, 250.0, 260.0, 270.0, 280.0,
                           290.0, 300.0, 310.0, 320.0, 330.0, 340.0, 350.0, 360.0, 370.0, 380.0, 390.0, 400.0, 410.0, 420.0, 430.0, 440.0,
                           445.0, 450.0, 455.0, 460.0, 465.0, 470.0, 475.0, 480.0, 485.0, 490.0, 495.0, 500.0, 505.0, 510.0};

    // vector container for different {B, S} sector contributions at different temperatures
    std::vector<Eigen::MatrixXd> SectorMatrixTotal{Ts.size(), Eigen::MatrixXd::Zero(4, 5)};
    // vector container for baryon number sector contributions (B = 1, 2)
    std::vector<std::vector<double>> BaryonVectorTotal{Ts.size(), std::vector<double>{0., 0.}};
    // hadron matrix for {B, S} sector contributions at fixed temperature
    Eigen::MatrixXd SectorMatrixTemperature;
    // loop for temperatures
    for (int TIndex = 0; TIndex < static_cast<int>(Ts.size()); TIndex++)
    {
        // temperature [MeV --> GeV]
        double T = Ts[TIndex] / 1000.;
        // set matrix elements to zero
        SectorMatrixTemperature = Eigen::MatrixXd::Zero(4, 5);
        // HRG partial pressure rewritten to acces kth contribution
        // partial pressure at mu = 0
        // loop for hadrons
        for (int hadronIndex = 0; hadronIndex < static_cast<int>(hadronList.size()); hadronIndex++)
        {
            // hadron
            Hadron hadron = hadronList[hadronIndex];
            // baryon number and strangeness
            int baryionNumber = hadron.getB();
            int strangeness = -hadron.getS();

            // particle data: mass, eta (boson / fermion), spin degeneracy
            double mass = hadron.getMass();
            int eta = EtaDetermination(hadron);
            int spinDeg = hadron.getSpinDegeneracy();

            // prefactor of Macdonald function
            double prefactor = spinDeg * sq(T * mass / M_PI) / 2;

            // loop for k index in partial pressure (baryion number only)
            for (int k = 1; k < kCut; k++)
            {
                // determine what sector to update
                int sectorB = k * baryionNumber;
                // only consider B = 1, 2
                if (sectorB > 2 || sectorB < 1)
                    break;

                // argument of Macdonald function
                double argument = k * mass / T;
                // calculate kth contribution and update matrix element
                BaryonVectorTotal[TIndex][sectorB - 1] += prefactor * std::pow(-eta, k + 1) / sq(k) * gsl_sf_bessel_Kn(2, argument);
            }

            // loop for k index in partial pressure (baryon strangeness together)
            for (int k = 1; k < kCut; k++)
            {
                // determine what sector to update
                int sectorB = k * baryionNumber;
                int sectorS = k * strangeness;
                // only consider 0 < B * k < 4 and -2 < S * k < 4 now
                if (sectorB > 3 || sectorB < 0 || sectorS > 3 || sectorS < -1)
                    break;
                // matrix cannot have negative index --> last index will belong to S = -1 sector
                if (sectorS == -1)
                    sectorS = 4;

                // argument of Macdonald function
                double argument = k * mass / T;
                // prefactor of Macdonald function
                double prefactor = spinDeg * sq(T * mass / M_PI) / 2;
                // calculate kth contribution and update matrix element
                SectorMatrixTemperature(sectorB, sectorS) += prefactor * std::pow(-eta, k + 1) / sq(k) * gsl_sf_bessel_Kn(2, argument);
            }
        }

        // add matrix to matrix container
        SectorMatrixTotal[TIndex] = SectorMatrixTemperature;
    }

    // write matrix to file
    // structure
    //
    // T, P_0-1, P_00, P_01, P_02, P_03, P_10, P_11, P_12, P_13, P_20, P_21, P_22, P_23, P_30, P_31, P_32, P_33
    /*
    std::ofstream file;
    file.open("HRGSectors.txt", std::ofstream::app);
    for (int iData = 0; iData < static_cast<int>(SectorMatrixTotal.size()); iData++)
    {
        file << Ts[iData] << " "
             << SectorMatrixTotal[iData](0, 4) << " ";
        for (int row = 0; row <= 3; row++)
        {
            for (int col = 0; col <= 3; col++)
            {
                file << SectorMatrixTotal[iData](row, col) << " ";
            }
        }
        file << std::endl;
    }
    file.close();
*/
    // write baryon vector to file
    // structure
    //
    // T, P_1, P_2
    std::ofstream fileBaryon;
    fileBaryon.open("BaryonSectors.txt", std::ofstream::app);
    for (int iData = 0; iData < static_cast<int>(BaryonVectorTotal.size()); iData++)
    {
        fileBaryon << Ts[iData] << " " << BaryonVectorTotal[iData][0] << " " << BaryonVectorTotal[iData][1] << std::endl;
    }
    fileBaryon.close();
}

/*
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
