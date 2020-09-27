// implementing hadron class to simplify HRG quantity calculation

// used headers/libraries
#include <iostream>
#include <string>

// hadron class
class Hadron
{
    // variables
    // name
    std::string name;
    // mass (in GeV)
    double mass;
    // type (boson/fermion)
    std::string type;
    // baryon number
    int B;
    // electric charge
    int Q;
    // strangeness
    int S;
    // spin degeneracy
    int g;

    //constructors
public:
    // default ~ setting everything to zero (physically bullshit but something)
    Hadron() : name{"None"}, mass{0.}, type{"boson"}, B{0}, Q{0}, S{0}, g{0} {}
    // parameterized default
    Hadron(std::string name, double mass, std::string type, int B, int Q, int S, int g) : name{name}, mass{mass}, type{type}, B{B}, Q{Q}, S{S}, g{g} {}
    // copy
    Hadron(Hadron const &) = default;
    // copy assigment
    Hadron &operator=(Hadron const &) = default;

    // inner functions
    // get name of hadron
    std::string getName() const
    {
        return (*this).name;
    }
    // get mass of hadron
    double getMass() const
    {
        return (*this).mass;
    }
    // get type of hadron
    std::string getType() const
    {
        return (*this).type;
    }
    // get baryon number of hadron
    int getB() const
    {
        return (*this).B;
    }
    // get electric charge of hadron
    int getQ() const
    {
        return (*this).Q;
    }
    // get strangeness of hadron
    int getS() const
    {
        return (*this).S;
    }
    // get spin degeneracy of hadron
    int getSpinDegeneracy() const
    {
        return (*this).g;
    }
    // write hadron data to screen
    void hadronData() const
    {
        std::cout << (*this).getName() << std::endl;
        std::cout << "mass: " << (*this).getMass() << " GeV" << std::endl;
        std::cout << "type: " << (*this).getType() << std::endl;
        std::cout << "baryon number: " << (*this).getB() << std::endl;
        std::cout << "electric charge: " << (*this).getQ() << std::endl;
        std::cout << "strangeness: " << (*this).getS() << std::endl;
        std::cout << "spin degeneracy: " << (*this).getSpinDegeneracy() << std::endl;
    }
};
