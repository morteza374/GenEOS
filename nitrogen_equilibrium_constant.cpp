#include <iostream>
#include <cmath>

typedef double Real;

double K4(const Real& p, const Real& newt, const Real omf[], Real& newk) {
    double y, d[13];

    d[0] = p / 1000000;  // Pressure (Pa)
    d[1] = newt;  // Temperature (C)
    d[2] = omf[1];  // Mole fraction of water in two-phase mixture
    d[3] = omf[2];  // Mole fraction of methane in two-phase mixture
    d[4] = omf[3];  // Mole fraction of carbon dioxide in two-phase mixture
    d[5] = omf[4];  // Mole fraction of nitrogen in two-phase mixture
    d[6] = omf[5];  // Mole fraction of hydrogen sulphide in two-phase mixture
    d[7] = omf[6];  // Mole fraction of NaCl in two-phase mixture
    d[8] = omf[7];  // Mole fraction of KCl in two-phase mixture
    d[9] = omf[8];  // Mole fraction of CaCl2 in two-phase mixture
    d[10] = omf[9];  // Mole fraction of MgCl2 in two-phase mixture

    const double G2C5 = 1.67885372478408;
    const double G2C6 = -7.3702728075115;
    const double G3C1 = -1.07575217067748;
    const double G4C2 = -0.794497735526597;
    const double G4C4 = 9.19351823242659;
    const double G4C7 = 7.07973015589618;
    const double G5C6 = -5.60064394054994;
    const double G5C4 = -6.87862299487289;
    const double G6C8 = -0.145897135650236;
    const double G7C4 = -7.60565278701384;
    const double G7C7 = -9.86017840205084;
    const double G7C0 = -0.884851044673633;
    const double G8C8 = -0.427982260844974;
    const double G8C3 = 11.8674617343996;

    y = pow((pow(exp(exp(d[7])), 2) * (1.0 / ((((0.0) + (1.0)) / 2.0)))), 3);
    y /= ((pow(((((d[5] + sqrt(d[8])) / 2.0) + ((G2C6 / G2C5) / G2C5)) / 2.0), 3) + sqrt(d[9])) / 2.0);
    y /= ((G3C1 + exp(d[8])) / 2.0);
    y /= ((d[9] * sqrt((((G4C4 * G4C7) / d[2]) + pow(d[6], 2)))) - pow(G4C2, 2));
    y /= ((d[0] * G5C6) + ((-((d[5] + d[5]))) + (1.0 / (((G5C4 + d[1]) / 2.0)))));
    y /= pow(((d[10] + G6C8) / 2.0), 3);
    y /= (((sqrt(d[2]) + ((G7C4 + G7C7) * (d[9] + d[10]))) + exp(((G7C0 - d[3]) - d[8])))/2.0);
    y /= (((G8C8 + ((((G8C8 + d[9]) + (1.0)) / 2.0) - (1.0 / (((G8C3 + d[1]) / 2.0)))))/2.0) + d[7]);

    newk = 0.0000105 * y + 82.87;

    return newk;
}

int main() {
    // Sample usage
    Real p = 20000000;
    Real newt = 70.0;
    Real omf[] = { 0.9, 0.05, 0.0, 0.05, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
    Real newk;
    double result = K4(p, newt, omf, newk);
    std::cout << "Result: " << result << std::endl;

    return 0;
}
