#include <iostream>
#include <cmath>

typedef double Real;

double K2(const Real &p, const Real &newt, const Real omf[], Real &newk) {
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

    const double G2C5 = 2.22339771721549;
    const double G2C6 = -2.14637233854183;
    const double G3C1 = -1.07790690658374;
    const double G4C2 = -0.794709937672379;
    const double G4C4 = 9.25857853236488;
    const double G4C9 = 3.02574550033259;
    const double G5C8 = -7.32716149082919;
    const double G5C7 = -4.92416150395215;
    const double G6C2 = -1.10201172067956e-03;
    const double G6C8 = -0.146189368197262;
    const double G7C8 = -9.06125064851833;
    const double G7C0 = -0.732014600985748;
    const double G8C8 = -0.419130189333237;
    const double G8C6 = -57.7562791833247;
    const double G8C3 = 11.7499738632408;

    y = (-(pow(exp((d[8] * exp(pow((d[2] + (-(d[3]))), 3)))), 2)));
    y /= (pow(((((((d[7] + G2C6) / 2.0) + d[5]) / 2.0) + ((d[3] + G2C6) - G2C5)) / 2.0), 3) + sqrt(d[0]));
    y /= ((G3C1 + exp(d[8])) / 2.0);
    y /= ((d[9] * sqrt((((G4C9 - d[2]) * d[0]) + pow(G4C4, 2)))) - pow(G4C2, 2));
    y /= ((d[0] * G5C8) + exp(((1.0 / (((d[4] + d[1]) / 2.0))) + ((((G5C7 + d[2]) / 2.0) + d[6]) / 2.0))));
    y /= pow((-(((G6C2 + ((G6C8 + d[10]) / 2.0)) / 2.0))), 3);
    y /= (((sqrt(d[2]) + ((G7C8 + G7C8) * (d[9] + d[10]))) + exp((pow(G7C0, 3) - d[5]))) / 2.0);
    y /= (((((((d[9] + G8C8) + (1.0)) / 2.0) - (1.0 / (((d[1] + G8C3) / 2.0)))) + G8C8) / 2.0) + d[7]);

    newk = 0.040489 * y + 98.83;
    
    return newk;
}

int main() {
    // Sample usage
    Real p = 20000000;
    Real newt = 70.0;
    Real omf[] = { 0.9, 0.05, 0.05, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
    Real newk;
    double result = K2(p, newt, omf, newk);
    std::cout << "Result: " << result << std::endl;

    return 0;
}
