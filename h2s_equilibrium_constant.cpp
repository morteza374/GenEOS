#include <iostream>
#include <cmath>

typedef double Real;

double K5(const Real& p, const Real& newt, const Real omf[], Real& newk) {
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

    const double G1C7 = 3.33536790063173;
    const double G1C6 = 5.94302407910398;
    const double G2C0 = -7.01362340355867;
    const double G3C0 = 4.18216557671013;
    const double G3C2 = -7.15920981856555;
    const double G3C8 = 4.65926084170049;
    const double G4C1 = -0.759923389690848;
    const double G4C7 = 7.49086935868265;
    const double G4C3 = 8.43005236884671;
    const double G4C5 = -9.50941801202429;
    const double G5C4 = 0.936272689231344;
    const double G6C6 = 8.64558854945524;
    const double G6C8 = -0.144408844877855;
    const double G7C7 = 4.36078981902524;
    const double G7C8 = -10.6495893859506;
    const double G8C1 = 6.26367671004059;

    y = exp((((G1C7 + sqrt(d[8])) / 2.0) + exp((((G1C6 + d[7]) / 2.0) + d[8]))));
    y /= (exp(((d[6] + (d[2] - (((d[0] * d[6]) + G2C0) / 2.0))) / 2.0)) + sqrt(d[0]));
    y /= (((((G3C0 - G3C2) * (d[9] + d[7])) + (((G3C8 + d[6]) + (d[7] + G3C2)) / 2.0)) / 2.0) + d[10]);
    y /= ((d[9] * (((d[6] + G4C7) + ((G4C3 + G4C5) / 2.0)) + (d[0] * d[6]))) - pow(G4C1, 2));
    y /= (d[9] - (d[0] / sqrt((pow(((d[5] / d[0]) + G5C4), 3) - d[5]))));
    y /= ((1.0 / G6C6) * pow((((d[10] * exp((d[4] + d[6]))) + G6C8) / 2.0), 3));
    y /= ((((1.0) + ((d[2] + G7C8) * (d[9] + d[10]))) / 2.0) / exp(((d[9] * d[4]) * d[1])));
    y /= (1.0 / (((G8C1 + ((((d[1] + d[0]) / 2.0) + (1.0 / (d[1]))) / 2.0)) / 2.0) - ((-(d[0])) + d[6])));

    newk = 3.07E-13 * y + 2.98;

    return newk;
}

int main() {
    // Sample usage
    Real p = 20000000;
    Real newt = 30.0;
    Real omf[] = { 0.9, 0.05, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.05 };
    Real newk;
    double result = K5(p, newt, omf, newk);
    std::cout << "Result: " << result << std::endl;

    return 0;
}
