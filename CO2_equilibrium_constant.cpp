#include <iostream>
#include <cmath>

typedef double Real;

double K3(const Real& p, const Real& newt, const Real omf[], Real& newk) {
    double y, d[13];

    d[0] = p / 1000000;  // Pressure (Pa)
    d[1] = newt;  // Temperature (C)
    d[2] = omf[1];  // Mole fraction of water in two-phase mixture
    d[3] = omf[2];  // Mole fraction of methane in two-phase mixture
    d[4] = omf[3];  // Mole fraction of carbon dioxide in two-phase mixture
    d[5] = omf[4];  // Mole fraction of nitrogen in two-phase mixture
    d[6] = omf[5];  // Mole fraction of hydrogen sulfide in two-phase mixture
    d[7] = omf[6];  // Mole fraction of NaCl in two-phase mixture
    d[8] = omf[7];  // Mole fraction of KCl in two-phase mixture
    d[9] = omf[8];  // Mole fraction of CaCl2 in two-phase mixture
    d[10] = omf[9];  // Mole fraction of MgCl2 in two-phase mixture

    const double G1C5 = -4.37535902279733;
    const double G1C9 = 7.77214880825221;
    const double G2C5 = -5.47728293710135;
    const double G3C1 = -1.0671278375179;
    const double G3C5 = 2.79641102328562;
    const double G4C2 = -0.785190099381902;
    const double G5C6 = -8.36603900265511;
    const double G5C7 = 0.728626256797193;
    const double G5C9 = -1.60924100466933;
    const double G6C8 = -1.47651261879235;
    const double G7C5 = -5.13412884914701;
    const double G7C3 = 4.56887319587694;
    const double G7C0 = -3.8067540910062;
    const double G8C8 = -0.418711059143904;
    const double G8C1 = 7.49136353832348;

    y = (((((d[1]+d[0])+(d[0]+G1C9))/2.0)-G1C5)+(((d[0]-d[2])-d[6])+sqrt(d[9])));
    y /= (d[7]+(((((G2C5+d[5])-(d[4]+d[4]))+sqrt(d[9]))+exp((d[4]+d[5])))/2.0));
    y /= (((1.0)+d[9])+G3C1);
    y /= ((pow(((d[9]+G4C2)/2.0),3)+(((((d[8]+d[7])/2.0)+(d[7]+d[8]))/2.0)+pow(d[7],2)))/2.0);
    y /= (((((exp((1.0))+(((0.0)+G5C7)/2.0))/2.0)*d[0])+(1.0/(G5C6)))/2.0);
    y /= ((((((G6C8+((d[5]+d[6])/2.0))/2.0)+(d[10]+d[10]))/2.0)+(d[9]+((d[4]+d[10])/2.0)))/2.0);
    y /= (((((exp(d[9])*G7C3)+exp(G7C0))+(G7C5+(d[8]+d[8])))/2.0)+d[10]);
    y /= (((((1.0/(((G8C1+(((1.0/(G8C8))+(0.0))/2.0))/2.0)))+G8C8)/2.0)+d[10])/2.0);

    newk = 0.0000295 * y + 10.24;

    return newk;
}

int main() {
    // Sample usage
    Real p = 20000000;
    Real newt = 50.0;
    Real omf[] = { 0.9, 0.05, 0.05, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
    Real newk;
    double result = K3(p, newt, omf, newk);
    std::cout << "Result: " << result << std::endl;

    return 0;
}
