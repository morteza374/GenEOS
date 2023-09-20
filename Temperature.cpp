#include <iostream>

typedef double Real;

double Temperature(const Real& p, const Real& h, const Real omf[], Real& ct) {
    double y, d[13];

    d[0] = p;  // Pressure (Pa)
    d[1] = h / 1000000.0;  // Enthalpy (J/Kg/K)
    d[2] = omf[1];  // Mole fraction of water in two-phase mixture
    d[3] = omf[2];  // Mole fraction of methane in two-phase mixture
    d[4] = omf[3];  // Mole fraction of carbon dioxide in two-phase mixture
    d[5] = omf[4];  // Mole fraction of nitrogen in two-phase mixture
    d[6] = omf[5];  // Mole fraction of hydrogen sulphide in two-phase mixture
    d[7] = omf[6];  // Mole fraction of NaCl in two-phase mixture
    d[8] = omf[7];  // Mole fraction of KCl in two-phase mixture
    d[9] = omf[8];  // Mole fraction of CaCl2 in two-phase mixture
    d[10] = omf[9];  // Mole fraction of MgCl2 in two-phase mixture

    if (p >= 10000000) {
        const double G1C1 = 20468.3987895087;
        const double G1C8 = -0.632584680896207;
        const double G1C4 = -0.126694048571721;
        const double G1C0 = -2.34224648339043;
        const double G2C6 = -4.6741958733112;
        const double G2C5 = 26.4014636827949;
        const double G2C1 = 21.0155662431996;
        const double G3C6 = 1.38009587889663;
        const double G3C7 = 11.0276506089214;
        const double G4C0 = 6.40107789479513;
        const double G4C7 = -15.3514731948368;

        y = ((d[10] + ((d[2] * G1C8) + (d[7] * G1C4))) * (((d[1] + G1C0) + d[3]) * (G1C1 + d[0])));
        y /= ((((d[1] + d[7]) + (d[7] + d[7])) * ((G2C5 + G2C1) * d[8])) + ((d[8] + G2C6) + d[6]));
        y /= (d[0] + ((((d[9] * G3C7) * d[0]) * (d[2] * G3C6)) * ((d[9] + d[1]) + (d[3] + d[6]))));
        y /= (d[1] * ((d[1] + (((d[10] * G4C7) * d[10]) * (G4C0 + d[1]))) + d[2]));

        ct = -0.00080757 * y + 0.00057108;
    }
    if (p < 10000000) {
        const double G1C1 = -32318.1483269907;
        const double G1C8 = -13.8547740182964;
        const double G1C7 = -3.88915202762322;
        const double G2C3 = 6.31536021447302;
        const double G2C5 = 4.0254239988184;
        const double G2C2 = -1.73940616972321;
        const double G3C4 = 7.38868080049002;
        const double G3C2 = 6.19566639338576;
        const double G3C8 = 5.98974791711173;
        const double G3C0 = -5.72206137923993;
        const double G4C2 = 2.98313090164621;

        y = (((d[4] + (G1C8 * d[8])) + ((G1C7 + d[3]) + d[3])) * ((G1C1 + d[2]) + (d[2] * d[0])));
        y /= ((((d[6] + d[6]) * (d[3] + d[1])) * ((d[2] * d[2]) * G2C3)) + ((d[4] + G2C5) + (d[7] * G2C2)));
        y /= (d[0] + ((((G3C4 + G3C2) * d[0]) * ((G3C8 + G3C0) + d[3])) * d[9]));
        y /= (((((d[2] * d[5]) + d[8]) * (d[8] * G4C2)) + d[1]) * (d[1] + d[2]));

        ct = -0.0002703 * y + 0.00050082;
    }

    return ct*h;
}

int main() {
    // Sample usage
    Real p = 10000000;
    Real h = 600000;
    Real omf[] = { 0.9, 0.05, 0.05, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
    Real ct;
    double result = Temperature(p, h, omf, ct);
    std::cout << "Result: " << result << std::endl;

    return 0;
}
