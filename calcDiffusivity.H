#include<cmath>
#include "functions.H"
//Apply temperature Correction
const scalar tol = 1e-12;
scalar error = tol + 1;
double T0 = Tmelt;
const int max_it = 100;
int it = 0;
scalar k = kS;
scalar rho = rhoS;

const cellList& cells_ = mesh.cells();

forAll(cells_, celli)
{
    const scalar& hc = h[celli];
    double hNew {};
    double TNew {} ;
    // Calculate temperature corresponding to enthalpy
    while (error > tol && it < max_it)
    { 
	hNew = calcH( T0, Tmelt, deltaT, CpS, CpL, hf, alphaLiq, Cp); 
        TNew = T0 - (hNew - hc) / Cp;
        error = fabs(TNew - T0);
        T0 = TNew;
        it++;
    }
    // update corrected temperature
    T[celli] = TNew;
    // update liquid fraction
    alphaL[celli] = alphaLiq; 
    // update thermal diffusivity
    k = (1 - alphaLiq)*kS + alphaLiq*kL;
    rho = (1 - alphaLiq)*rhoS + alphaLiq*rhoL;
    alphaT[celli] = k/(rho*Cp); 
    // reset error and iterator
    it = 0;
    error = tol + 1;
    T0 = Tmelt;
}

