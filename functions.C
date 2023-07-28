#include <cmath>
#include "functions.H"


double calcH(double T0, const double Tmelt, const double deltaT, const double CpS, const double CpL, const double hf, double& alphaLiq, double& Cp)
{
    const double Tsol = Tmelt - deltaT/2;
    const double Tliq = Tmelt + deltaT/2;
	
    alphaLiq = 0.5*tanh((T0 - Tmelt)/deltaT) + 0.5;
    double a = tanh((T0 - Tmelt)/deltaT); 	
    double dAlphaL = 0.5*(1 - a*a)/deltaT;
    double hsol = CpS*T0;
    double hliq = CpL*T0 + (CpS-CpL)*Tsol + hf ;
    //Cp = (1 - alphaLiq)*CpS + alphaLiq*CpL - dAlphaL*CpS*T0 + dAlphaL*(CpL*T0 + (CpS-CpL)*Tsol + hf);
    Cp = (1 - alphaLiq)*CpS + alphaLiq*CpL - dAlphaL*hsol + dAlphaL*hliq ;
    return  (alphaLiq)*hliq + (1 - (alphaLiq))*hsol;}

