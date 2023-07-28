#include "functions.H"
/* Reading initial conditions */
/* Reading material properties */
Info<< "Reading solidProperties\n" << endl;
IOdictionary solidProperties
(
    IOobject
    (
        "solidProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    )
);


Info<< "Reading field T\n" << endl;
volScalarField T
(
    IOobject
    (
        "T",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
//    dimensionedScalar("T",dimensionSet(0,0,0,1,0,0,0),solidProperties.lookup("Tinit"))
);



Info<< "Reading initial liquid fraction alphaL \n" << endl;

volScalarField alphaL
(
    IOobject
    (
        "alphaL",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("alphaL",dimensionSet(0,0,0,0,0,0,0),solidProperties.lookup("alphaL"))
);

Info<< "Reading initial thermal diffusivity alphaT \n" << endl;

volScalarField alphaT
(
    IOobject
    (
        "alphaT",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("alphaT",dimensionSet(0,2,-1,0,0,0,0),solidProperties.lookup("alphaT"))
);


const scalar Tmelt
(
   readScalar(solidProperties.lookup("Tmelt"))
);

const scalar deltaT
(
   readScalar(solidProperties.lookup("deltaT"))
);
const scalar CpS
(
   readScalar(solidProperties.lookup("CpS"))
);
const scalar CpL
(
   readScalar(solidProperties.lookup("CpL"))
);
const scalar hf
(
   readScalar(solidProperties.lookup("hf"))
);
const scalar rhoS
(
   readScalar(solidProperties.lookup("rhoS"))
);
const scalar rhoL
(
   readScalar(solidProperties.lookup("rhoL"))
);

const scalar kS
(
   readScalar(solidProperties.lookup("kS"))
);
const scalar kL
(
   readScalar(solidProperties.lookup("kL"))
);

Info<< "Creating initial Enthalpy h \n" << endl;
volScalarField h

(
    IOobject
    (
        "h",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("h",dimensionSet(0,2,-2,0,0,0,0), 0.0)
);

Info<< "Updating initial Enthalpy h \n" << endl;
double T0 {};
double alphaLiq {};
double Cp {};
const cellList& cells_ = mesh.cells();
forAll(cells_, celli)
{
        T0 = T[celli];
	h[celli] = calcH( T0, Tmelt, deltaT, CpS, CpL, hf, alphaLiq, Cp); 
}
forAll(T.boundaryField(), patchi)
{
h.boundaryFieldRef().set(patchi, T.boundaryField()[patchi]);
    if (isA<fixedValueFvPatchScalarField>(T.boundaryField()[patchi]))
	    {
forAll(mesh.boundaryMesh()[patchi], facei)
{
    T0 = T.boundaryField()[patchi][facei];
    h.boundaryFieldRef()[patchi][facei] = calcH( T0, Tmelt, deltaT, CpS, CpL, hf, alphaLiq, Cp);
}
}
/*
else 
{
h.boundaryFieldRef().set(patchi, T.boundaryField()[patchi]);
}
*/
}
