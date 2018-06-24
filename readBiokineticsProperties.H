Info<< "Reading biokineticsProperties" << endl;

IOdictionary biokineticsProperties
(
    IOobject
    (
        "biokineticsProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    )
);

Info << "  Reading stoichiometric parameters" << endl;

Info << "    Reading autotrophic yield, YA" << endl;

dimensionedScalar YA
(   
    biokineticsProperties.lookup("YA")
);

Info << "    Reading heterotrophic yield, YH" << endl;

dimensionedScalar YH
(
    biokineticsProperties.lookup("YH")
);

Info << "    Reading fraction of biomass leading to particulate products, fP" << endl;

dimensionedScalar fP
(
    biokineticsProperties.lookup("fP")
);

Info << "    Reading mass of nitrogen per mass of COD in biomass, iXB" << endl;

dimensionedScalar iXB
(
    biokineticsProperties.lookup("iXB")
);

Info << "    Reading mass of nitrogen per mass of COD in products from biomass, iXP" << endl;

dimensionedScalar iXP
(
    biokineticsProperties.lookup("iXP")
);

Info << "  Reading kinetic parameters" << endl;

Info << "    Reading maximum specific growth rate for heterotrophic biomass, muH" << endl;

dimensionedScalar muH
(
    biokineticsProperties.lookup("muH")
);

Info << "    Reading half-saturation coefficient for heterotrophic biomass, KS" << endl;

dimensionedScalar KS
(
    biokineticsProperties.lookup("KS")
);

Info << "    Reading oxygen half-saturation coefficient for heterotrophic biomass, KOH" << endl;

dimensionedScalar KOH
(
    biokineticsProperties.lookup("KOH")
);

Info << "    Reading nitrate half-saturation coefficient for denitrifying heterotrophic biomass, KNO" << endl;

dimensionedScalar KNO
(
    biokineticsProperties.lookup("KNO")
);

Info << "    Reading decay coefficient for heterotrophic biomass, bH" << endl;

dimensionedScalar bH
(
    biokineticsProperties.lookup("bH")
);

Info << "    Reading maximum specific growth rate for autotrophic biomass, muA" << endl;

dimensionedScalar muA
(
    biokineticsProperties.lookup("muA")
);

Info << "    Reading ammonia half-saturation coefficient for autotrophic biomass, KNH" << endl;

dimensionedScalar KNH
(
    biokineticsProperties.lookup("KNH")
);

Info << "    Reading oxygen half-saturation coefficient for autotrophic biomass, KOA" << endl;

dimensionedScalar KOA
(
    biokineticsProperties.lookup("KOA")
);

Info << "    Reading decay coefficient for autotrophic biomass, bA" << endl;

dimensionedScalar bA
(
    biokineticsProperties.lookup("bA")
);

Info << "    Reading correction factor for muH under anoxic conditions, etag" << endl;

dimensionedScalar etag
(
    biokineticsProperties.lookup("etag")
);

Info << "    Reading ammonification rate, ka" << endl;

dimensionedScalar ka
(
    biokineticsProperties.lookup("ka")
);

Info << "    Reading maximum specific hydrolysis rate, kh" << endl;

dimensionedScalar kh
(
    biokineticsProperties.lookup("kh")
);

Info << "    Reading half-saturation coefficient for hydrolysis of slowly biodegradable substrate, KX" << endl;

dimensionedScalar KX
(
    biokineticsProperties.lookup("KX")
);

Info << "    Reading correction factor for hydrolysis under anoxic conditions, etah" << endl;

dimensionedScalar etah
(
    biokineticsProperties.lookup("etah")
);

Info << "    Reading oxygen mass transfer coefficient, kLa" << endl;

dimensionedScalar kLa
(
    biokineticsProperties.lookup("kLa")
);

Info << "    Reading oxygen saturation concentration, SOsat" << endl;

dimensionedScalar SOsat
(
    biokineticsProperties.lookup("SOsat")
);


Info << "    Reading turbulent Schmidt number, ScT" << endl;
scalar ScT
(
    readScalar(biokineticsProperties.lookup("ScT"))
);

Info << "    Reading liquid phase name, liquidName" << endl;
word liquidName
(
    biokineticsProperties.lookup("liquidName")
);

Info << endl;
