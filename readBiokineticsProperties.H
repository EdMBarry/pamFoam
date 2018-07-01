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

Info << "    Reading photoheterotrophic yield, YPBPH" << endl;

scalar YPBPH
(
    readScalar(biokineticsProperties.subDict("stoichiometry").lookup("YPBPH"))
);

Info << "    Reading mass fraction of soluble organic substrate SS in XS" << endl;
scalar fSSXS
(
    readScalar(biokineticsProperties.subDict("stoichiometry").lookup("fSSXS"))
);

Info << "    Reading mass fraction of acetate substrate SAC in XS" << endl;
scalar fSACXS
(
    readScalar(biokineticsProperties.subDict("stoichiometry").lookup("fSACXS"))
);

Info << "    Reading mass fraction of inorganic nitrogen SIN in XS" << endl;
scalar fSINXS
(
    readScalar(biokineticsProperties.subDict("stoichiometry").lookup("fSINXS"))
);

Info << "    Reading mass fraction of inorganic phosphorus SIP in XS" << endl;
scalar fSIPXS
(
    readScalar(biokineticsProperties.subDict("stoichiometry").lookup("fSIPXS"))
);

Info << "    Reading mass fraction of soluble inerts SI from XS" << endl;
scalar fSIXS
(
    readScalar(biokineticsProperties.subDict("stoichiometry").lookup("fSIXS"))
);

Info << "    Reading mass fraction of solid inerts XI from XS" << endl;
scalar fXIXS
(
    readScalar(biokineticsProperties.subDict("stoichiometry").lookup("fXIXS"))
);

Info << "    Reading mass fraction of inorganic nitrogen in biomass" << endl;
scalar fNB
(
    readScalar(biokineticsProperties.subDict("stoichiometry").lookup("fNB"))
);

Info << "    Reading mass fraction of inorganic phosphorus in biomass" << endl;
scalar fPB
(
    readScalar(biokineticsProperties.subDict("stoichiometry").lookup("fPB"))
);

Info << "    Reading mass fraction of inorganic nitrogen "
     << "released at decay" << endl;
scalar fINDEC
(
    readScalar(biokineticsProperties.subDict("stoichiometry").lookup("fINDEC"))
);

Info << "    Reading mass fraction of inorganic phosphorus "
     << "released at decay" << endl;
scalar fIPDEC
(
    readScalar(biokineticsProperties.subDict("stoichiometry").lookup("fIPDEC"))
);

//=============================================================================

Info << "  Reading kinetic parameters" << endl;

Info << "    Reading maximum specific growth rate for "
     << "photoheterotrophic biomass on soluble organics SS, (muPH)" << endl;

dimensionedScalar muPH
(
    biokineticsProperties.subDict("kinetics").lookup("muPH")
);

Info << "    Reading maximum specific growth rate for "
     << "photoheterotrophic biomass on acetate SAC, (muAC)" << endl;
dimensionedScalar muAC
(
    biokineticsProperties.subDict("kinetics").lookup("muAC")
);

Info << "    Reading decay rate of XPB (kd)" << endl;
dimensionedScalar kd
(
    biokineticsProperties.subDict("kinetics").lookup("kd")
);

Info << "    Reading hydrolysis rate of XS (kh)" << endl;
dimensionedScalar kh
(
    biokineticsProperties.subDict("kinetics").lookup("kh")
);

Info << "    Reading half-saturation coefficient on soluble organics, KSS" << endl;

dimensionedScalar KSS
(
    biokineticsProperties.subDict("kinetics").lookup("KSS")
);

Info << "    Reading half-saturation coefficient on acetate, KSAC" << endl;

dimensionedScalar KSAC
(
    biokineticsProperties.subDict("kinetics").lookup("KSAC")
);

Info << "    Reading ammonium half-saturation coefficient on nitrogen uptake, KSIN" << endl;

dimensionedScalar KSIN
(
    biokineticsProperties.subDict("kinetics").lookup("KSIN")
);

Info << "    Reading free ammonia inhibition coefficient, KFA" << endl;

dimensionedScalar KFA
(
    biokineticsProperties.subDict("kinetics").lookup("KFA")
);

Info << "    Reading phosphate half-saturation coefficient KSIP" << endl;

dimensionedScalar KSIP
(
    biokineticsProperties.subDict("kinetics").lookup("KSIP")
);



Info << "    Reading half-saturation coefficient, radiative field KG" << endl;

dimensionedScalar KG
(
    biokineticsProperties.subDict("kinetics").lookup("KG")
);

//=============================================================================
Info << "  Reading species diffusivities" << endl;

dimensionedScalar DSSValue
(
    biokineticsProperties.subDict("diffusivities").lookup("DSS")
);

dimensionedScalar DXSValue
(
    biokineticsProperties.subDict("diffusivities").lookup("DXS")
);

dimensionedScalar DXPBValue
(
    biokineticsProperties.subDict("diffusivities").lookup("DXPB")
);

dimensionedScalar DXIValue
(
    biokineticsProperties.subDict("diffusivities").lookup("DXI")
);

dimensionedScalar DSIValue
(
    biokineticsProperties.subDict("diffusivities").lookup("DSI")
);

dimensionedScalar DSACValue
(
    biokineticsProperties.subDict("diffusivities").lookup("DSAC")
);

dimensionedScalar DSINValue
(
    biokineticsProperties.subDict("diffusivities").lookup("DSIN")
);


dimensionedScalar DSIPValue
(
    biokineticsProperties.subDict("diffusivities").lookup("DSIP")
);

Info << "  Reading turbulent Schmidt number, ScT" << endl;
scalar ScT
(
    readScalar(biokineticsProperties.lookup("ScT"))
);

Info << "  Reading liquid phase name, liquidName" << endl;
word liquidName
(
    biokineticsProperties.lookup("liquidName")
);

Info << "  Reading gas phase name, gasName" << endl;
word gasName
(
    biokineticsProperties.lookup("gasName")
);

Info << endl;
