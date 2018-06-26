Info << "Reading scalar fields for ASM" << endl;

Info << "  Reading field SS (Readily biodegradable substrate)" << endl;

volScalarField SS
(
    IOobject
    (
        "SS",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info << "  Reading field XS (slowly biodegradable substrate)" << endl;

volScalarField XS
(
    IOobject
    (
        "XS",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info << "  Reading field XBH (active heterotrophic biomass)" << endl;

volScalarField XBH
(
    IOobject
    (
        "XBH",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info << "  Reading field XBA (active autotrophic biomass)" << endl;

volScalarField XBA
(
    IOobject
    (
        "XBA",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info << "  Reading field XP (particulate products arising from biomass decay)" << endl;

volScalarField XP
(
    IOobject
    (
        "XP",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info << "  Reading field SO (oxygen)" << endl;

volScalarField SO
(
    IOobject
    (
        "SO",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info << "  Reading field SNO (nitrate and nitrite)" << endl;

volScalarField SNO
(
    IOobject
    (
        "SNO",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info << "  Reading field SNH (ammonia and ammonium nitrogen)" << endl;

volScalarField SNH
(
    IOobject
    (
        "SNH",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info << "  Reading field SND (soluble biodegradable organic nitrogen)" << endl;

volScalarField SND
(
    IOobject
    (
        "SND",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info << "  Reading field XND (particulate biodegradable organic nitrogen)" << endl;

volScalarField XND
(
    IOobject
    (
        "XND",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info << "Creating scalar fields for process rates" << endl;

volScalarField rhoHeterotrophAerobicGrowth
(
    IOobject
    (
        "rhoHeterotrophAerobicGrowth",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("rhoHeterotrophAerobicGrowth", dimMass/dimVolume/dimTime, 0)
);

volScalarField rhoHeterotrophAnoxicGrowth
(
    IOobject
    (
        "rhoHeterotrophAnoxicGrowth",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("rhoHeterotrophAnoxicGrowth", dimMass/dimVolume/dimTime, 0)
);

volScalarField rhoAutotrophAerobicGrowth
(
    IOobject
    (
        "rhoAutotrophAerobicGrowth",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("rhoAutotrophAerobicGrowth", dimMass/dimVolume/dimTime, 0)
);

volScalarField rhoHeterotrophDecay
(
    IOobject
    (
        "rhoHeterotrophDecay",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("rhoHeterotrophDecay", dimMass/dimVolume/dimTime, 0)
);

volScalarField rhoAutotrophDecay
(
    IOobject
    (
        "rhoAutotrophDecay",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("rhoAutotrophDecay", dimMass/dimVolume/dimTime, 0)
);

volScalarField rhoAmmonification
(
    IOobject
    (
        "rhoAmmonification",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("rhoAmmonification", dimMass/dimVolume/dimTime, 0)
);

volScalarField rhoOrganicHydrolysis
(
    IOobject
    (
        "rhoOrganicHydrolysis",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("rhoOrganicHydrolysis", dimMass/dimVolume/dimTime, 0)
);

volScalarField rhoOrganicNitrogenHydrolysis
(
    IOobject
    (
        "rhoOrganicNitrogenHydrolysis",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("rhoOrganicNitrogenHydrolysis", dimMass/dimVolume/dimTime, 0)
);

Info << "Creating scalar fields for species diffusivities" << endl;

volScalarField DSND
(
    IOobject
    (
        "DSND",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    DSNDValue
);

volScalarField DSNH
(
    IOobject
    (
        "DSNH",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    DSNHValue
);

volScalarField DSNO
(
    IOobject
    (
        "DSNO",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    DSNOValue
);

volScalarField DSO
(
    IOobject
    (
        "DSO",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    DSOValue
);

volScalarField DSS
(
    IOobject
    (
        "DSS",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    DSSValue
);

volScalarField DXBA
(
    IOobject
    (
        "DXBA",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    DXBAValue
);

volScalarField DXBH
(
    IOobject
    (
        "DXBH",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    DXBHValue
);

volScalarField DXND
(
    IOobject
    (
        "DXND",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    DXNDValue
);

volScalarField DXP
(
    IOobject
    (
        "DXP",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    DXPValue
);

volScalarField DXS
(
    IOobject
    (
        "DXS",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    DXSValue
);

Info << "Creating scalar field for dissipation function" << endl;

volScalarField dissipationCoeff
(
    IOobject
    (
        "dissipationCoeff",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("dissipationCoeff", dimless/dimTime, 0)
);

Info << endl;
