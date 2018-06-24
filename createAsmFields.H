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

Info << endl;
