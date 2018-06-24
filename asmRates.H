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
    muH * SS/(KS + SS) * SO/(KOH + SO) * XBH
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
    muH * SS/(KS + SS) * KOH/(KOH + SO) * SNO/(KNO + SNO) * etag * XBH
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
    muA * SNH/(KNH + SNH) * SO/(KOA + SO) * XBA
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
    bH * XBH
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
    bA * XBA
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
    ka * SND * XBH
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
    kh * XS/XBH/(KX + XS/XBH) * (SO/(KOH + SO) + etah * KOH/(KOH + SO) * SNO/(KNO + SNO)) * XBH
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
    rhoOrganicHydrolysis * XND / XS
);
