// Scalar transport equation for inorganic nitrogen

fvScalarMatrix SICEqn
(
    fvm::ddt(alphaLiquid, SIC)
    + fvm::div(liquidPhase.alphaPhi(), SIC)
    - fvm::laplacian(alphaLiquid*SIC, SIC)
    ==
    - alphaLiquid*fNB*YPBPH*rhoAcetateUptake
    - alphaLiquid*fNB*YPBPH*rhoPhotoheterotrophicUptake
    + alphaLiquid*fSINXS*rhoHydrolysis
    + alphaLiquid*fINDEC*rhoDecay
    - fvm::Sp(dissipationCoeff, SIC)
);

SICEqn.relax();
SICEqn.solve();