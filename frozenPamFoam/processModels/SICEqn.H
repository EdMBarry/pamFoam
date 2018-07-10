// Scalar transport equation for inorganic nitrogen

fvScalarMatrix SICEqn
(
    fvm::ddt(alphaLiquid, SIC)
    + fvm::div(liquidPhase.alphaPhi(), SIC)
    - fvm::laplacian(alphaLiquid*DSIC, SIC)
    ==
    - fvm::Sp(dissipationCoeff, SIC)
    + alphaLiquid*fSICAC*rhoAcetateUptake
    + alphaLiquid*fSICSS*rhoPhotoheterotrophicUptake
    + alphaLiquid*fSICXS*rhoHydrolysis
    - alphaLiquid*fSICAU*rhoAutotrophicUptake
    + alphaLiquid*fICDEC*rhoDecay
);

SICEqn.relax();
SICEqn.solve();
