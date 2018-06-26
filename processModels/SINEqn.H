// Scalar transport equation for inorganic nitrogen

fvScalarMatrix SINEqn
(
    fvm::ddt(alphaLiquid, SIN)
    + fvm::div(liquidPhase.alphaPhi(), SIN)
    - fvm::laplacian(alphaLiquid*SIN, SIN)
    ==
    + alphaLiquid*fINXS * rhoHydrolysis
    - alphaLiquid*fNB*YPBPH*rhoAcetateUptake
    - alphaLiquid*fNB*YPBPH*rhoPhotoheterotrophicUptake
    + alphaLiquid*fINDEC*rhoDecay
    - fvm::Sp(dissipationCoeff, SIN)
);

SINEqn.relax();
SINEqn.solve();