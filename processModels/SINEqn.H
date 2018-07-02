// Scalar transport equation for inorganic nitrogen

fvScalarMatrix SINEqn
(
    fvm::ddt(alphaLiquid, SIN)
    + fvm::div(liquidPhase.alphaPhi(), SIN)
    - fvm::laplacian(alphaLiquid*DSIN, SIN)
    ==
    - alphaLiquid*fNB*YPBPH*rhoAcetateUptake
    - alphaLiquid*fNB*YPBPH*rhoPhotoheterotrophicUptake
    - alphaLiquid*fNB*YPBCH*rhoChemoheterotrophicUptake
    - alphaLiquid*fNB*YPBAU*rhoAutotrophicUptake
    + alphaLiquid*fSINXS*rhoHydrolysis
    + alphaLiquid*fINDEC*rhoDecay
    - fvm::Sp(dissipationCoeff, SIN)
);

SINEqn.relax();
SINEqn.solve();