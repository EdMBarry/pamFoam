// Scalar transport equation for hydrogen substrate

fvScalarMatrix SH2Eqn
(
    fvm::ddt(alphaLiquid, SH2)
    + fvm::div(liquidPhase.alphaPhi(), SH2)
    - fvm::laplacian(alphaLiquid*DSH2, SH2)
    ==
    - fvm::Sp(dissipationCoeff, SH2)
    + alphaLiquid * fSH2XS * rhoHydrolysis
    + (1 - YPBCH) * fSH2CH * rhoChemoheterotrophicUptake
    - alphaLiquid * fSH2AU * rhoAutotrophicUptake
);

SH2Eqn.relax();
SH2Eqn.solve();
SH2 *= pos(SH2);