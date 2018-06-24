fvScalarMatrix SNDEqn
(
    fvm::ddt(max(alphaLiquid, 1.0e-6), SND)
    + fvm::div(liquidPhase.alphaPhi(), SND)
    - fvm::laplacian(alphaLiquid*turbulence->nut()/ScT, SND)
    ==
    - alphaLiquid*rhoAmmonification
    + alphaLiquid*rhoOrganicNitrogenHydrolysis
);
