// Scalar transport equation for ammonia and ammonium nitrogen

fvScalarMatrix SNHEqn
(
    fvm::ddt(alphaLiquid, SNH)
    + fvm::div(liquidPhase.alphaPhi(), SNH)
    - fvm::laplacian(alphaLiquid*DSNH, SNH)
    ==
    - alphaLiquid*iXB*rhoHeterotrophAerobicGrowth
    - alphaLiquid*iXB*rhoHeterotrophAnoxicGrowth
    - alphaLiquid*(iXB+1/YA)*rhoAutotrophAerobicGrowth
    + alphaLiquid*rhoAmmonification
    - fvm::Sp(dissipationCoeff, SNH)
);

SNHEqn.relax();
SNHEqn.solve();
