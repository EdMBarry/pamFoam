fvScalarMatrix SNHEqn
(
    fvm::ddt(max(alphaLiquid, 1.0e-6), SNH)
    + fvm::div(liquidPhase.alphaPhi(), SNH)
    - fvm::laplacian(alphaLiquid*turbulence->nut()/ScT, SNH)
    ==
    - alphaLiquid*iXB*rhoHeterotrophAerobicGrowth
    - alphaLiquid*iXB*rhoHeterotrophAnoxicGrowth
    - alphaLiquid*(iXB+1/YA)*rhoAutotrophAerobicGrowth
    + alphaLiquid*rhoAmmonification
);
