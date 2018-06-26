// Scalar transport equation for readily biodegradable substrate

fvScalarMatrix SSEqn
(
    fvm::ddt(alphaLiquid, SS)
    + fvm::div(liquidPhase.alphaPhi(), SS)
    - fvm::laplacian(alphaLiquid*DSS, SS)
    ==
    - alphaLiquid*rhoHeterotrophAerobicGrowth/YH
    - alphaLiquid*rhoHeterotrophAnoxicGrowth/YH
    + alphaLiquid*rhoOrganicHydrolysis
    - fvm::Sp(dissipationCoeff, SS)
);

SSEqn.relax();
SSEqn.solve();
