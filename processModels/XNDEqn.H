// Scalar transport equation for particulate biodegradable organic nitrogen

fvScalarMatrix XNDEqn
(
    fvm::ddt(alphaLiquid, XND)
    + fvm::div(liquidPhase.alphaPhi(), XND)
    - fvm::laplacian(alphaLiquid*DXND, XND)
    ==
      alphaLiquid*(iXB-fP*iXP)*rhoHeterotrophDecay
    + alphaLiquid*(iXB-fP*iXP)*rhoAutotrophDecay
    - alphaLiquid*rhoOrganicNitrogenHydrolysis
    - fvm::Sp(dissipationCoeff, XND)
);

XNDEqn.relax();
XNDEqn.solve();
