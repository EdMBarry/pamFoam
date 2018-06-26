// Scalar transport equation for oxygen

fvScalarMatrix SOEqn
(
    fvm::ddt(alphaLiquid, SO)
    + fvm::div(liquidPhase.alphaPhi(), SO)
    - fvm::laplacian(alphaLiquid*DSO, SO)
    ==
    - alphaLiquid*(1-YH)/YH*rhoHeterotrophAerobicGrowth
    - alphaLiquid*(4.57-YA)/YA*rhoAutotrophAerobicGrowth
    + alphaLiquid*oxygenTransfer->kLa()*(SOsat - SO)
    - fvm::Sp(dissipationCoeff, SO)
);

SOEqn.relax();
SOEqn.solve();
