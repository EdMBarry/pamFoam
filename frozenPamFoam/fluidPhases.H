phaseModel& liquidPhase = fluid.phases()[liquidName];
phaseModel& gasPhase = fluid.phases()[gasName];
const volScalarField& alphaLiquid = liquidPhase;
const volScalarField& alphaGas = gasPhase;
