Info << endl;
Info << "Creating oxygen transfer model" << endl;

const dictionary& oxygenTransferModelDict = biokineticsProperties.subDict("oxygenTransferModel");
autoPtr<oxygenTransferModel> oxygenTransfer 
(
    oxygenTransferModel::New(oxygenTransferModelDict, liquidPhase, gasPhase)
);

Info << endl;
