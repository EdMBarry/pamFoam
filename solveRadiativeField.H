double radInterval
(
    runTime.controlDict().lookupOrDefault("m", 900.0)
);

double tol
(
    runTime.controlDict().lookupOrDefault("radTol", 0.1)
);

double timeValue = runTime.value() + VSMALL;

Info<< "Time Value is " << timeValue << nl << endl;

Info << "The modulus of time and m is " << Foam::mag(fmod(timeValue, radInterval)) << nl << endl;

if (Foam::mag(fmod(timeValue, radInterval)) <= tol)
{
    Info<< "Solving radiative field" << nl << endl;
    photoBioModel->correct();
}
