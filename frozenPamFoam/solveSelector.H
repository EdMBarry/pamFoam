bool solvePam = false;
if (biokineticsProperties.found("solvePam"))
{
    solvePam = Switch(biokineticsProperties.lookup("solvePam"));
}

bool solveFlow = false;
if (biokineticsProperties.found("solveFlow"))
{
    solveFlow = Switch(biokineticsProperties.lookup("solveFlow"));
}

bool solveRad = false;
if (biokineticsProperties.found("solveRad"))
{
    solveRad = Switch(biokineticsProperties.lookup("solveRad"));
}
