# asmFoam - A multiphase CFD solver with coupled activated sludge biokinetics

`asmFoam` is a multiphase CFD solver that is coupled with an activated sludge model, based on the ASM1 model of Henze et al. (1987). There are currently two oxygen transfer models available, based on: (i) constant, spatially uniform values of `kL` (mass transfer coefficient) and `a` (interfacial area), expressedd as the product, kLa  (ii) dynamically calculated values of `kL`, based on Higbie's penetration theory, and `a`, based on local bubble size and gas volume fraction.

The model is implemented using the [OpenFOAM](https://openfoam.org) libaries for discretization and solution of the governing equations.

## Getting Started

### Prerequisites

Install OpenFOAM and set up paths.

Install ParaView.

### Installation

Build `asmFoam`

## Running the Code

Setting up a case for the `asmFoam` solver proceeds similarly to typical `OpenFOAM` cases, with the exception of the additional configuration required by those models implemented specifically for `asmFoam`. It is assumed that the user is generally familiar with using `OpenFOAM` and setting up cases, so a detailed explanation will not be provided here.

### Initial and Boundary Conditions

To use the `asmFoam` solver, boundary and initial conditions must be specified on all additional fields that are solved, using appropriate files in the `0` subdirectory of the main case directory. The additional fields introduced by the `asmFoam` solver are: `SND`, `SNH`, `SNO`, `SO`, `SS`, `XBA`, `XBH`, `XND`, `XP`, and `XS`. The units of all these fields are kilograms per cubic metre.

### ASM Configuration (biokineticsProperties dictionary)

A dictionary file called `biokineticsProperties` must be provided in the `constant` subdirectory to configure the activated sludge model. This file requires a number of dictionaries and parameters to be defined, as outlined below.

The stoichiometry coefficients are defined as follows:

```
stoichiometry
{
    YA     0.24;
    YH     0.67;
    fP     0.08;
    iXB    0.086;
    iXP    0.06;
}

```

A description of each parameter is provided in the table below. Note that all of the stoichometric coefficients are dimensionless.

| Parameter | Description |
| --------- | ----------- |
| YA        | Autotrophic yield |
| YH        | Heterotrophic yield |
| fP        | Fraction of biomass yielding particulate products |
| iXB       | Reading mass of nitrogen per mass of COD in biomass |
| iXP       | Reading mass of nitrogen per mass of COD in products from biomass |

### No Flow Verification Problem

All test cases are available in the [asmFoamTests](https://bitbucket.org/cdegroot/asmfoamtests/) Git repository. The most basic example is called ``NoFlow`` and can be found in the ``Verification`` directory of the test repository. As the name suggests, there is no flow involved in this problem, which means that the correct numerical solution to the problem can be found by simply solving the ordinary differential equations representing the ASM1 model.

## Authors

- [Chris DeGroot](http://www.eng.uwo.ca/mechanical/faculty/degroot_c/index.html)

## License

- This project is licensed under the GNU General Public License, version 3 (GPL-3.0). See LICENSE.txt for details.

## Acknowledgements

- This work is supported by a research grant from the Natural Sciences and Engineering Research Council of Canada (NSERC) through the Discovery Grants program.
