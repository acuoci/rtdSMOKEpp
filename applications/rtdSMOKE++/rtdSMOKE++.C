/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application
    rtdSMOKE++

Description
    Solves transient transport equation for a passive tracer to estimate
    the residence time distribution (RTD) function

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "fvOptions.H"
#include "simpleControl.H"
#include <vector>
#include <iomanip>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

double Moment0(const std::vector<double>& t, const std::vector<double>& E)
{
	double m0 = 0.;
	for (unsigned int i=0;i<t.size()-1;i++)
		m0 += 0.5*(E[i+1]+E[i])*(t[i+1]-t[i]);
	return m0;
}

double Moment1(const std::vector<double>& t, const std::vector<double>& E)
{
	double m1 = 0.;
	for (unsigned int i=0;i<t.size()-1;i++)
		m1 += 0.5*(E[i+1]*t[i+1]+E[i]*t[i])*(t[i+1]-t[i]);
	return m1;
}

double Moment2(const std::vector<double>& t, const std::vector<double>& E)
{
	double m2 = 0.;
	for (unsigned int i=0;i<t.size()-1;i++)
		m2 += 0.5*(E[i+1]*t[i+1]*t[i+1]+E[i]*t[i]*t[i])*(t[i+1]-t[i]);
	return m2;
}

int main(int argc, char *argv[])
{
    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createMesh.H"

    simpleControl simple(mesh);

    #include "createFields.H"
    #include "createRTD.H"

    Info<< "\nSolving tracer transport equation\n" << endl;

    #include "CourantNo.H"

    while (simple.loop(runTime))
    {
        Info<< "Time = " << runTime.timeName() << nl << endl;

        while (simple.correctNonOrthogonal())
        {
            fvScalarMatrix FEqn
            (
                fvm::ddt(rho, F)
              + fvm::div(phi, F)
              - fvm::laplacian(rho*alpha, F)
             ==
                fvOptions(rho, F)
            );

            FEqn.relax();
            fvOptions.constrain(FEqn);
            FEqn.solve();
            fvOptions.correct(F);
        }

        runTime.write();

	// Post-processing
	#include "analyzeOutletPatches.H"
    }

    
    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
