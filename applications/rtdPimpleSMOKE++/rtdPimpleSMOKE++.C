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
    rtdPimpleSMOKE++

Description
    Solves transient transport equation for a passive tracer to estimate
    the residence time distribution (RTD) function

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "fvOptions.H"
#include "pimpleControl.H"
#include "buildGlobalBoundaryList.H"
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
    #include "postProcess.H"

    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "createControl.H"
    #include "createTimeControls.H"
    
    const wordList globalBoundaryList = buildGlobalBoundaryList(mesh);

    #include "createFields.H"
    #include "createRTD.H"

    Info<< "\nSolving tracer transport equation\n" << endl;

    #include "compressibleCourantNo.H"
    #include "setInitialDeltaT.H"

    Info<< nl << "Starting time loop" << nl << endl;

    while (pimple.run(runTime))
    {
	// Set the time step 
	#include "readTimeControls.H"
	#include "compressibleCourantNo.H"
	#include "setDeltaT.H"

	// Advance in time
	runTime++;
	Info<< "Time = " << runTime.timeName() << nl << endl;

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

	// Post-processing
	#include "analyzeOutletPatches.H"

        runTime.write();

	// Write CPU time on the screen
	Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s" << "  ClockTime = " << runTime.elapsedClockTime() << " s" << nl << endl;
    }

    
    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
