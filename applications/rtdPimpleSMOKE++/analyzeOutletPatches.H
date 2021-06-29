{
	const scalar deltat = runTime.deltaT().value();

	if (deltat > 1e-10)
	{
		const scalar FOutMeanOld = history_F[history_F.size()-1];

		scalar mInTot = 0.;
		scalar FIn = 0.;
		scalar mOutTot = 0.;
		scalar FOut = 0.;

		forAll (globalBoundaryList, patchi)
		{
			for(int iP=0;iP<listInletPatches.size();iP++)
			{
				scalar flow = 0.;	
				scalar F_flow = 0.;

				if(mesh.boundary()[patchi].name() == listInletPatches[iP])
				{
					flow =
					sum
					(
						(
						rho.boundaryField()[patchi]
						* U.boundaryField()[patchi]
						& mesh.boundary()[patchi].Sf()
						)
					);

					F_flow =
					sum
					(
						(
						rho.boundaryField()[patchi] * F.boundaryField()[patchi]
						* U.boundaryField()[patchi]
						& mesh.boundary()[patchi].Sf()
						)
					);
				}

				reduce(flow, sumOp<scalar>());
				reduce(F_flow, sumOp<scalar>());

				mInTot += mag(flow);
				FIn += mag(F_flow);
			}

			for(int oP=0;oP<listOutletPatches.size();oP++)
			{
				scalar flow = 0.;	
				scalar F_flow = 0.;

				if(mesh.boundary()[patchi].name() == listOutletPatches[oP])
				{
					flow =
					sum
					(
						(
						rho.boundaryField()[patchi]
						* U.boundaryField()[patchi]
						& mesh.boundary()[patchi].Sf()
						)
					);

					F_flow =
					sum
					(
						(
						rho.boundaryField()[patchi] * F.boundaryField()[patchi]
						* U.boundaryField()[patchi]
						& mesh.boundary()[patchi].Sf()
						)
					);
				}

				reduce(flow, sumOp<scalar>());
				reduce(F_flow, sumOp<scalar>());

				mOutTot += mag(flow);
				FOut += mag(F_flow);
			}
		}
	
		// Print on file
		if (Pstream::master())
		{
			scalar FInMean = 0.;
			if (mInTot > 1e-32) FInMean = FIn/mInTot;
	
			scalar FOutMean = 0.;
			if (mOutTot > 1e-32) FOutMean = FOut/mOutTot;

			// Summary on the screen
			Info << endl;
			Info << "INLET: " << " m[kg/s]=" << mInTot << " FMean=" << FInMean << endl;
			Info << "OUTLET: " << " m[kg/s]=" << mOutTot << " FMean=" << FOutMean << endl;
			Info << endl;

			const double tf = runTime.value();
			const double E = (FOutMean-FOutMeanOld)/deltat;

			history_time.push_back(tf);
			history_F.push_back(FOutMean);
			history_E.push_back(E);

			const double area = Moment0(history_time, history_E);

			const double tMean = Moment1(history_time, history_E);
			const double sigma2 = Moment2(history_time, history_E) - tMean*tMean;
			
			double Pe = 0.;
			if (sigma2 > 1e-12)
			{
				Pe = 2.*tMean*tMean/sigma2;
				for (unsigned int k=0;k<5;k++)
					Pe = 2./(sigma2/tMean/tMean + 2./Pe/Pe*(1.-std::exp(-Pe)));
			}

			fRTD() << std::left << std::setw(16) << tf;
			fRTD() << std::left << std::setw(16) << FOutMean;
			fRTD() << std::left << std::setw(16) << E;
			fRTD() << std::left << std::setw(16) << area;
			fRTD() << std::left << std::setw(16) << tMean;
			fRTD() << std::left << std::setw(16) << sigma2;
			fRTD() << std::left << std::setw(16) << Pe;
			fRTD() << std::endl;

		}

	}
}


