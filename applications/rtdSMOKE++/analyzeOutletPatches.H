{
	// Inlet properties
	List<scalar> mIn(listInletPatches.size());
	List<scalar> FIn(listInletPatches.size());
	for (label i = 0; i<listInletPatches.size();i++)
	{
		const label patch = mesh.boundaryMesh().findPatchID(listInletPatches[i]); 
		mIn[i] = -1*sum(phi.boundaryField()[patch]); 
		FIn[i] = -1*sum(phi.boundaryField()[patch]*F.boundaryField()[patch])/mIn[i];
   	}

	// Outlet properties
	List<scalar> mOut(listOutletPatches.size());
	List<scalar> FOut(listOutletPatches.size());
	for (label i = 0; i<listOutletPatches.size();i++)
	{
		const label patch = mesh.boundaryMesh().findPatchID(listOutletPatches[i]); 
		mOut[i] = sum(phi.boundaryField()[patch]); 
		FOut[i] = sum(phi.boundaryField()[patch]*F.boundaryField()[patch])/mOut[i];
   	}

	// Storing old values
	const scalar FOutMeanOld = FOutMean;
	const scalar deltat = runTime.deltaT().value();

	// Average inlet properties
	mInTot = 0.;
	FInMean = 0.;
	for (label i = 0; i<listInletPatches.size();i++)
	{
		mInTot += mIn[i];
		FInMean += mIn[i]*FIn[i];
	}
	FInMean /= mInTot;

	// Average outlet properties
	mOutTot = 0.;
	FOutMean = 0.;
	for (label i = 0; i<listOutletPatches.size();i++)
	{
		mOutTot += mOut[i];
		FOutMean += mOut[i]*FOut[i];
	}
	FOutMean /= mOutTot;

	// Summary on the screen
	Info << endl;
	Info << "INLET: " << " m[kg/s]=" << mInTot << " FMean=" << FInMean << endl;
	for (label i = 0; i<listInletPatches.size();i++)
		Info << " - " << listInletPatches[i] << ": m[kg/s]=" << mIn[i] << " FMean=" << FIn[i] << endl;
	Info << "OUTLET: " << " m[kg/s]=" << mOutTot << " FMean=" << FOutMean << endl;
	for (label i = 0; i<listOutletPatches.size();i++)
		Info << " - " << listOutletPatches[i] << ": m[kg/s]=" << mOut[i] << " FMean=" << FOut[i] << endl;
	Info << endl;

	// Print on file
	if (deltat > 1e-10)
	{
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

		fRTD << std::left << std::setw(16) << tf;
		fRTD << std::left << std::setw(16) << FOutMean;
		fRTD << std::left << std::setw(16) << E;
		fRTD << std::left << std::setw(16) << area;
		fRTD << std::left << std::setw(16) << tMean;
		fRTD << std::left << std::setw(16) << sigma2;
		fRTD << std::left << std::setw(16) << Pe;
		fRTD << std::endl;
	}
}


