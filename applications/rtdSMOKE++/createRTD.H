// Reading RTD Dictionary
Info<< "Reading RTD dictionary\n" << endl;
IOdictionary RTDDictionary
(
	IOobject
	(
		"RTD",
		U.time().constant(),
		U.db(),
		IOobject::MUST_READ,
		IOobject::NO_WRITE
	)
);

List<word>  listInletPatches(RTDDictionary.lookup("inletPatches"));
List<word>  listOutletPatches(RTDDictionary.lookup("outletPatches"));
Switch diffusivity_from_fields = Switch(RTDDictionary.lookup(word("diffusivityFromFields")));


// Check patches
for (label i = 0; i<listInletPatches.size();i++)
{
	label patch = mesh.boundaryMesh().findPatchID(listInletPatches[i]); 
    	if (patch==-1) 
		FatalError << "Failure to find patch named " << listInletPatches[i] << endl;
}
for (label i = 0; i<listOutletPatches.size();i++)
{
	label patch = mesh.boundaryMesh().findPatchID(listOutletPatches[i]); 
    	if (patch==-1) 
		FatalError << "Failure to find patch named " << listOutletPatches[i] << endl;
}		

// Summary on the screen
Info << endl;
Info << "List of inlet patches: " << endl;
for (label i = 0; i<listInletPatches.size();i++)
	Info << " - " << listInletPatches[i] << endl;
Info << "List of outlet patches: " << endl;
for (label i = 0; i<listOutletPatches.size();i++)
	Info << " - " << listOutletPatches[i] << endl;
Info << endl;

// Create global variables
scalar mInTot = 0.;
scalar mOutTot = 0.;
scalar FInMean = 0.;
scalar FOutMean = 0.;
std::vector<double> history_time;
std::vector<double> history_F;
std::vector<double> history_E;

// Open output file
std::ofstream fRTD("RTD.out", std::ios::out);
fRTD.setf(std::ios::scientific);
fRTD << std::left << std::setw(16) << "time[s]";
fRTD << std::left << std::setw(16) << "F[-]";
fRTD << std::left << std::setw(16) << "E[1/s]";
fRTD << std::left << std::setw(16) << "area[-]";
fRTD << std::left << std::setw(16) << "tmean[s]";
fRTD << std::left << std::setw(16) << "sigma2[s2]";
fRTD << std::left << std::setw(16) << "Pe[-]";
fRTD << std::endl;

// Evaluation of passive scalar diffusion
if (diffusivity_from_fields == false)
{
	const double muA = -1.853769e+01;
	const double muB =  2.232971e+00;
	const double muC = -2.102014e-01;
	const double muD =  9.264622e-03;

	const scalarField& TCells = T.internalField();
	const scalarField& pCells = p.internalField();
	const scalarField& rhoCells = rho.internalField();
	scalarField& alphaCells = alpha.ref(); 

	forAll(TCells, celli)
	{
		const double logT = std::log(TCells[celli]);
		const double MW = rhoCells[celli]*8314.*TCells[celli]/pCells[celli];
		const double Cc = 28.014/MW;
		alphaCells[celli] = Cc * std::exp( muA+muB*logT+muC*(logT*logT)+muD*(logT*logT*logT) ) / rhoCells[celli];
	}

	const volScalarField::Boundary& Tf = T.boundaryField();
	const volScalarField::Boundary& pf = p.boundaryField();
	const volScalarField::Boundary& rhof = rho.boundaryField();
	volScalarField::Boundary& alphaf = alpha.boundaryFieldRef();

	forAll(T.boundaryField(), patchi)
	{
		const fvPatchScalarField& pT = Tf[patchi];
		const fvPatchScalarField& pp = pf[patchi];
		const fvPatchScalarField& prho = rhof[patchi];
		fvPatchScalarField& palpha = alphaf[patchi];

		forAll(pT, facei)
		{
			const double logT = std::log(pT[facei]);
			const double MW = prho[facei]*8314.*pT[facei]/pp[facei];
			const double Cc = 28.014/MW;
			palpha[facei] = Cc * std::exp( muA+muB*logT+muC*(logT*logT)+muD*(logT*logT*logT) ) / prho[facei];
		}
	}
}
else
{
	const volScalarField alpha_
	(
    		IOobject
    		(
        		"mix:alpha",
        		runTime.timeName(),
       		 	mesh,
       	 		IOobject::MUST_READ,
        		IOobject::AUTO_WRITE
    		),
    		mesh
	);

	alpha = alpha_;
}

