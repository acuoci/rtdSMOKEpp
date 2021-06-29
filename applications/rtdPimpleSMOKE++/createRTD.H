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
std::vector<double> history_time;
std::vector<double> history_F;
std::vector<double> history_E;

history_time.push_back(0.);
history_F.push_back(0.);
history_E.push_back(0.);

// Open output file
autoPtr<std::ofstream> fRTD;
fRTD.reset(new std::ofstream("RTD.out", std::ios::out));
fRTD().setf(std::ios::scientific);

fRTD() << std::left << std::setw(16) << "time[s]";
fRTD() << std::left << std::setw(16) << "F[-]";
fRTD() << std::left << std::setw(16) << "E[1/s]";
fRTD() << std::left << std::setw(16) << "area[-]";
fRTD() << std::left << std::setw(16) << "tmean[s]";
fRTD() << std::left << std::setw(16) << "sigma2[s2]";
fRTD() << std::left << std::setw(16) << "Pe[-]";
fRTD() << std::endl;

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
	IOobject header
	(
		"mix:alpha",
		mesh.time().timeName(),
		mesh,
		IOobject::NO_READ
	);

	// Check if field exists and can be read
	if (header.typeHeaderOk<volScalarField>(true))
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
	else
	{

		const volScalarField Cp
		(
    			IOobject
    			(
        			"mix:Cp",
        			runTime.timeName(),
       		 		mesh,
       	 			IOobject::MUST_READ,
        			IOobject::AUTO_WRITE
    			),
    			mesh
		);

		const volScalarField lambda
		(
    			IOobject
    			(
        			"mix:lambda",
        			runTime.timeName(),
       		 		mesh,
       	 			IOobject::MUST_READ,
        			IOobject::AUTO_WRITE
    			),
    			mesh
		);

		alpha = lambda/rho/Cp;
	}
	
}

