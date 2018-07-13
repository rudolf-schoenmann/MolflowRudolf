#include "Simulation.h"
#include "IntersectAABB_shared.h"

SuperStructure::SuperStructure()
{
	aabbTree = NULL;
}

SuperStructure::~SuperStructure()
{
	SAFE_DELETE(aabbTree);
}

Simulation::Simulation()
{
	totalDesorbed = 0;

	loadOK = false;
	wp.sMode = MC_MODE;
	currentParticle.lastHitFacet = NULL;

	hasVolatile = false;

	memset(&tmpGlobalResult, 0, sizeof(GlobalHitBuffer));

	sh.nbSuper = 0;
	acDensity =
		acMatrix =
		acDensity =
		acDesorb =
		acAbsorb =
		acTArea =
		acRho = 
		acTMatrix =
		acTDensity = acArea =
		NULL;
		
		acLines =
		acTLines = NULL;
}
