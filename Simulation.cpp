/*
Program:     ContaminationFlow
Description: Monte Carlo simulator for satellite contanimation studies
Authors:     Rudolf Schönmann / Hoai My Van
Copyright:   TU Munich
Forked from: Molflow (CERN) (https://cern.ch/molflow)

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

Full license text: https://www.gnu.org/licenses/old-licenses/gpl-2.0.en.html
*/
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
