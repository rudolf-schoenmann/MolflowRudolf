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

	//memset(&tmpGlobalResult, 0, sizeof(GlobalHitBuffer));
	tmpGlobalResult.globalHits.nbMCHit = 0;
	tmpGlobalResult.globalHits.nbHitEquiv = 0.0;
	tmpGlobalResult.globalHits.nbAbsEquiv = 0.0;
	tmpGlobalResult.globalHits.nbDesorbed = 0;
	tmpGlobalResult.globalHits.sum_1_per_ort_velocity = 0;
	tmpGlobalResult.globalHits.sum_1_per_velocity = 0;
	tmpGlobalResult.globalHits.sum_v_ort = 0;
	tmpGlobalResult.globalHits.covering = boost::multiprecision::uint128_t(0);
	tmpGlobalResult.hitCacheSize = 0;
	tmpGlobalResult.lastHitIndex = 0;
	for (size_t i = 0; i < HITCACHESIZE; i++) {
		tmpGlobalResult.hitCache[i].pos.x = 0.0;
		tmpGlobalResult.hitCache[i].pos.y = 0.0;
		tmpGlobalResult.hitCache[i].pos.z = 0.0;
		tmpGlobalResult.hitCache[i].type = 0;
	}
	for (size_t i = 0; i < LEAKCACHESIZE; i++) {
		tmpGlobalResult.leakCache[i].dir.x = 0.0;
		tmpGlobalResult.leakCache[i].dir.y = 0.0;
		tmpGlobalResult.leakCache[i].dir.z = 0.0;
		tmpGlobalResult.leakCache[i].pos.x = 0.0;
		tmpGlobalResult.leakCache[i].pos.y = 0.0;
		tmpGlobalResult.leakCache[i].pos.z = 0.0;
	}
	tmpGlobalResult.lastLeakIndex = 0;
	tmpGlobalResult.leakCacheSize = 0;
	tmpGlobalResult.nbLeakTotal = 0;
	for (size_t i = 0; i < 3; i++) {
		tmpGlobalResult.texture_limits[i].max.all = 0.0;
		tmpGlobalResult.texture_limits[i].max.moments_only = 0.0;
		tmpGlobalResult.texture_limits[i].min.all = 0.0;
		tmpGlobalResult.texture_limits[i].min.moments_only = 0.0;
	}
	tmpGlobalResult.distTraveled_total = 0.0;
	tmpGlobalResult.distTraveledTotal_fullHitsOnly = 0.0;

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
