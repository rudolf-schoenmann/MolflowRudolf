/*
Program:     MolFlow+ / Synrad+
Description: Monte Carlo simulator for ultra-high vacuum and synchrotron radiation
Authors:     Jean-Luc PONS / Roberto KERSEVAN / Marton ADY
Copyright:   E.S.R.F / CERN
Website:     https://cern.ch/molflow

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
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <sstream>
#include "Simulation.h"
#include "IntersectAABB_shared.h"
#include "Random.h"
#include "GLApp/MathTools.h"
#include <tuple> //std::tie

extern SIMULATION *sHandle;

// Compute area of all the desorption facet

void CalcTotalOutgassing() {
	int i, j, k;
	size_t tSize;
	SubprocessFacet *f;
	double scale;
	//float scale_precomputed;

	// Update texture increment for MC
	//scale_precomputed=(float)(40.0/(sqrt(8.0*8.31/(PI*sHandle->gasMass*0.001))));
	for (j = 0; j < sHandle->nbSuper; j++) {
		for (i = 0; i < sHandle->str[j].nbFacet; i++) {
			f = sHandle->str[j].facets[i];
			if (f->inc) {
				tSize = f->sh.texWidth*f->sh.texHeight;
				scale = 1.0;
				if (f->sh.is2sided) scale = scale * 0.5;
				for (k = 0; k < tSize; k++) f->inc[k] *= scale;
				f->fullSizeInc *= scale;
			}
		}
	}

}

//void PolarToCartesian(SubprocessFacet *iFacet, double theta, double phi, bool reverse) {
//
//	Vector3d U, V, N;
//	double u, v, n;
//
//	// Polar in (nU,nV,N) to Cartesian(x,y,z) transformation  ( nU = U/|U| , nV = V/|V| )
//	// tetha is the angle to the normal of the facet N, phi to U
//	// ! See Geometry::InitializeGeometry() for further informations on the (U,V,N) basis !
//	// (nU,nV,N) and (x,y,z) are both left handed
//
//	//This should be a speed-up routine, but I didn't experience any speed difference so I commented it out. Marton
//	/*#ifdef WIN
//	_asm {                    // FPU stack
//	fld qword ptr [theta]
//	fsincos                 // cos(t)        sin(t)
//	fld qword ptr [phi]
//	fsincos                 // cos(p)        sin(p) cos(t) sin(t)
//	fmul st(0),st(3)        // cos(p)*sin(t) sin(p) cos(t) sin(t)
//	fstp qword ptr [u]      // sin(p)        cos(t) sin(t)
//	fmul st(0),st(2)        // sin(p)*sin(t) cos(t) sin(t)
//	fstp qword ptr [v]      // cos(t) sin(t)
//	fstp qword ptr [n]      // sin(t)
//	fstp qword ptr [dummy]  // Flush the sin(t)
//	}
//	#else*/
//	u = sin(theta)*cos(phi);
//	v = sin(theta)*sin(phi);
//	n = cos(theta);
//	//#endif
//
//	// Get the (nU,nV,N) orthonormal basis of the facet
//	U = iFacet->sh.nU;
//	V = iFacet->sh.nV;
//	N = iFacet->sh.N;
//	if (reverse) {
//		N.x = N.x*(-1.0);
//		N.y = N.y*(-1.0);
//		N.z = N.z*(-1.0);
//	}
//
//	// Basis change (nU,nV,N) -> (x,y,z)
//	sHandle->pDir.x = u*U.x + v*V.x + n*N.x;
//	sHandle->pDir.y = u*U.y + v*V.y + n*N.y;
//	sHandle->pDir.z = u*U.z + v*V.z + n*N.z;
//
//}
//
//void CartesianToPolar(SubprocessFacet *iFacet, double *theta, double *phi) {
//
//	// Get polar coordinates of the incoming particule direction in the (U,V,N) facet space.
//	// Note: The facet is parallel to (U,V), we use its (nU,nV,N) orthonormal basis here.
//	// (nU,nV,N) and (x,y,z) are both left handed
//
//	// Cartesian(x,y,z) to polar in (nU,nV,N) transformation
//
//	// Basis change (x,y,z) -> (nU,nV,N)
//	// We use the fact that (nU,nV,N) belongs to SO(3)
//	double u = Dot(sHandle->pDir, iFacet->sh.nU);
//	double v = Dot(sHandle->pDir, iFacet->sh.nV);
//	double n = Dot(sHandle->pDir, iFacet->sh.N);
//
//	/*
//	// (u,v,n) -> (theta,phi)
//	double rho = sqrt(v*v + u*u);
//	*theta = acos(n);              // Angle to normal (PI/2 .. PI for frontal and 0..PI/2 for back incidence)
//	*phi = asin(v / rho);			// Returns -PI/2 ... +PI/2
//	if (u < 0.0) *phi = PI - *phi;  // Angle to U, -PI/2 .. 3PI/2
//	*/
//
//	*theta = acos(n);
//	*phi = atan2(v, u); // -PI..PI
//}

void UpdateMCHits(Dataport *dpHit,int prIdx, size_t nbMoments, DWORD timeout) {

	BYTE *buffer;
	GlobalHitBuffer *gHits;
	TEXTURE_MIN_MAX texture_limits_old[3];
	int i, j, s, x, y;
#ifdef _DEBUG
	double t0, t1;
	t0 = GetTick();
#endif
	SetState(NULL, "Waiting for 'hits' dataport access...", false, true);
	sHandle->lastHitUpdateOK = AccessDataportTimed(dpHit, timeout);
	SetState(NULL, "Updating MC hits...", false, true);
	if (!sHandle->lastHitUpdateOK) return; //Timeout, will try again later

	buffer = (BYTE*)dpHit->buff;
	gHits = (GlobalHitBuffer *)buffer;

	// Global hits and leaks: adding local hits to shared memory
	gHits->total.hit.nbMCHit += sHandle->tmpGlobalCount.hit.nbMCHit;
	gHits->total.hit.nbHitEquiv += sHandle->tmpGlobalCount.hit.nbHitEquiv;
	gHits->total.hit.nbAbsEquiv += sHandle->tmpGlobalCount.hit.nbAbsEquiv;
	gHits->total.hit.nbDesorbed += sHandle->tmpGlobalCount.hit.nbDesorbed;
	gHits->distTraveled_total += sHandle->distTraveledSinceUpdate_total;
	gHits->distTraveledTotal_fullHitsOnly += sHandle->distTraveledSinceUpdate_fullHitsOnly;

	//Memorize current limits, then do a min/max search
	for (i = 0; i < 3; i++) {
		texture_limits_old[i] = gHits->texture_limits[i];
		gHits->texture_limits[i].min.all = gHits->texture_limits[i].min.moments_only = HITMAX;
		gHits->texture_limits[i].max.all = gHits->texture_limits[i].max.moments_only = 0;
	}

	gHits->sMode = MC_MODE;
	//for(i=0;i<BOUNCEMAX;i++) gHits->wallHits[i] += sHandle->wallHits[i];

	// Leak
	for (size_t leakIndex = 0; leakIndex < sHandle->leakCacheSize; leakIndex++)
		gHits->leakCache[(leakIndex + gHits->lastLeakIndex) % LEAKCACHESIZE] = sHandle->leakCache[leakIndex];
	gHits->nbLeakTotal += sHandle->nbLeakSinceUpdate;
	gHits->lastLeakIndex = (gHits->lastLeakIndex + sHandle->leakCacheSize) % LEAKCACHESIZE;
	gHits->leakCacheSize = Min(LEAKCACHESIZE, gHits->leakCacheSize + sHandle->leakCacheSize);

	// HHit (Only prIdx 0)
	if (prIdx == 0) {
		for (size_t hitIndex = 0; hitIndex < sHandle->hitCacheSize; hitIndex++)
			gHits->hitCache[(hitIndex + gHits->lastHitIndex) % HITCACHESIZE] = sHandle->hitCache[hitIndex];

		if (sHandle->hitCacheSize > 0) {
			gHits->lastHitIndex = (gHits->lastHitIndex + sHandle->hitCacheSize) % HITCACHESIZE;

			//if (gHits->lastHitIndex < (HITCACHESIZE - 1)) {
				//gHits->lastHitIndex++;
				gHits->hitCache[gHits->lastHitIndex].type = HIT_LAST; //Penup (border between blocks of consecutive hits in the hit cache)
			//}

			gHits->hitCacheSize = Min(HITCACHESIZE, gHits->hitCacheSize + sHandle->hitCacheSize);
		}
	}

	size_t facetHitsSize = (1 + nbMoments) * sizeof(FacetHitBuffer);
	// Facets
	for (s = 0; s < sHandle->nbSuper; s++) {
		for (i = 0; i < sHandle->str[s].nbFacet; i++) {

			SubprocessFacet *f = sHandle->str[s].facets[i];
			if (f->hitted) {

				for (int m = 0; m < (1 + nbMoments); m++) {
					FacetHitBuffer *fFit = (FacetHitBuffer *)(buffer + f->sh.hitOffset + m * sizeof(FacetHitBuffer));
					fFit->hit.nbAbsEquiv += f->counter[m].hit.nbAbsEquiv;
					fFit->hit.nbDesorbed += f->counter[m].hit.nbDesorbed;
					fFit->hit.nbMCHit += f->counter[m].hit.nbMCHit;
					fFit->hit.nbHitEquiv += f->counter[m].hit.nbHitEquiv;
					fFit->hit.sum_1_per_ort_velocity += f->counter[m].hit.sum_1_per_ort_velocity;
					fFit->hit.sum_v_ort += f->counter[m].hit.sum_v_ort;
					fFit->hit.sum_1_per_velocity += f->counter[m].hit.sum_1_per_velocity;
				}

				if (f->sh.isProfile) {
					for (int m = 0; m < (1 + nbMoments); m++) {
						ProfileSlice *shProfile = (ProfileSlice *)(buffer + f->sh.hitOffset + facetHitsSize + m*f->profileSize);
						for (j = 0; j < PROFILE_SIZE; j++) {
							shProfile[j] += f->profile[m][j];
						}
					}
				}

				if (f->sh.isTextured) {
					for (int m = 0; m < (1 + nbMoments); m++) {
						TextureCell *shTexture = (TextureCell *)(buffer + (f->sh.hitOffset + facetHitsSize + f->profileSize*(1 + nbMoments) + m*f->textureSize));
						//double dCoef = gHits->total.hit.nbDesorbed * 1E4 * sHandle->gasMass / 1000 / 6E23 * MAGIC_CORRECTION_FACTOR;  //1E4 is conversion from m2 to cm2
						double timeCorrection = m == 0 ? sHandle->finalOutgassingRate : (sHandle->totalDesorbedMolecules) / sHandle->timeWindowSize;
						//Timecorrection is required to compare constant flow texture values with moment values (for autoscaling)

						for (y = 0; y < f->sh.texHeight; y++) {
							for (x = 0; x < f->sh.texWidth; x++) {
								size_t add = x + y*f->sh.texWidth;

								//Add temporary hit counts
								shTexture[add] += f->texture[m][add];

								double val[3];  //pre-calculated autoscaling values (Pressure, imp.rate, density)

								val[0] = shTexture[add].sum_v_ort_per_area*timeCorrection; //pressure without dCoef_pressure
								val[1] = shTexture[add].countEquiv*f->inc[add] * timeCorrection; //imp.rate without dCoef
								val[2] = f->inc[add] * shTexture[add].sum_1_per_ort_velocity* timeCorrection; //particle density without dCoef

								//Global autoscale
								for (int v = 0; v < 3; v++) {
									if (val[v] > gHits->texture_limits[v].max.all && f->largeEnough[add])
										gHits->texture_limits[v].max.all = val[v];

									if (val[v] > 0.0 && val[v] < gHits->texture_limits[v].min.all && f->largeEnough[add])
										gHits->texture_limits[v].min.all = val[v];

									//Autoscale ignoring constant flow (moments only)
									if (m != 0) {
										if (val[v] > gHits->texture_limits[v].max.moments_only && f->largeEnough[add])
											gHits->texture_limits[v].max.moments_only = val[v];

										if (val[v] > 0.0 && val[v] < gHits->texture_limits[v].min.moments_only && f->largeEnough[add])
											gHits->texture_limits[v].min.moments_only = val[v];
									}
								}
							}
						}
					}
				}

				if (f->sh.countDirection) {
					for (int m = 0; m < (1 + nbMoments); m++) {
						DirectionCell *shDir = (DirectionCell *)(buffer + (f->sh.hitOffset + facetHitsSize + f->profileSize*(1 + nbMoments) + f->textureSize*(1 + nbMoments) + f->directionSize*m));
						for (y = 0; y < f->sh.texHeight; y++) {
							for (x = 0; x < f->sh.texWidth; x++) {
								size_t add = x + y*f->sh.texWidth;
								shDir[add].dir.x += f->direction[m][add].dir.x;
								shDir[add].dir.y += f->direction[m][add].dir.y;
								shDir[add].dir.z += f->direction[m][add].dir.z;
								//shDir[add].sumSpeed += f->direction[m][add].sumSpeed;
								shDir[add].count += f->direction[m][add].count;
							}
						}
					}
				}

				if (f->sh.anglemapParams.record) {
					size_t *shAngleMap = (size_t *)(buffer + f->sh.hitOffset + facetHitsSize + f->profileSize*(1 + nbMoments) + f->textureSize*(1 + nbMoments) + f->directionSize*(1 + nbMoments));
					for (y = 0; y < (f->sh.anglemapParams.thetaLowerRes + f->sh.anglemapParams.thetaHigherRes); y++) {
						for (x = 0; x < f->sh.anglemapParams.phiWidth; x++) {
							size_t add = x + y*f->sh.anglemapParams.phiWidth;
							shAngleMap[add] += f->angleMap.pdf[add];
						}
					}
				}
			} // End if(hitted)
		} // End nbFacet
	} // End nbSuper

	//if there were no textures:
	for (int v = 0; v < 3; v++) {
		if (gHits->texture_limits[v].min.all == HITMAX) gHits->texture_limits[v].min.all = texture_limits_old[v].min.all;
		if (gHits->texture_limits[v].min.moments_only == HITMAX) gHits->texture_limits[v].min.moments_only = texture_limits_old[v].min.moments_only;
		if (gHits->texture_limits[v].max.all == 0.0) gHits->texture_limits[v].max.all = texture_limits_old[v].max.all;
		if (gHits->texture_limits[v].max.moments_only == 0.0) gHits->texture_limits[v].max.moments_only = texture_limits_old[v].max.moments_only;
	}

	ReleaseDataport(dpHit);

	ResetTmpCounters();
	extern char* GetSimuStatus();
	SetState(NULL, GetSimuStatus(), false, true);

#ifdef _DEBUG
	t1 = GetTick();
	printf("Update hits: %f us\n", (t1 - t0)*1000000.0);
#endif

}



void UpdateLog(Dataport * dpLog, DWORD timeout)
{
	if (sHandle->tmpParticleLog.size()) {
		double t0, t1;
		t0 = GetTick();
		SetState(NULL, "Waiting for 'dpLog' dataport access...", false, true);
		sHandle->lastLogUpdateOK = AccessDataportTimed(dpLog, timeout);
		SetState(NULL, "Updating Log...", false, true);
		if (!sHandle->lastLogUpdateOK) return;

		size_t* logBuff = (size_t*)dpLog->buff;
		size_t recordedLogSize = *logBuff;
		ParticleLoggerItem* logBuff2 = (ParticleLoggerItem*)(logBuff + 1);

		size_t writeNb;
		if (recordedLogSize > sHandle->ontheflyParams.logLimit) writeNb = 0;
		else writeNb = Min(sHandle->tmpParticleLog.size(), sHandle->ontheflyParams.logLimit - recordedLogSize);
		memcpy(&logBuff2[recordedLogSize], &sHandle->tmpParticleLog[0], writeNb * sizeof(ParticleLoggerItem)); //Knowing that vector memories are contigious
		(*logBuff) += writeNb;
		ReleaseDataport(dpLog);
		sHandle->tmpParticleLog.clear();
		extern char* GetSimuStatus();
		SetState(NULL, GetSimuStatus(), false, true);

#ifdef _DEBUG
		t1 = GetTick();
		printf("Update hits: %f us\n", (t1 - t0)*1000000.0);
#endif
	}
}

// Compute particle teleport

void PerformTeleport(SubprocessFacet *iFacet) {

	
	//Search destination
	SubprocessFacet *destination;
	bool found = false;
	bool revert = false;
	int destIndex;
	if (iFacet->sh.teleportDest == -1) {
		destIndex = sHandle->teleportedFrom;
		if (destIndex == -1) {
			/*char err[128];
			sprintf(err, "Facet %d tried to teleport to the facet where the particle came from, but there is no such facet.", iFacet->globalId + 1);
			SetErrorSub(err);*/
			RecordHit(HIT_REF);
			sHandle->lastHitFacet = iFacet;
			return; //LEAK
		}
	}
	else destIndex = iFacet->sh.teleportDest - 1;

	//Look in which superstructure is the destination facet:
	for (int i = 0; i < sHandle->nbSuper && (!found); i++) {
		for (int j = 0; j < sHandle->str[i].nbFacet && (!found); j++) {
			if (destIndex == sHandle->str[i].facets[j]->globalId) {
				destination = sHandle->str[i].facets[j];
				sHandle->curStruct = destination->sh.superIdx; //change current superstructure
				sHandle->teleportedFrom = (int)iFacet->globalId; //memorize where the particle came from
				found = true;
			}
		}
	}
	if (!found) {
		/*char err[128];
		sprintf(err, "Teleport destination of facet %d not found (facet %d does not exist)", iFacet->globalId + 1, iFacet->sh.teleportDest);
		SetErrorSub(err);*/
		RecordHit(HIT_REF);
		sHandle->lastHitFacet = iFacet;
		return; //LEAK
	}
	// Count this hit as a transparent pass
	RecordHit(HIT_TELEPORTSOURCE);
	if (iFacet->texture && iFacet->sh.countTrans) RecordHitOnTexture(iFacet, sHandle->flightTimeCurrentParticle, true, 2.0, 2.0);
	if (iFacet->direction && iFacet->sh.countDirection) RecordDirectionVector(iFacet, sHandle->flightTimeCurrentParticle);
	ProfileFacet(iFacet, sHandle->flightTimeCurrentParticle, true, 2.0, 2.0);
	LogHit(iFacet);
	if (iFacet->sh.anglemapParams.record) RecordAngleMap(iFacet);

	// Relaunch particle from new facet
	double inPhi, inTheta;
	std::tie(inTheta, inPhi) = CartesianToPolar(iFacet->sh.nU, iFacet->sh.nV, iFacet->sh.N);
	PolarToCartesian(destination, inTheta, inPhi, false);
	// Move particle to teleport destination point
	double u = iFacet->colU;
	double v = iFacet->colV;
	sHandle->pPos = destination->sh.O + u*destination->sh.U + v*destination->sh.V;
	RecordHit(HIT_TELEPORTDEST);
	int nbTry = 0;
	if (!IsInFacet(*destination, u, v)) { //source and destination facets not the same shape, would generate leak
		// Choose a new starting point
		RecordHit(HIT_ABS);
		bool found = false;
		while (!found && nbTry < 1000) {
			u = rnd();
			v = rnd();
			if (IsInFacet(*destination, u, v)) {
				found = true;
				sHandle->pPos = destination->sh.O + u*destination->sh.U + v*destination->sh.V;
				RecordHit(HIT_DES);
			}
		}
		nbTry++;
	}

	sHandle->lastHitFacet = destination;

	//Count hits on teleport facets
	/*iFacet->sh.counter.hit.nbAbsEquiv++;
	destination->sh.counter.hit.nbDesorbed++;*/

	double ortVelocity = sHandle->velocityCurrentParticle*abs(Dot(sHandle->pDir, iFacet->sh.N));
	//We count a teleport as a local hit, but not as a global one since that would affect the MFP calculation
	/*iFacet->sh.counter.hit.nbMCHit++;
	iFacet->sh.counter.hit.sum_1_per_ort_velocity += 2.0 / ortVelocity;
	iFacet->sh.counter.hit.sum_v_ort += 2.0*(sHandle->useMaxwellDistribution ? 1.0 : 1.1781)*ortVelocity;*/
	IncreaseFacetCounter(iFacet, sHandle->flightTimeCurrentParticle, 1, 0, 0, 2.0 / ortVelocity, 2.0*(sHandle->useMaxwellDistribution ? 1.0 : 1.1781)*ortVelocity);
	iFacet->hitted = true;
	/*destination->sh.counter.hit.sum_1_per_ort_velocity += 2.0 / sHandle->velocityCurrentParticle;
	destination->sh.counter.hit.sum_v_ort += sHandle->velocityCurrentParticle*abs(DOT3(
	sHandle->pDir.x, sHandle->pDir.y, sHandle->pDir.z,
	destination->sh.N.x, destination->sh.N.y, destination->sh.N.z));*/
}

// Perform nbStep simulation steps (a step is a bounce)

bool SimulationMCStep(size_t nbStep) {

	// Perform simulation steps
	for (size_t i = 0; i < nbStep; i++) {
		
		//Prepare output values
		SubprocessFacet   *collidedFacet;
		double   d;
		bool     found;

		std::tie(found,collidedFacet,d) = Intersect(sHandle->pPos, sHandle->pDir, /*sHandle->lastHit,*/THitCache);

		if (found) {

			// Move particle to intersection point
			sHandle->pPos = sHandle->pPos + d*sHandle->pDir;
			//sHandle->distTraveledCurrentParticle += d;

			double lastFLightTime = sHandle->flightTimeCurrentParticle; //memorize for partial hits
			sHandle->flightTimeCurrentParticle += d / 100.0 / sHandle->velocityCurrentParticle; //conversion from cm to m

			if ((!sHandle->calcConstantFlow && (sHandle->flightTimeCurrentParticle > sHandle->latestMoment))
				|| (sHandle->enableDecay && (sHandle->particleDecayMoment < sHandle->flightTimeCurrentParticle))) {
				//hit time over the measured period - we create a new particle
				//OR particle has decayed
				double remainderFlightPath = sHandle->velocityCurrentParticle*100.0*
					Min(sHandle->latestMoment - lastFLightTime, sHandle->particleDecayMoment - lastFLightTime); //distance until the point in space where the particle decayed
				sHandle->distTraveledSinceUpdate_total += remainderFlightPath * sHandle->oriRatio;
				RecordHit(HIT_LAST);
				//sHandle->distTraveledSinceUpdate += sHandle->distTraveledCurrentParticle;
				if (!StartFromSource())
					// desorptionLimit reached
					return false;
			}
			else { //hit within measured time, particle still alive
				if (collidedFacet->sh.teleportDest != 0) {
					sHandle->distTraveledSinceUpdate_total += d * sHandle->oriRatio;
					PerformTeleport(collidedFacet);
				}
				/*else if ((GetOpacityAt(collidedFacet, sHandle->flightTimeCurrentParticle) < 1.0) && (rnd() > GetOpacityAt(collidedFacet, sHandle->flightTimeCurrentParticle))) {
					//Transparent pass
					sHandle->distTraveledSinceUpdate_total += d;
					PerformTransparentPass(collidedFacet);
				}*/
				else {
					sHandle->distTraveledSinceUpdate_total += d * sHandle->oriRatio;
					sHandle->distTraveledSinceUpdate_fullHitsOnly += d * sHandle->oriRatio;
					double stickingProbability = GetStickingAt(collidedFacet, sHandle->flightTimeCurrentParticle);
					if (!sHandle->ontheflyParams.lowFluxMode) { //Regular stick or bounce
						if (stickingProbability == 1.0 || ((stickingProbability > 0.0) && (rnd() < (stickingProbability)))) {
							//Absorbed
							RecordAbsorb(collidedFacet);
							//sHandle->distTraveledSinceUpdate += sHandle->distTraveledCurrentParticle;
							if (!StartFromSource())
								// desorptionLimit reached
								return false;
						}
						else {
							//Reflected
							PerformBounce(collidedFacet);
						}
					}
					else { //Low flux mode
						if (stickingProbability > 0.0) {
							double oriRatioBeforeCollision = sHandle->oriRatio; //Local copy
							sHandle->oriRatio *= (stickingProbability); //Sticking part
							RecordAbsorb(collidedFacet);
							sHandle->oriRatio = oriRatioBeforeCollision * (1.0 - stickingProbability); //Reflected part
						}
						else
							sHandle->oriRatio *= (1.0 - stickingProbability);
						if (sHandle->oriRatio > sHandle->ontheflyParams.lowFluxCutoff) {							
							PerformBounce(collidedFacet);
						}
						else { //eliminate remainder and create new particle
							if (!StartFromSource())
								// desorptionLimit reached
								return false;
						}
					}
				}
			} //end hit within measured time
		} //end intersection found
		else {
			// No intersection found: Leak
			sHandle->nbLeakSinceUpdate++;
			RecordLeakPos();
			if (!StartFromSource())
				// desorptionLimit reached
				return false;
		}
	}
	return true;
}

// Launch a ray from a source facet. The ray 
// direction is chosen according to the desorption type.

bool StartFromSource() {
	bool found = false;
	bool foundInMap = false;
	bool reverse;
	size_t mapPositionW, mapPositionH;
	SubprocessFacet *src = NULL;
	double srcRnd;
	double sumA = 0.0;
	int i = 0, j = 0;
	int nbTry = 0;

	// Check end of simulation
	if (sHandle->ontheflyParams.desorptionLimit > 0) {
		if (sHandle->totalDesorbed >= sHandle->ontheflyParams.desorptionLimit/sHandle->ontheflyParams.nbProcess) {
			sHandle->lastHitFacet = NULL;
			return false;
		}
	}

	// Select source
	srcRnd = rnd() * sHandle->totalDesorbedMolecules;

	while (!found && j < sHandle->nbSuper) { //Go through superstructures
		i = 0;
		while (!found && i < sHandle->str[j].nbFacet) { //Go through facets in a structure
			SubprocessFacet *f = sHandle->str[j].facets[i];
			if (f->sh.desorbType != DES_NONE) { //there is some kind of outgassing
				if (f->sh.useOutgassingFile) { //Using SynRad-generated outgassing map
					if (f->sh.totalOutgassing > 0.0) {
						found = (srcRnd >= sumA) && (srcRnd < (sumA + sHandle->latestMoment * f->sh.totalOutgassing / (1.38E-23*f->sh.temperature)));
						if (found) {
							//look for exact position in map
							double rndRemainder = (srcRnd - sumA) / sHandle->latestMoment*(1.38E-23*f->sh.temperature); //remainder, should be less than f->sh.totalOutgassing
							/*double sumB = 0.0;
							for (w = 0; w < f->sh.outgassingMapWidth && !foundInMap; w++) {
								for (h = 0; h < f->sh.outgassingMapHeight && !foundInMap; h++) {
									double cellOutgassing = f->outgassingMap[h*f->sh.outgassingMapWidth + w];
									if (cellOutgassing > 0.0) {
										foundInMap = (rndRemainder >= sumB) && (rndRemainder < (sumB + cellOutgassing));
										if (foundInMap) mapPositionW = w; mapPositionH = h;
										sumB += cellOutgassing;
									}
								}
							}*/
							double lookupValue = rndRemainder;
							int outgLowerIndex = my_lower_bound(lookupValue, f->outgassingMap, f->sh.outgassingMapWidth*f->sh.outgassingMapHeight); //returns line number AFTER WHICH LINE lookup value resides in ( -1 .. size-2 )
							outgLowerIndex++;
							mapPositionH = (size_t)((double)outgLowerIndex / (double)f->sh.outgassingMapWidth);
							mapPositionW = (size_t)outgLowerIndex - mapPositionH * f->sh.outgassingMapWidth;
							foundInMap = true;
							/*if (!foundInMap) {
								SetErrorSub("Starting point not found in imported desorption map");
								return false;
							}*/
						}
						sumA += sHandle->latestMoment * f->sh.totalOutgassing / (1.38E-23*f->sh.temperature);
					}
				} //end outgassing file block
				else { //constant or time-dependent outgassing
					double facetOutgassing =
						(f->sh.outgassing_paramId >= 0)
						? sHandle->IDs[f->sh.IDid].back().second / (1.38E-23*f->sh.temperature)
						: sHandle->latestMoment*f->sh.outgassing / (1.38E-23*f->sh.temperature);
					found = (srcRnd >= sumA) && (srcRnd < (sumA + facetOutgassing));
					sumA += facetOutgassing;
				} //end constant or time-dependent outgassing block
			} //end 'there is some kind of outgassing'
			if (!found) i++;
			if (f->sh.is2sided) reverse = rnd() > 0.5;
			else reverse = false;
		}
		if (!found) j++;
	}
	if (!found) {
		SetErrorSub("No starting point, aborting");
		return false;
	}
	src = sHandle->str[j].facets[i];

	sHandle->lastHitFacet = src;
	//sHandle->distTraveledCurrentParticle = 0.0;  //for mean free path calculations
	//sHandle->flightTimeCurrentParticle = sHandle->desorptionStartTime + (sHandle->desorptionStopTime - sHandle->desorptionStartTime)*rnd();
	sHandle->flightTimeCurrentParticle = GenerateDesorptionTime(src);
	if (sHandle->useMaxwellDistribution) sHandle->velocityCurrentParticle = GenerateRandomVelocity(src->sh.CDFid);
	else sHandle->velocityCurrentParticle = 145.469*sqrt(src->sh.temperature / sHandle->gasMass);  //sqrt(8*R/PI/1000)=145.47
	sHandle->oriRatio = 1.0;
	if (sHandle->enableDecay) { //decaying gas
		sHandle->particleDecayMoment = sHandle->flightTimeCurrentParticle + sHandle->halfLife*1.44269*-log(rnd()); //1.44269=1/ln2
		//Exponential distribution PDF: probability of 't' life = 1/TAU*exp(-t/TAU) where TAU = half_life/ln2
		//Exponential distribution CDF: probability of life shorter than 't" = 1-exp(-t/TAU)
		//Equation: rnd()=1-exp(-t/TAU)
		//Solution: t=TAU*-log(1-rnd()) and 1-rnd()=rnd() therefore t=half_life/ln2*-log(rnd())
	}
	else {
		sHandle->particleDecayMoment = 1e100; //never decay
	}
	//sHandle->temperature = src->sh.temperature; //Thermalize particle

	found = false; //Starting point within facet

	// Choose a starting point
	while (!found && nbTry < 1000) {
		double u, v;

		if (foundInMap) {
			if (mapPositionW < (src->sh.outgassingMapWidth - 1)) {
				//Somewhere in the middle of the facet
				u = ((double)mapPositionW + rnd()) / src->outgassingMapWidthD;
			}
			else {
				//Last element, prevent from going out of facet
				u = ((double)mapPositionW + rnd() * (src->outgassingMapWidthD - (src->sh.outgassingMapWidth - 1))) / src->outgassingMapWidthD;
			}
			if (mapPositionH < (src->sh.outgassingMapHeight - 1)) {
				//Somewhere in the middle of the facet
				v = ((double)mapPositionH + rnd()) / src->outgassingMapHeightD;
			}
			else {
				//Last element, prevent from going out of facet
				v = ((double)mapPositionH + rnd() * (src->outgassingMapHeightD - (src->sh.outgassingMapHeight - 1))) / src->outgassingMapHeightD;
			}
		}
		else {
			u = rnd();
			v = rnd();
		}
		if (IsInFacet(*src, u, v)) {

			// (U,V) -> (x,y,z)
			sHandle->pPos = src->sh.O + u*src->sh.U + v*src->sh.V;
			src->colU = u;
			src->colV = v;
			found = true;

		}
		nbTry++;
	}

	if (!found) {
		// Get the center, if the center is not included in the facet, a leak is generated.
		if (foundInMap) {
			//double uLength = sqrt(pow(src->sh.U.x, 2) + pow(src->sh.U.y, 2) + pow(src->sh.U.z, 2));
			//double vLength = sqrt(pow(src->sh.V.x, 2) + pow(src->sh.V.y, 2) + pow(src->sh.V.z, 2));
			double u = ((double)mapPositionW + 0.5) / src->outgassingMapWidthD;
			double v = ((double)mapPositionH + 0.5) / src->outgassingMapHeightD;
			sHandle->pPos = src->sh.O + u*src->sh.U + v*src->sh.V;
			src->colU = u;
			src->colV = v;
		}
		else {
			src->colU = 0.5;
			src->colV = 0.5;
			sHandle->pPos = sHandle->str[j].facets[i]->sh.center;
		}

	}

	if (src->sh.isMoving && sHandle->motionType) RecordHit(HIT_MOVING);
	else RecordHit(HIT_DES); //create blue hit point for created particle

	//See docs/theta_gen.png for further details on angular distribution generation
	switch (src->sh.desorbType) {
	case DES_UNIFORM:
		PolarToCartesian(src, acos(rnd()), rnd()*2.0*PI, reverse);
		break;
	case DES_NONE: //for file-based
	case DES_COSINE:
		PolarToCartesian(src, acos(sqrt(rnd())), rnd()*2.0*PI, reverse);
		break;
	case DES_COSINE_N:
		PolarToCartesian(src, acos(pow(rnd(), 1.0 / (src->sh.desorbTypeN + 1.0))), rnd()*2.0*PI, reverse);
		break;
	case DES_ANGLEMAP:
		{
		double theta,phi;
		int thetaLowerIndex;
		double thetaOvershoot;

		std::tie(theta,thetaLowerIndex,thetaOvershoot) = src->angleMap.GenerateThetaFromAngleMap(src->sh.anglemapParams);
		phi = src->angleMap.GeneratePhiFromAngleMap(thetaLowerIndex, thetaOvershoot, src->sh.anglemapParams);
			/*

			size_t angleMapSum = src->angleMapLineSums[(src->sh.anglemapParams.thetaLowerRes + src->sh.anglemapParams.thetaHigherRes) - 1];
			if (angleMapSum == 0) {
				std::stringstream tmp;
				tmp << "Facet " << src->globalId + 1 << ": angle map has all 0 values";
				SetErrorSub(tmp.str().c_str());
				return false;
			}
			double lookupValue = rnd()*(double)angleMapSum; //last element of cumulative distr. = sum
			int thetaLowerIndex = my_lower_bound(lookupValue, src->angleMapLineSums, (src->sh.anglemapParams.thetaLowerRes + src->sh.anglemapParams.thetaHigherRes)); //returns line number AFTER WHICH LINE lookup value resides in ( -1 .. size-2 )
			double thetaOvershoot;
			double theta_a, theta_b, theta_c, theta_cumulative_A, theta_cumulative_B, theta_cumulative_C;
			bool theta_hasThreePoints;
			if (thetaLowerIndex == -1) {//In the first line
				thetaOvershoot = lookupValue / (double)src->angleMapLineSums[0]; 
			}
			else {
				//(lower index can't be last element)
				thetaOvershoot = (lookupValue - (double)src->angleMapLineSums[thetaLowerIndex]) / (double)(src->angleMapLineSums[thetaLowerIndex + 1] - src->angleMapLineSums[thetaLowerIndex]);
			}

			theta_a = GetTheta(src, (double)thetaLowerIndex); theta_cumulative_A = (thetaLowerIndex == -1) ? 0.0 : (double)src->angleMapLineSums[thetaLowerIndex + 0];
			theta_b = GetTheta(src, (double)(thetaLowerIndex+1)); theta_cumulative_B = (double)src->angleMapLineSums[thetaLowerIndex + 1];
			if ((thetaLowerIndex + 2)<(src->sh.anglemapParams.thetaLowerRes + src->sh.anglemapParams.thetaHigherRes)) {
				theta_c = GetTheta(src, (double)(thetaLowerIndex+2)); theta_cumulative_C = (double)src->angleMapLineSums[thetaLowerIndex + 2];
				theta_hasThreePoints = true;
			}
			else {
				theta_hasThreePoints = false;
			}
		
			double theta;
			theta_hasThreePoints = false; //debug
			if (!theta_hasThreePoints) theta = GetTheta(src,thetaLowerIndex + thetaOvershoot);
			else theta = InverseQuadraticInterpolation(lookupValue, theta_a, theta_b, theta_c, theta_cumulative_A, theta_cumulative_B, theta_cumulative_C);

			//Theta determined, let's find phi
			double weigh;
			size_t distr1index, distr2index;
			if (thetaOvershoot < 0.5) {
				distr1index = thetaLowerIndex;
				distr2index = thetaLowerIndex + 1;
				weigh = thetaOvershoot + 0.5;
			}
			else {
				distr1index = thetaLowerIndex + 1;
				distr2index = thetaLowerIndex + 2;
				weigh = thetaOvershoot - 0.5;
			}
			if (distr1index == -1) distr1index++; //In case we interpolate in the first theta range and thetaOvershoot<0.5
			if (distr2index == (src->sh.anglemapParams.thetaLowerRes + src->sh.anglemapParams.thetaHigherRes)) distr2index--; //In case we interpolate in the last theta range and thetaOvershoot>=0.5

			size_t distr1sum = src->angleMapLineSums[distr1index]-((distr1index!=distr2index && distr1index>0)?src->angleMapLineSums[distr1index-1]:0);
			size_t distr2sum = src->angleMapLineSums[distr2index] - ((distr1index != distr2index) ? src->angleMapLineSums[distr2index - 1] : 0);
			double weighedSum=Weigh((double)distr1sum, (double)distr2sum, weigh);
			double phiLookup = rnd() * weighedSum;
			int phiLowerIndex = weighed_lower_bound_X(phiLookup, weigh, &(src->angleMap_pdf[src->sh.anglemapParams.phiWidth*distr1index]), &(src->angleMap_pdf[src->sh.anglemapParams.phiWidth*distr2index]), src->sh.anglemapParams.phiWidth);
		
			////////////////

			double phiOvershoot;
			double phi_a, phi_b, phi_c, phi_cumulative_A, phi_cumulative_B, phi_cumulative_C;
			bool phi_hasThreePoints;
			if (phiLowerIndex == -1) {
				phiOvershoot = phiLookup / Weigh((double)src->angleMap_pdf[src->sh.anglemapParams.phiWidth*distr1index], (double)src->angleMap_pdf[src->sh.anglemapParams.phiWidth*distr2index],weigh);
			}
			else {
				//(lower index can't be last element)
				phiOvershoot = (phiLookup - Weigh((double)src->angleMap_pdf[src->sh.anglemapParams.phiWidth*distr1index + phiLowerIndex], (double)src->angleMap_pdf[src->sh.anglemapParams.phiWidth*distr2index+phiLowerIndex], weigh))
					/ (Weigh((double)src->angleMap_pdf[src->sh.anglemapParams.phiWidth*distr1index + phiLowerIndex + 1], (double)src->angleMap_pdf[src->sh.anglemapParams.phiWidth*distr2index + phiLowerIndex + 1], weigh)
					- Weigh((double)src->angleMap_pdf[src->sh.anglemapParams.phiWidth*distr1index + phiLowerIndex], (double)src->angleMap_pdf[src->sh.anglemapParams.phiWidth*distr2index + phiLowerIndex], weigh));
			}

			phi_a = 2.0*PI*(double)(phiLowerIndex + 1) / (double)src->sh.anglemapParams.phiWidth - PI; phi_cumulative_A = (phiLowerIndex == -1) ? 0.0 : Weigh((double)src->angleMap_pdf[src->sh.anglemapParams.phiWidth*distr1index + phiLowerIndex], (double)src->angleMap_pdf[src->sh.anglemapParams.phiWidth*distr2index + phiLowerIndex], weigh);
			phi_b = 2.0*PI*(double)(phiLowerIndex + 2) / (double)src->sh.anglemapParams.phiWidth - PI; phi_cumulative_B = Weigh((double)src->angleMap_pdf[src->sh.anglemapParams.phiWidth*distr1index + phiLowerIndex + 1], (double)src->angleMap_pdf[src->sh.anglemapParams.phiWidth*distr2index + phiLowerIndex + 1], weigh);
			if ((phiLowerIndex + 2)<(src->sh.anglemapParams.phiWidth)) {
				phi_c = 2.0*PI*(double)(phiLowerIndex + 3) / (double)src->sh.anglemapParams.phiWidth - PI; phi_cumulative_C = Weigh((double)src->angleMap_pdf[src->sh.anglemapParams.phiWidth*distr1index + phiLowerIndex + 2], (double)src->angleMap_pdf[src->sh.anglemapParams.phiWidth*distr2index + phiLowerIndex + 2], weigh);
				phi_hasThreePoints = true;
			}
			else {
				phi_hasThreePoints = false;
			}

			double phi;
			phi_hasThreePoints = false; //debug
			if (!phi_hasThreePoints) phi = 2.0*PI*((double)(phiLowerIndex + 1) + phiOvershoot) / (double)src->sh.anglemapParams.phiWidth - PI;
			else phi = InverseQuadraticInterpolation(phiLookup, phi_a, phi_b, phi_c, phi_cumulative_A, phi_cumulative_B, phi_cumulative_C);

			/////////////////////////////
			*/
			PolarToCartesian(src, PI - theta, phi, false); //angle map contains incident angle (between N and source dir) and theta is dir (between N and dest dir)
			_ASSERTE(sHandle->pDir.x == sHandle->pDir.x);
			
		}
	}

	// Current structure
	sHandle->curStruct = src->sh.superIdx;
	sHandle->teleportedFrom = -1;

	// Count

	src->hitted = true;
	sHandle->totalDesorbed++;
	sHandle->tmpGlobalCount.hit.nbDesorbed++;
	//sHandle->nbPHit = 0;

	if (src->sh.isMoving) {
		TreatMovingFacet();
	}

	double ortVelocity = sHandle->velocityCurrentParticle*abs(Dot(sHandle->pDir, src->sh.N));
	/*src->sh.counter.hit.nbDesorbed++;
	src->sh.counter.hit.sum_1_per_ort_velocity += 2.0 / ortVelocity; //was 2.0 / ortV
	src->sh.counter.hit.sum_v_ort += (sHandle->useMaxwellDistribution ? 1.0 : 1.1781)*ortVelocity;*/
	IncreaseFacetCounter(src, sHandle->flightTimeCurrentParticle, 0, 1, 0, 2.0 / ortVelocity, (sHandle->useMaxwellDistribution ? 1.0 : 1.1781)*ortVelocity);
	//Desorption doesn't contribute to angular profiles, nor to angle maps
	ProfileFacet(src, sHandle->flightTimeCurrentParticle, false, 2.0, 1.0); //was 2.0, 1.0
	LogHit(src);
	if (src->texture && src->sh.countDes) RecordHitOnTexture(src, sHandle->flightTimeCurrentParticle, true, 2.0, 1.0); //was 2.0, 1.0
	//if (src->direction && src->sh.countDirection) RecordDirectionVector(src, sHandle->flightTimeCurrentParticle);

	// Reset volatile state
	if (sHandle->hasVolatile) {
		for (j = 0; j < sHandle->nbSuper; j++) {
			for (i = 0; i < sHandle->str[j].nbFacet; i++) {
				sHandle->str[j].facets[i]->ready = true;
			}
		}
	}

	found = false;
	return true;
}

std::tuple<double, int, double> Anglemap::GenerateThetaFromAngleMap(const AnglemapParams& anglemapParams)
{
	double lookupValue = rnd();
	int thetaLowerIndex = my_lower_bound(lookupValue, theta_CDF, (anglemapParams.thetaLowerRes + anglemapParams.thetaHigherRes)); //returns line number AFTER WHICH LINE lookup value resides in ( -1 .. size-2 )
	double theta, thetaOvershoot;

	if (thetaLowerIndex == -1) { //first half section
		thetaOvershoot = 0.5 + 0.5 * lookupValue / theta_CDF[0]; //between 0.5 and 1
		theta = GetTheta((double)thetaLowerIndex + 0.5 + thetaOvershoot, anglemapParams); //between 0 and the first section end
		return std::tie(theta, thetaLowerIndex, thetaOvershoot);
	}
	else if (thetaLowerIndex == (anglemapParams.thetaLowerRes + anglemapParams.thetaHigherRes - 1)) { //last half section //can this happen?
		thetaOvershoot = 0.5 * (lookupValue - theta_CDF[thetaLowerIndex])
			/ (1.0 - theta_CDF[thetaLowerIndex]); //between 0 and 0.5
		theta = GetTheta((double)thetaLowerIndex + 0.5 + thetaOvershoot, anglemapParams); //between 0 and the first section end
		return std::tie(theta, thetaLowerIndex, thetaOvershoot);
	}
	else { //regular section
		if (/*true || */phi_CDFsums[thetaLowerIndex] == phi_CDFsums[thetaLowerIndex + 1]) {
			//The pdf's slope is 0, linear interpolation
			thetaOvershoot = (lookupValue - theta_CDF[thetaLowerIndex]) / (theta_CDF[thetaLowerIndex + 1] - theta_CDF[thetaLowerIndex]);
			theta = GetTheta((double)thetaLowerIndex + 0.5 + thetaOvershoot,anglemapParams);
		}
		else {
			//2nd degree interpolation
			// y(x) = ax^2 + bx + c
			// c: CDF value at lower index
			// b: pdf value at lower index
			// a: pdf slope at lower index / 2
			// dy := y - c
			// dx := x - [x at lower index]
			// dy = ax^2 + bx
			// dx = ( -b + sqrt(b^2 +4*a*dy) ) / (2a)
			double thetaStep = GetTheta((double)thetaLowerIndex + 1.5, anglemapParams) - GetTheta((double)thetaLowerIndex + 0.5, anglemapParams);
			double c = theta_CDF[thetaLowerIndex]; //CDF value at lower index
			double b = (double)phi_CDFsums[thetaLowerIndex] / (double)theta_CDFsum / thetaStep; //pdf value at lower index
			double a = 0.5 * ((double)(phi_CDFsums[thetaLowerIndex + 1]) - (double)phi_CDFsums[thetaLowerIndex])/(double)theta_CDFsum/Sqr(thetaStep); //pdf slope at lower index
			double dy = lookupValue - c;

			double dx = (-b + sqrt(Sqr(b) + 4 * a*dy)) / (2 * a); //Since b>=0 it's the + branch of the +- that is valid for us

			thetaOvershoot = dx/thetaStep;
			theta = GetTheta((double)thetaLowerIndex + 0.5 + thetaOvershoot,anglemapParams);
		}
	}
	_ASSERTE(theta == theta);
	return std::tie(theta, thetaLowerIndex, thetaOvershoot);
}

double Anglemap::GeneratePhiFromAngleMap(const int & thetaLowerIndex, const double & thetaOvershoot, const AnglemapParams & anglemapParams)
{
	double lookupValue = rnd();
	if (anglemapParams.phiWidth == 1) return -PI + 2.0 * PI * lookupValue; //special case, uniform phi distribution
	int phiLowerIndex;
	double weigh; //0: take previous theta line, 1: take next theta line, 0..1: interpolate in-between
	if (thetaLowerIndex == -1) { //first theta half section
		lookupValue += phi_CDFs[0]; //periodic BCs over -PI...PI, can be larger than 1
		phiLowerIndex = my_lower_bound(lookupValue, &phi_CDFs[0], anglemapParams.phiWidth); //take entirely the phi ditro belonging to first theta
		weigh = thetaOvershoot; // [0.5 - 1], will subtract 0.5 when evaluating thetaIndex
	}
	else if (thetaLowerIndex == (anglemapParams.thetaLowerRes + anglemapParams.thetaHigherRes - 1)) { //last theta half section
		lookupValue += phi_CDFs[thetaLowerIndex*anglemapParams.phiWidth]; //periodic BCs over -PI...PI, can be larger than 1
		phiLowerIndex = my_lower_bound(lookupValue, &phi_CDFs[thetaLowerIndex*anglemapParams.phiWidth], anglemapParams.phiWidth); //take entirely the phi ditro belonging to latest theta
		weigh = thetaOvershoot; // [0 - 0.5], will add 0.5 when evaluating thetaIndex
	}
	else {
		//Here we do a weighing both by the hit sum of the previous and next lines (w1 and w2) and also the weighs of the two lines based on thetaOvershoot (w3 and w4)
		// w1: sum of hits in previous line
		// w2: sum of hits in next line
		// w3: weigh of previous line (1 - thetaOvershoot)
		// w4: weigh of next line     (thetaOvershoot)
		// result: previous value weight: w1*w3 / (w1*w3 + w2*w4)
		//         next     value weight: w2*w4 / (w1*w3 + w2*w4) <- this will be the input for weighed_lower_bound
		
		double div,weigh;
		div = ((double)phi_CDFsums[thetaLowerIndex] * (1.0 - thetaOvershoot) + (double)phi_CDFsums[thetaLowerIndex + 1] * thetaOvershoot); // (w1*w3 + w2*w4)
		if (div>0.0) {
			weigh = (thetaOvershoot * (double)phi_CDFsums[thetaLowerIndex + 1]) / div;    //      w2*w4 / (w1*w3 + w2*w4)
		}
		else {
			weigh = thetaOvershoot;
		}
		lookupValue += Weigh((double)phi_CDFs[thetaLowerIndex*anglemapParams.phiWidth], (double)phi_CDFs[(thetaLowerIndex + 1) * anglemapParams.phiWidth], weigh);
		phiLowerIndex = weighed_lower_bound_X(lookupValue, weigh, &phi_CDFs[thetaLowerIndex*anglemapParams.phiWidth], &phi_CDFs[(thetaLowerIndex + 1) * anglemapParams.phiWidth], anglemapParams.phiWidth);
	}

	double phi, phiOvershoot;
	double thetaIndex = (double)thetaLowerIndex + 0.5 + weigh;
	if (phiLowerIndex == -1) { //first half section
		__debugbreak(); //should not happen since we shifted the lookup value with first value
		phiOvershoot = 0.5 + 0.5 * lookupValue / GetPhiCDFValue(thetaIndex, 0, anglemapParams); //between 0.5 and 1
		phi = GetPhi((double)phiLowerIndex + 0.5 + phiOvershoot, anglemapParams); //between 0 and the first section end
	}
	/*else if (phiLowerIndex == (anglemapParams.phiWidth - 1)) { //last half section
		phiOvershoot = 0.5 * (lookupValue - GetPhiCDFValue(thetaIndex, phiLowerIndex, anglemapParams) )
			/ (1.0 - GetPhiCDFValue(thetaIndex, phiLowerIndex, anglemapParams)); //between 0 and 0.5
		phi = GetPhi((double)phiLowerIndex + 0.5 + phiOvershoot, anglemapParams); //between 0 and the first section end
	}*/
	else { //regular or last section 
		if (/*true ||*/ GetPhipdfValue(thetaIndex, phiLowerIndex, anglemapParams) == GetPhipdfValue(thetaIndex, phiLowerIndex+1, anglemapParams)) {
			//The pdf's slope is 0, linear interpolation
			phiOvershoot = (lookupValue - GetPhiCDFValue(thetaIndex, phiLowerIndex, anglemapParams))
				/ (GetPhiCDFValue(thetaIndex, phiLowerIndex + 1, anglemapParams) - GetPhiCDFValue(thetaIndex, phiLowerIndex, anglemapParams));
			phi = GetPhi((double)phiLowerIndex + 0.5 + phiOvershoot, anglemapParams);
		}
		else {

			//2nd degree interpolation
			// y(x) = ax^2 + bx + c
			// c: CDF value at lower index
			// b: pdf value at lower index
			// a: pdf slope at lower index / 2
			// dy := y - c
			// dx := x - [x at lower index]
			// dy = ax^2 + bx
			// dx = ( -b + sqrt(b^2 +4*a*dy) ) / (2a)
			double phiStep = 2.0 * PI / (double)anglemapParams.phiWidth;
			double c = GetPhiCDFValue(thetaIndex, phiLowerIndex, anglemapParams); //CDF value at lower index
			double b = GetPhipdfValue(thetaIndex,phiLowerIndex,anglemapParams) / GetPhiCDFSum(thetaIndex,anglemapParams) / phiStep; //pdf value at lower index
			double a = 0.5 * (GetPhipdfValue(thetaIndex, phiLowerIndex + 1, anglemapParams) - GetPhipdfValue(thetaIndex, phiLowerIndex, anglemapParams)) / GetPhiCDFSum(thetaIndex, anglemapParams) / Sqr(phiStep); //pdf slope at lower index
			double dy = lookupValue - c;

			double D = Sqr(b) + 4 * a*dy; //Discriminant. In rare cases it might be slightly negative, then fall back to linear interpolation:
			if (D < 0) {
				phiOvershoot = (lookupValue - GetPhiCDFValue(thetaIndex, phiLowerIndex, anglemapParams))
					/ (GetPhiCDFValue(thetaIndex, (int)IDX(phiLowerIndex + 1, anglemapParams.phiWidth), anglemapParams) - GetPhiCDFValue(thetaIndex, phiLowerIndex, anglemapParams));
			}
			else {
				double dx = (-b + sqrt(Sqr(b) + 4 * a*dy)) / (2 * a); //Since b>=0 it's the + branch of the +- that is valid for us
				phiOvershoot = dx / phiStep;
			}
			phi = GetPhi((double)phiLowerIndex + 0.5 + phiOvershoot, anglemapParams);
		}
	}
	_ASSERTE(phi == phi);
	_ASSERTE(phi > -PI && phi < PI);
	return phi;
}

double Anglemap::GetTheta(const double& thetaIndex, const AnglemapParams& anglemapParams)
{
	if ((size_t)(thetaIndex) < anglemapParams.thetaLowerRes) { // 0 < theta < limit
		return anglemapParams.thetaLimit * (thetaIndex) / (double)anglemapParams.thetaLowerRes;
	}
	else { // limit < theta < PI/2
		return anglemapParams.thetaLimit + (PI/2.0- anglemapParams.thetaLimit) * (thetaIndex - (double)anglemapParams.thetaLowerRes) / (double)anglemapParams.thetaHigherRes;
	}
}

double Anglemap::GetPhi(const double & phiIndex, const AnglemapParams & anglemapParams)
//makes phiIndex circular and converts from index to -pi...pi
{
	double width = (double)anglemapParams.phiWidth;
	double correctedIndex = (phiIndex < width) ? phiIndex : phiIndex - width;
	return - PI + 2.0 * PI * correctedIndex / width;
}

double Anglemap::GetPhipdfValue(const double & thetaIndex, const int & phiLowerIndex, const AnglemapParams& anglemapParams)
//phiLowerIndex is circularized
{
	if (thetaIndex < 0.5) {
		return (double)pdf[IDX(phiLowerIndex, anglemapParams.phiWidth)];
	}
	else if (thetaIndex >(double)(anglemapParams.thetaLowerRes + anglemapParams.thetaHigherRes) - 0.5) {
		return (double)pdf[anglemapParams.phiWidth * (anglemapParams.thetaLowerRes + anglemapParams.thetaHigherRes - 1) + IDX(phiLowerIndex, anglemapParams.phiWidth)];
	}
	else {
		size_t thetaLowerIndex = (size_t)(thetaIndex - 0.5);
		double thetaOvershoot = thetaIndex - 0.5 - (double)thetaLowerIndex;
		double valueFromLowerpdf = (double)pdf[anglemapParams.phiWidth * thetaLowerIndex + IDX(phiLowerIndex, anglemapParams.phiWidth)];
		double valueFromHigherpdf = (double)pdf[anglemapParams.phiWidth * (thetaLowerIndex + 1) + IDX(phiLowerIndex, anglemapParams.phiWidth)];
		return Weigh(valueFromLowerpdf, valueFromHigherpdf, thetaOvershoot);
	}
}

double Anglemap::GetPhiCDFValue(const double & thetaIndex, const int & phiLowerIndex, const AnglemapParams& anglemapParams)
{
	if (thetaIndex < 0.5) {
		return (phiLowerIndex < anglemapParams.phiWidth) ? phi_CDFs[phiLowerIndex] : 1.0 + phi_CDFs[0];
	}
	else if (thetaIndex >(double)(anglemapParams.thetaLowerRes + anglemapParams.thetaHigherRes) - 0.5) {
		return (phiLowerIndex < anglemapParams.phiWidth) ? phi_CDFs[anglemapParams.phiWidth * (anglemapParams.thetaLowerRes + anglemapParams.thetaHigherRes - 1) + phiLowerIndex] : 1.0 + phi_CDFs[anglemapParams.phiWidth * (anglemapParams.thetaLowerRes + anglemapParams.thetaHigherRes - 1)];
	}
	else {
		size_t thetaLowerIndex = (size_t)(thetaIndex - 0.5);
		double thetaOvershoot = thetaIndex - 0.5 - (double)thetaLowerIndex;
		double valueFromLowerCDF = (phiLowerIndex < anglemapParams.phiWidth) ? phi_CDFs[anglemapParams.phiWidth * thetaLowerIndex + phiLowerIndex] : 1.0 + phi_CDFs[anglemapParams.phiWidth * (thetaLowerIndex)];
		double valueFromHigherCDF = (phiLowerIndex < anglemapParams.phiWidth) ? phi_CDFs[anglemapParams.phiWidth * (thetaLowerIndex + 1) + phiLowerIndex] : 1.0 + phi_CDFs[anglemapParams.phiWidth * (thetaLowerIndex + 1)];
		return Weigh(valueFromLowerCDF, valueFromHigherCDF, thetaOvershoot);
	}

}

double Anglemap::GetPhiCDFSum(const double & thetaIndex, const AnglemapParams& anglemapParams)
{
	if (thetaIndex < 0.5) {
		return (double)phi_CDFsums[0];
	}
	else if (thetaIndex >(double)(anglemapParams.thetaLowerRes + anglemapParams.thetaHigherRes) - 0.5) {
		return (double)phi_CDFsums[anglemapParams.thetaLowerRes + anglemapParams.thetaHigherRes - 1];
	}
	else {
		size_t thetaLowerIndex = (size_t)(thetaIndex - 0.5);
		double thetaOvershoot = thetaIndex - 0.5 - (double)thetaLowerIndex;
		double valueFromLowerSum = (double)phi_CDFsums[thetaLowerIndex];
		double valueFromHigherSum = (double)phi_CDFsums[thetaLowerIndex + 1];
		return Weigh(valueFromLowerSum, valueFromHigherSum, thetaOvershoot);
	}
}

void PerformBounce(SubprocessFacet *iFacet) {

	bool revert = false;

	sHandle->tmpGlobalCount.hit.nbMCHit++; //global
	sHandle->tmpGlobalCount.hit.nbHitEquiv += sHandle->oriRatio;

	// Handle super structure link facet. Can be 
	if (iFacet->sh.superDest) {
		IncreaseFacetCounter(iFacet, sHandle->flightTimeCurrentParticle, 1, 0, 0, 0, 0);
		sHandle->curStruct = iFacet->sh.superDest - 1;
		if (iFacet->sh.isMoving) { //A very special case where link facets can be used as transparent but moving facets
			RecordHit(HIT_MOVING);
			TreatMovingFacet();
		}
		else {
			// Count this hit as a transparent pass
			RecordHit(HIT_TRANS);
		}
		LogHit(iFacet);
		ProfileFacet(iFacet, sHandle->flightTimeCurrentParticle, true, 2.0, 2.0);
		if (iFacet->sh.anglemapParams.record) RecordAngleMap(iFacet);
		if (iFacet->texture && iFacet->sh.countTrans) RecordHitOnTexture(iFacet, sHandle->flightTimeCurrentParticle, true, 2.0, 2.0);
		if (iFacet->direction && iFacet->sh.countDirection) RecordDirectionVector(iFacet, sHandle->flightTimeCurrentParticle);
		return;

	}

	// Handle volatile facet
	if (iFacet->sh.isVolatile) {

		if (iFacet->ready) {
			IncreaseFacetCounter(iFacet, sHandle->flightTimeCurrentParticle, 0, 0, 1, 0, 0);
			iFacet->ready = false;
			LogHit(iFacet);
			ProfileFacet(iFacet, sHandle->flightTimeCurrentParticle, true, 2.0, 1.0);
			if (iFacet->texture && iFacet->sh.countAbs) RecordHitOnTexture(iFacet, sHandle->flightTimeCurrentParticle, true, 2.0, 1.0);
			if (iFacet->direction && iFacet->sh.countDirection) RecordDirectionVector(iFacet, sHandle->flightTimeCurrentParticle);
		}
		return;

	}

	if (iFacet->sh.is2sided) {
		// We may need to revert normal in case of 2 sided hit
		revert = Dot(sHandle->pDir, iFacet->sh.N) > 0.0;
	}

	//Texture/Profile incoming hit
	

	//Register (orthogonal) velocity
	double ortVelocity = sHandle->velocityCurrentParticle*abs(Dot(sHandle->pDir, iFacet->sh.N));

	/*iFacet->sh.counter.hit.nbMCHit++; //hit facet
	iFacet->sh.counter.hit.sum_1_per_ort_velocity += 1.0 / ortVelocity;
	iFacet->sh.counter.hit.sum_v_ort += (sHandle->useMaxwellDistribution ? 1.0 : 1.1781)*ortVelocity;*/

	IncreaseFacetCounter(iFacet, sHandle->flightTimeCurrentParticle, 1, 0, 0, 1.0 / ortVelocity, (sHandle->useMaxwellDistribution ? 1.0 : 1.1781)*ortVelocity);
	if (iFacet->texture && iFacet->sh.countRefl) RecordHitOnTexture(iFacet, sHandle->flightTimeCurrentParticle, true, 1.0, 1.0);
	if (iFacet->direction && iFacet->sh.countDirection) RecordDirectionVector(iFacet, sHandle->flightTimeCurrentParticle);
	LogHit(iFacet);
	ProfileFacet(iFacet, sHandle->flightTimeCurrentParticle, true, 1.0, 1.0);
	if (iFacet->sh.anglemapParams.record) RecordAngleMap(iFacet);

	// Relaunch particle
	UpdateVelocity(iFacet);
	//Sojourn time
	if (iFacet->sh.enableSojournTime) {
		double A = exp(-iFacet->sh.sojournE / (8.31*iFacet->sh.temperature));
		sHandle->flightTimeCurrentParticle += -log(rnd()) / (A*iFacet->sh.sojournFreq);
	}

	if (iFacet->sh.reflection.diffusePart > 0.999999) { //Speedup branch for most common, diffuse case
		PolarToCartesian(iFacet, acos(sqrt(rnd())), rnd()*2.0*PI, revert);
	}
	else {
		double reflTypeRnd = rnd();
		if (reflTypeRnd < iFacet->sh.reflection.diffusePart)
		{
			//diffuse reflection
			//See docs/theta_gen.png for further details on angular distribution generation
			PolarToCartesian(iFacet, acos(sqrt(rnd())), rnd()*2.0*PI, revert);
		}
		else  if (reflTypeRnd < (iFacet->sh.reflection.diffusePart + iFacet->sh.reflection.specularPart))
		{
			//specular reflection
			double inTheta, inPhi;
			std::tie(inTheta, inPhi) = CartesianToPolar(iFacet->sh.nU, iFacet->sh.nV, iFacet->sh.N);
			PolarToCartesian(iFacet, PI - inTheta, inPhi, false);

		}
		else {
			//Cos^N reflection
			PolarToCartesian(iFacet, acos(pow(rnd(), 1.0 / (iFacet->sh.reflection.cosineExponent + 1.0))), rnd()*2.0*PI, revert);
		}
	}

	if (iFacet->sh.isMoving) {
		TreatMovingFacet();
	}

	//Texture/Profile outgoing particle
	//Register outgoing velocity
	ortVelocity = sHandle->velocityCurrentParticle*abs(Dot(sHandle->pDir, iFacet->sh.N));

	/*iFacet->sh.counter.hit.sum_1_per_ort_velocity += 1.0 / ortVelocity;
	iFacet->sh.counter.hit.sum_v_ort += (sHandle->useMaxwellDistribution ? 1.0 : 1.1781)*ortVelocity;*/
	IncreaseFacetCounter(iFacet, sHandle->flightTimeCurrentParticle, 0, 0, 0, 1.0 / ortVelocity, (sHandle->useMaxwellDistribution ? 1.0 : 1.1781)*ortVelocity);
	if (iFacet->texture && iFacet->sh.countRefl) RecordHitOnTexture(iFacet, sHandle->flightTimeCurrentParticle, false, 1.0, 1.0); //count again for outward velocity
	ProfileFacet(iFacet, sHandle->flightTimeCurrentParticle, false, 1.0, 1.0);
	//no direction count on outgoing, neither angle map

	if (iFacet->sh.isMoving && sHandle->motionType) RecordHit(HIT_MOVING);
	else RecordHit(HIT_REF);
	sHandle->lastHitFacet = iFacet;
	//sHandle->nbPHit++;
}

void PerformTransparentPass(SubprocessFacet *iFacet) { //disabled, caused finding hits with the same facet
	/*double directionFactor = abs(DOT3(
		sHandle->pDir.x, sHandle->pDir.y, sHandle->pDir.z,
		iFacet->sh.N.x, iFacet->sh.N.y, iFacet->sh.N.z));
	iFacet->sh.counter.hit.nbMCHit++;
	iFacet->sh.counter.hit.sum_1_per_ort_velocity += 2.0 / (sHandle->velocityCurrentParticle*directionFactor);
	iFacet->sh.counter.hit.sum_v_ort += 2.0*(sHandle->useMaxwellDistribution ? 1.0 : 1.1781)*sHandle->velocityCurrentParticle*directionFactor;
	iFacet->hitted = true;
	if (iFacet->texture && iFacet->sh.countTrans) RecordHitOnTexture(iFacet, sHandle->flightTimeCurrentParticle + iFacet->colDist / 100.0 / sHandle->velocityCurrentParticle,
		true, 2.0, 2.0);
	if (iFacet->direction && iFacet->sh.countDirection) RecordDirectionVector(iFacet, sHandle->flightTimeCurrentParticle + iFacet->colDist / 100.0 / sHandle->velocityCurrentParticle);
	ProfileFacet(iFacet, sHandle->flightTimeCurrentParticle + iFacet->colDist / 100.0 / sHandle->velocityCurrentParticle,
		true, 2.0, 2.0);
	RecordHit(HIT_TRANS);
	sHandle->lastHit = iFacet;*/
}

void RecordAbsorb(SubprocessFacet *iFacet) {
	sHandle->tmpGlobalCount.hit.nbMCHit++; //global	
	sHandle->tmpGlobalCount.hit.nbHitEquiv += sHandle->oriRatio;
	sHandle->tmpGlobalCount.hit.nbAbsEquiv += sHandle->oriRatio;
	RecordHit(HIT_ABS);
	double ortVelocity = sHandle->velocityCurrentParticle*abs(Dot(sHandle->pDir, iFacet->sh.N));
	IncreaseFacetCounter(iFacet, sHandle->flightTimeCurrentParticle, 1, 0, 1, 2.0 / ortVelocity, (sHandle->useMaxwellDistribution ? 1.0 : 1.1781)*ortVelocity);
	LogHit(iFacet);
	ProfileFacet(iFacet, sHandle->flightTimeCurrentParticle, true, 2.0, 1.0); //was 2.0, 1.0
	if (iFacet->sh.anglemapParams.record) RecordAngleMap(iFacet);
	if (iFacet->texture && iFacet->sh.countAbs) RecordHitOnTexture(iFacet, sHandle->flightTimeCurrentParticle, true, 2.0, 1.0); //was 2.0, 1.0
	if (iFacet->direction && iFacet->sh.countDirection) RecordDirectionVector(iFacet, sHandle->flightTimeCurrentParticle);
}

void RecordHitOnTexture(SubprocessFacet *f, double time, bool countHit, double velocity_factor, double ortSpeedFactor) {

	size_t tu = (size_t)(f->colU * f->sh.texWidthD);
	size_t tv = (size_t)(f->colV * f->sh.texHeightD);
	size_t add = tu + tv*(f->sh.texWidth);
	double ortVelocity = (sHandle->useMaxwellDistribution ? 1.0 : 1.1781)*sHandle->velocityCurrentParticle*abs(Dot(sHandle->pDir, f->sh.N)); //surface-orthogonal velocity component

	for (size_t m = 0; m <= sHandle->nbMoments; m++)
		if (m == 0 || abs(time - sHandle->moments[m - 1]) < sHandle->timeWindowSize / 2.0) {
			if (countHit) f->texture[m][add].countEquiv += sHandle->oriRatio;
			f->texture[m][add].sum_1_per_ort_velocity += sHandle->oriRatio * velocity_factor / ortVelocity;
			f->texture[m][add].sum_v_ort_per_area += sHandle->oriRatio * ortSpeedFactor*ortVelocity*f->inc[add]; // sum ortho_velocity[m/s] / cell_area[cm2]
		}
}

void RecordDirectionVector(SubprocessFacet *f, double time) {
	size_t tu = (size_t)(f->colU * f->sh.texWidthD);
	size_t tv = (size_t)(f->colV * f->sh.texHeightD);
	size_t add = tu + tv*(f->sh.texWidth);

	for (size_t m = 0; m <= sHandle->moments.size(); m++) {
		if (m == 0 || abs(time - sHandle->moments[m - 1]) < sHandle->timeWindowSize / 2.0) {
			f->direction[m][add].dir = f->direction[m][add].dir + sHandle->oriRatio * sHandle->pDir * sHandle->velocityCurrentParticle;
			f->direction[m][add].count++;
		}
	}

}

void ProfileFacet(SubprocessFacet *f, double time, bool countHit, double velocity_factor, double ortSpeedFactor) {

	size_t nbMoments = sHandle->moments.size();

	if (countHit && f->sh.profileType == PROFILE_ANGULAR) {
			double dot = Dot(f->sh.N, sHandle->pDir);
			double theta = acos(abs(dot));     // Angle to normal (PI/2 => PI)
			size_t pos = (size_t)(theta / (PI / 2)*((double)PROFILE_SIZE)); // To Grad
			Saturate(pos, 0, PROFILE_SIZE - 1);
			for (size_t m = 0; m <= nbMoments; m++) {
				if (m == 0 || abs(time - sHandle->moments[m - 1]) < sHandle->timeWindowSize / 2.0) {
					f->profile[m][pos].countEquiv += sHandle->oriRatio;
				}
			}
	}
	else if (Contains({PROFILE_U,PROFILE_V},f->sh.profileType)) {
		size_t pos = (size_t)((f->sh.profileType == PROFILE_U ? f->colU : f->colV)*(double)PROFILE_SIZE);
		if (pos >= 0 && pos < PROFILE_SIZE) {
			for (size_t m = 0; m <= nbMoments; m++) {
				if (m == 0 || abs(time - sHandle->moments[m - 1]) < sHandle->timeWindowSize / 2.0) {
					if (countHit) f->profile[m][pos].countEquiv += sHandle->oriRatio;
					double ortVelocity = sHandle->velocityCurrentParticle*abs(Dot(f->sh.N, sHandle->pDir));
					f->profile[m][pos].sum_1_per_ort_velocity += sHandle->oriRatio * velocity_factor / ortVelocity;
					f->profile[m][pos].sum_v_ort += sHandle->oriRatio * ortSpeedFactor*(sHandle->useMaxwellDistribution ? 1.0 : 1.1781)*ortVelocity;
				}
			}
		}
	}
	else if (countHit && Contains({ PROFILE_VELOCITY,PROFILE_ORT_VELOCITY,PROFILE_TAN_VELOCITY },f->sh.profileType)) {
		double dot;
		if (f->sh.profileType == PROFILE_VELOCITY) {
			dot = 1.0;
		}
		else if (f->sh.profileType == PROFILE_ORT_VELOCITY) {
			dot = abs(Dot(f->sh.N, sHandle->pDir));  //cos(theta) as "dot" value
		}
		else { //Tangential
			dot = sqrt(1 - Sqr(abs(Dot(f->sh.N, sHandle->pDir))));  //tangential
		}
		size_t pos = (size_t)(dot*sHandle->velocityCurrentParticle / f->sh.maxSpeed*(double)PROFILE_SIZE); //"dot" default value is 1.0
		if (pos >= 0 && pos < PROFILE_SIZE) {
			for (size_t m = 0; m <= nbMoments; m++) {
				if (m == 0 || abs(time - sHandle->moments[m - 1]) < sHandle->timeWindowSize / 2.0) {
					f->profile[m][pos].countEquiv += sHandle->oriRatio;
				}
			}
		}
	}	
}

void LogHit(SubprocessFacet * f)
{
	if (sHandle->ontheflyParams.enableLogging &&
		sHandle->ontheflyParams.logFacetId == f->globalId &&
		sHandle->tmpParticleLog.size() < (sHandle->ontheflyParams.logLimit / sHandle->ontheflyParams.nbProcess)) {
		ParticleLoggerItem log;
		log.facetHitPosition = Vector2d(f->colU, f->colV);
		std::tie(log.hitTheta,log.hitPhi) = CartesianToPolar(f->sh.nU, f->sh.nV, f->sh.N);
		log.oriRatio = sHandle->oriRatio;
		log.particleDecayMoment = sHandle->particleDecayMoment;
		log.time = sHandle->flightTimeCurrentParticle;
		log.velocity = sHandle->velocityCurrentParticle;
		sHandle->tmpParticleLog.push_back(log);
	}
}

void RecordAngleMap(SubprocessFacet* collidedFacet) {
	double inTheta, inPhi;
	std::tie(inTheta,inPhi) = CartesianToPolar(collidedFacet->sh.nU, collidedFacet->sh.nV, collidedFacet->sh.N);
	if (inTheta > PI / 2.0) inTheta = abs(PI - inTheta); //theta is originally respective to N, but we'd like the angle between 0 and PI/2
	bool countTheta = true;
	size_t thetaIndex;
	if (inTheta < collidedFacet->sh.anglemapParams.thetaLimit) {
		if (collidedFacet->sh.anglemapParams.thetaLowerRes > 0) {
			thetaIndex = (size_t)(inTheta / collidedFacet->sh.anglemapParams.thetaLimit*(double)collidedFacet->sh.anglemapParams.thetaLowerRes);
		}
		else {
			countTheta = false;
		}
	}
	else {
		if (collidedFacet->sh.anglemapParams.thetaHigherRes > 0) {
			thetaIndex = collidedFacet->sh.anglemapParams.thetaLowerRes + (size_t)((inTheta - collidedFacet->sh.anglemapParams.thetaLimit)
				/ (PI / 2.0 - collidedFacet->sh.anglemapParams.thetaLimit)*(double)collidedFacet->sh.anglemapParams.thetaHigherRes);
		}
		else {
			countTheta = false;
		}
	}
	if (countTheta) {
		size_t phiIndex = (size_t)((inPhi + 3.1415926) / (2.0*PI)*(double)collidedFacet->sh.anglemapParams.phiWidth); //Phi: -PI..PI , and shifting by a number slightly smaller than PI to store on interval [0,2PI[
		collidedFacet->angleMap.pdf[thetaIndex*collidedFacet->sh.anglemapParams.phiWidth + phiIndex]++;
	}
}

void UpdateVelocity(SubprocessFacet *collidedFacet) {
	if (collidedFacet->sh.accomodationFactor > 0.9999) { //speedup for the most common case: perfect thermalization
		if (sHandle->useMaxwellDistribution) sHandle->velocityCurrentParticle = GenerateRandomVelocity(collidedFacet->sh.CDFid);
		else sHandle->velocityCurrentParticle = 145.469*sqrt(collidedFacet->sh.temperature / sHandle->gasMass);
	}
	else {
		double oldSpeed2 = pow(sHandle->velocityCurrentParticle, 2);
		double newSpeed2;
		if (sHandle->useMaxwellDistribution) newSpeed2 = pow(GenerateRandomVelocity(collidedFacet->sh.CDFid), 2);
		else newSpeed2 = /*145.469*/ 29369.939*(collidedFacet->sh.temperature / sHandle->gasMass);
		//sqrt(29369)=171.3766= sqrt(8*R*1000/PI)*3PI/8, that is, the constant part of the v_avg=sqrt(8RT/PI/m/0.001)) found in literature, multiplied by
		//the corrective factor of 3PI/8 that accounts for moving from volumetric speed distribution to wall collision speed distribution
		sHandle->velocityCurrentParticle = sqrt(oldSpeed2 + (newSpeed2 - oldSpeed2)*collidedFacet->sh.accomodationFactor);
	}
}

double GenerateRandomVelocity(int CDFId) {
	//return FastLookupY(rnd(),sHandle->CDFs[CDFId],false);
	double r = rnd();
	double v = InterpolateX(r, sHandle->CDFs[CDFId], false, true); //Allow extrapolate
	return v;
}

double GenerateDesorptionTime(SubprocessFacet *src) {
	if (src->sh.outgassing_paramId >= 0) { //time-dependent desorption
		return InterpolateX(rnd()*sHandle->IDs[src->sh.IDid].back().second, sHandle->IDs[src->sh.IDid], false, true); //allow extrapolate
	}
	else {
		return rnd()*sHandle->latestMoment; //continous desorption between 0 and latestMoment
	}
}

double GetStickingAt(SubprocessFacet *f, double time) {
	if (f->sh.sticking_paramId == -1) //constant sticking
		return f->sh.sticking;
	else return sHandle->parameters[f->sh.sticking_paramId].InterpolateY(time,false);
}

double GetOpacityAt(SubprocessFacet *f, double time) {
	if (f->sh.opacity_paramId == -1) //constant sticking
		return f->sh.opacity;
	else return sHandle->parameters[f->sh.opacity_paramId].InterpolateY(time, false);
}

void TreatMovingFacet() {
	Vector3d localVelocityToAdd;
	if (sHandle->motionType == 1) {
		localVelocityToAdd = sHandle->motionVector2;
	}
	else if (sHandle->motionType == 2) {
		Vector3d distanceVector = 0.01*(sHandle->pPos - sHandle->motionVector1); //distance from base, with cm->m conversion
		localVelocityToAdd = CrossProduct(sHandle->motionVector2, distanceVector);
	}
	Vector3d oldVelocity, newVelocity;
	oldVelocity = sHandle->pDir*sHandle->velocityCurrentParticle;
	newVelocity = oldVelocity + localVelocityToAdd;
	sHandle->pDir = newVelocity.Normalized();
	sHandle->velocityCurrentParticle = newVelocity.Norme();
}

void IncreaseFacetCounter(SubprocessFacet *f, double time, size_t hit, size_t desorb, size_t absorb, double sum_1_per_v, double sum_v_ort) {
	size_t nbMoments = sHandle->moments.size();
	for (size_t m = 0; m <= nbMoments; m++) {
		if (m == 0 || abs(time - sHandle->moments[m - 1]) < sHandle->timeWindowSize / 2.0) {
			f->counter[m].hit.nbMCHit += hit;
			double hitEquiv = static_cast<double>(hit)*sHandle->oriRatio;
			f->counter[m].hit.nbHitEquiv += hitEquiv;
			f->counter[m].hit.nbDesorbed += desorb;
			f->counter[m].hit.nbAbsEquiv += static_cast<double>(absorb)*sHandle->oriRatio;
			f->counter[m].hit.sum_1_per_ort_velocity += sHandle->oriRatio * sum_1_per_v;
			f->counter[m].hit.sum_v_ort += sHandle->oriRatio * sum_v_ort;
			f->counter[m].hit.sum_1_per_velocity += (hitEquiv + static_cast<double>(desorb)) / sHandle->velocityCurrentParticle;
		}
	}
}

void SubprocessFacet::ResetCounter() {
	memset(&counter[0], 0, counter.size() * sizeof(counter[0]));
}

void SubprocessFacet::ResizeCounter(size_t nbMoments) {
	//FacetHitBuffer zeroes;memset(&zeroes, 0, sizeof(zeroes)); //A new zero-value vector
	counter = std::vector<FacetHitBuffer>(nbMoments + 1); //Reserve a counter for each moment, plus an extra for const. flow
	//std::fill(counter.begin(), counter.end(), zeroes); //Initialize each moment with 0 values
	memset(&counter[0], 0, counter.size() * sizeof(counter[0]));
}

void SubprocessFacet::RegisterTransparentPass()
{
	double directionFactor = abs(Dot(sHandle->pDir, this->sh.N));
	IncreaseFacetCounter(this, sHandle->flightTimeCurrentParticle + this->colDist / 100.0 / sHandle->velocityCurrentParticle, 1, 0, 0, 2.0 / (sHandle->velocityCurrentParticle*directionFactor), 2.0*(sHandle->useMaxwellDistribution ? 1.0 : 1.1781)*sHandle->velocityCurrentParticle*directionFactor);

	this->hitted = true;
	if (this->texture && this->sh.countTrans) {
		RecordHitOnTexture(this, sHandle->flightTimeCurrentParticle + this->colDist / 100.0 / sHandle->velocityCurrentParticle,
			true, 2.0, 2.0);
	}
	if (this->direction && this->sh.countDirection) {
		RecordDirectionVector(this, sHandle->flightTimeCurrentParticle + this->colDist / 100.0 / sHandle->velocityCurrentParticle);
	}
	LogHit(this);
	ProfileFacet(this, sHandle->flightTimeCurrentParticle + this->colDist / 100.0 / sHandle->velocityCurrentParticle,
		true, 2.0, 2.0);
	if (this->sh.anglemapParams.record) RecordAngleMap(this);
}
