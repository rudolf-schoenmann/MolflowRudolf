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

extern Simulation *sHandle; //delcared in molflowSub.cpp

// Compute area of all the desorption facet

void CalcTotalOutgassing() {
	//float scale_precomputed;

	// Update texture increment for MC
	//scale_precomputed=(float)(40.0/(sqrt(8.0*8.31/(PI*sHandle->wp.gasMass*0.001))));
	for (size_t j = 0; j < sHandle->sh.nbSuper; j++) {
		for (SubprocessFacet& f : sHandle->structures[j].facets) {
			if (f.sh.is2sided) {
				f.fullSizeInc *= 0.5;
				for (auto& inc : f.textureCellIncrements)
					inc *= 0.5;
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
//	sHandle->currentParticle.direction.x = u*U.x + v*V.x + n*N.x;
//	sHandle->currentParticle.direction.y = u*U.y + v*V.y + n*N.y;
//	sHandle->currentParticle.direction.z = u*U.z + v*V.z + n*N.z;
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
//	double u = Dot(sHandle->currentParticle.direction, iFacet->sh.nU);
//	double v = Dot(sHandle->currentParticle.direction, iFacet->sh.nV);
//	double n = Dot(sHandle->currentParticle.direction, iFacet->sh.N);
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

void UpdateMCHits(Dataport *dpHit, int prIdx, size_t nbMoments, DWORD timeout) {

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
	gHits->globalHits.hit.nbMCHit += sHandle->tmpGlobalResult.globalHits.hit.nbMCHit;
	gHits->globalHits.hit.nbHitEquiv += sHandle->tmpGlobalResult.globalHits.hit.nbHitEquiv;
	gHits->globalHits.hit.nbAbsEquiv += sHandle->tmpGlobalResult.globalHits.hit.nbAbsEquiv;
	gHits->globalHits.hit.nbDesorbed += sHandle->tmpGlobalResult.globalHits.hit.nbDesorbed;
	gHits->distTraveled_total += sHandle->tmpGlobalResult.distTraveled_total;
	gHits->distTraveledTotal_fullHitsOnly += sHandle->tmpGlobalResult.distTraveledTotal_fullHitsOnly;

	//Memorize current limits, then do a min/max search
	for (i = 0; i < 3; i++) {
		texture_limits_old[i] = gHits->texture_limits[i];
		gHits->texture_limits[i].min.all = gHits->texture_limits[i].min.moments_only = HITMAX;
		gHits->texture_limits[i].max.all = gHits->texture_limits[i].max.moments_only = 0;
	}

	//sHandle->wp.sMode = MC_MODE;
	//for(i=0;i<BOUNCEMAX;i++) gHits->wallHits[i] += sHandle->wallHits[i];

	// Leak
	for (size_t leakIndex = 0; leakIndex < sHandle->tmpGlobalResult.leakCacheSize; leakIndex++)
		gHits->leakCache[(leakIndex + gHits->lastLeakIndex) % LEAKCACHESIZE] = sHandle->tmpGlobalResult.leakCache[leakIndex];
	gHits->nbLeakTotal += sHandle->tmpGlobalResult.nbLeakTotal;
	gHits->lastLeakIndex = (gHits->lastLeakIndex + sHandle->tmpGlobalResult.leakCacheSize) % LEAKCACHESIZE;
	gHits->leakCacheSize = Min(LEAKCACHESIZE, gHits->leakCacheSize + sHandle->tmpGlobalResult.leakCacheSize);

	// HHit (Only prIdx 0)
	if (prIdx == 0) {
		for (size_t hitIndex = 0; hitIndex < sHandle->tmpGlobalResult.hitCacheSize; hitIndex++)
			gHits->hitCache[(hitIndex + gHits->lastHitIndex) % HITCACHESIZE] = sHandle->tmpGlobalResult.hitCache[hitIndex];

		if (sHandle->tmpGlobalResult.hitCacheSize > 0) {
			gHits->lastHitIndex = (gHits->lastHitIndex + sHandle->tmpGlobalResult.hitCacheSize) % HITCACHESIZE;
			gHits->hitCache[gHits->lastHitIndex].type = HIT_LAST; //Penup (border between blocks of consecutive hits in the hit cache)
			gHits->hitCacheSize = Min(HITCACHESIZE, gHits->hitCacheSize + sHandle->tmpGlobalResult.hitCacheSize);
		}
	}

	//Global histograms
	
		for (int m = 0; m < (1 + nbMoments); m++) {
			BYTE *histCurrentMoment = buffer + sizeof(GlobalHitBuffer) + m * sHandle->wp.globalHistogramParams.GetDataSize();
			if (sHandle->wp.globalHistogramParams.recordBounce) {
				double* nbHitsHistogram = (double*)histCurrentMoment;
				for (size_t i = 0; i < sHandle->wp.globalHistogramParams.GetBounceHistogramSize(); i++) {
					nbHitsHistogram[i] += sHandle->tmpGlobalHistograms[m].nbHitsHistogram[i];
				}
			}
			if (sHandle->wp.globalHistogramParams.recordDistance) {
				double* distanceHistogram = (double*)(histCurrentMoment + sHandle->wp.globalHistogramParams.GetBouncesDataSize());
				for (size_t i = 0; i < (sHandle->wp.globalHistogramParams.GetDistanceHistogramSize()); i++) {
					distanceHistogram[i] += sHandle->tmpGlobalHistograms[m].distanceHistogram[i];
				}
			}
			if (sHandle->wp.globalHistogramParams.recordTime) {
				double* timeHistogram = (double*)(histCurrentMoment + sHandle->wp.globalHistogramParams.GetBouncesDataSize() + sHandle->wp.globalHistogramParams.GetDistanceDataSize());
				for (size_t i = 0; i < (sHandle->wp.globalHistogramParams.GetTimeHistogramSize()); i++) {
					timeHistogram[i] += sHandle->tmpGlobalHistograms[m].timeHistogram[i];
				}
			}
		}
	

	size_t facetHitsSize = (1 + nbMoments) * sizeof(FacetHitBuffer);
	// Facets
	for (s = 0; s < sHandle->sh.nbSuper; s++) {
		for (SubprocessFacet& f : sHandle->structures[s].facets) {
			if (f.hitted) {

				for (int m = 0; m < (1 + nbMoments); m++) {
					FacetHitBuffer *facetHitBuffer = (FacetHitBuffer *)(buffer + f.sh.hitOffset + m * sizeof(FacetHitBuffer));
					facetHitBuffer->hit.nbAbsEquiv += f.tmpCounter[m].hit.nbAbsEquiv;
					facetHitBuffer->hit.nbDesorbed += f.tmpCounter[m].hit.nbDesorbed;
					facetHitBuffer->hit.nbMCHit += f.tmpCounter[m].hit.nbMCHit;
					facetHitBuffer->hit.nbHitEquiv += f.tmpCounter[m].hit.nbHitEquiv;
					facetHitBuffer->hit.sum_1_per_ort_velocity += f.tmpCounter[m].hit.sum_1_per_ort_velocity;
					facetHitBuffer->hit.sum_v_ort += f.tmpCounter[m].hit.sum_v_ort;
					facetHitBuffer->hit.sum_1_per_velocity += f.tmpCounter[m].hit.sum_1_per_velocity;
				}

				if (f.sh.isProfile) {
					for (int m = 0; m < (1 + nbMoments); m++) {
						ProfileSlice *shProfile = (ProfileSlice *)(buffer + f.sh.hitOffset + facetHitsSize + m * f.profileSize);
						for (j = 0; j < PROFILE_SIZE; j++) {
							shProfile[j] += f.profile[m][j];
						}
					}
				}

				if (f.sh.isTextured) {
					for (int m = 0; m < (1 + nbMoments); m++) {
						TextureCell *shTexture = (TextureCell *)(buffer + (f.sh.hitOffset + facetHitsSize + f.profileSize*(1 + nbMoments) + m * f.textureSize));
						//double dCoef = gHits->globalHits.hit.nbDesorbed * 1E4 * sHandle->wp.gasMass / 1000 / 6E23 * MAGIC_CORRECTION_FACTOR;  //1E4 is conversion from m2 to cm2
						double timeCorrection = m == 0 ? sHandle->wp.finalOutgassingRate : (sHandle->wp.totalDesorbedMolecules) / sHandle->wp.timeWindowSize;
						//Timecorrection is required to compare constant flow texture values with moment values (for autoscaling)

						for (y = 0; y < f.sh.texHeight; y++) {
							for (x = 0; x < f.sh.texWidth; x++) {
								size_t add = x + y * f.sh.texWidth;

								//Add temporary hit counts
								shTexture[add] += f.texture[m][add];

								double val[3];  //pre-calculated autoscaling values (Pressure, imp.rate, density)

								val[0] = shTexture[add].sum_v_ort_per_area*timeCorrection; //pressure without dCoef_pressure
								val[1] = shTexture[add].countEquiv*f.textureCellIncrements[add] * timeCorrection; //imp.rate without dCoef
								val[2] = f.textureCellIncrements[add] * shTexture[add].sum_1_per_ort_velocity* timeCorrection; //particle density without dCoef

								//Global autoscale
								for (int v = 0; v < 3; v++) {
									if (val[v] > gHits->texture_limits[v].max.all && f.largeEnough[add])
										gHits->texture_limits[v].max.all = val[v];

									if (val[v] > 0.0 && val[v] < gHits->texture_limits[v].min.all && f.largeEnough[add])
										gHits->texture_limits[v].min.all = val[v];

									//Autoscale ignoring constant flow (moments only)
									if (m != 0) {
										if (val[v] > gHits->texture_limits[v].max.moments_only && f.largeEnough[add])
											gHits->texture_limits[v].max.moments_only = val[v];

										if (val[v] > 0.0 && val[v] < gHits->texture_limits[v].min.moments_only && f.largeEnough[add])
											gHits->texture_limits[v].min.moments_only = val[v];
									}
								}
							}
						}
					}
				}

				if (f.sh.countDirection) {
					for (int m = 0; m < (1 + nbMoments); m++) {
						DirectionCell *shDir = (DirectionCell *)(buffer + (f.sh.hitOffset + facetHitsSize + f.profileSize*(1 + nbMoments) + f.textureSize*(1 + nbMoments) + f.directionSize*m));
						for (y = 0; y < f.sh.texHeight; y++) {
							for (x = 0; x < f.sh.texWidth; x++) {
								size_t add = x + y * f.sh.texWidth;
								shDir[add].dir.x += f.direction[m][add].dir.x;
								shDir[add].dir.y += f.direction[m][add].dir.y;
								shDir[add].dir.z += f.direction[m][add].dir.z;
								//shDir[add].sumSpeed += f.direction[m][add].sumSpeed;
								shDir[add].count += f.direction[m][add].count;
							}
						}
					}
				}

				if (f.sh.anglemapParams.record) {
					size_t *shAngleMap = (size_t *)(buffer + f.sh.hitOffset + facetHitsSize + f.profileSize*(1 + nbMoments) + f.textureSize*(1 + nbMoments) + f.directionSize*(1 + nbMoments));
					for (y = 0; y < (f.sh.anglemapParams.thetaLowerRes + f.sh.anglemapParams.thetaHigherRes); y++) {
						for (x = 0; x < f.sh.anglemapParams.phiWidth; x++) {
							size_t add = x + y * f.sh.anglemapParams.phiWidth;
							shAngleMap[add] += f.angleMap.pdf[add];
						}
					}
				}

				//Facet histograms
				
					for (int m = 0; m < (1 + nbMoments); m++) {
						BYTE *histCurrentMoment = buffer + f.sh.hitOffset + facetHitsSize + f.profileSize*(1 + nbMoments) + f.textureSize*(1 + nbMoments) + f.directionSize*(1 + nbMoments) + f.sh.anglemapParams.GetRecordedDataSize() + m * f.sh.facetHistogramParams.GetDataSize();
						if (f.sh.facetHistogramParams.recordBounce) {
							double* nbHitsHistogram = (double*)histCurrentMoment;
							for (size_t i = 0; i < f.sh.facetHistogramParams.GetBounceHistogramSize(); i++) {
								nbHitsHistogram[i] += f.tmpHistograms[m].nbHitsHistogram[i];
							}
						}
						if (f.sh.facetHistogramParams.recordDistance) {
							double* distanceHistogram = (double*)(histCurrentMoment + f.sh.facetHistogramParams.GetBouncesDataSize());
							for (size_t i = 0; i < (f.sh.facetHistogramParams.GetDistanceHistogramSize()); i++) {
								distanceHistogram[i] += f.tmpHistograms[m].distanceHistogram[i];
							}
						}
						if (f.sh.facetHistogramParams.recordTime) {
							double* timeHistogram = (double*)(histCurrentMoment + f.sh.facetHistogramParams.GetBouncesDataSize() + f.sh.facetHistogramParams.GetDistanceDataSize());
							for (size_t i = 0; i < (f.sh.facetHistogramParams.GetTimeHistogramSize()); i++) {
								timeHistogram[i] += f.tmpHistograms[m].timeHistogram[i];
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
#ifdef _DEBUG
		double t0, t1;
		t0 = GetTick();
#endif
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
		destIndex = sHandle->currentParticle.teleportedFrom;
		if (destIndex == -1) {
			/*char err[128];
			sprintf(err, "Facet %d tried to teleport to the facet where the particle came from, but there is no such facet.", iFacet->globalId + 1);
			SetErrorSub(err);*/
			RecordHit(HIT_REF);
			sHandle->currentParticle.lastHitFacet = iFacet;
			return; //LEAK
		}
	}
	else destIndex = iFacet->sh.teleportDest - 1;

	//Look in which superstructure is the destination facet:
	for (size_t i = 0; i < sHandle->sh.nbSuper && (!found); i++) {
		for (size_t j = 0; j < sHandle->structures[i].facets.size() && (!found); j++) {
			if (destIndex == sHandle->structures[i].facets[j].globalId) {
				destination = &(sHandle->structures[i].facets[j]);
				if (destination->sh.superIdx != -1) {
					sHandle->currentParticle.structureId = destination->sh.superIdx; //change current superstructure, unless the target is a universal facet
				}
				sHandle->currentParticle.teleportedFrom = (int)iFacet->globalId; //memorize where the particle came from
				found = true;
			}
		}
	}
	if (!found) {
		/*char err[128];
		sprintf(err, "Teleport destination of facet %d not found (facet %d does not exist)", iFacet->globalId + 1, iFacet->sh.teleportDest);
		SetErrorSub(err);*/
		RecordHit(HIT_REF);
		sHandle->currentParticle.lastHitFacet = iFacet;
		return; //LEAK
	}
	// Count this hit as a transparent pass
	RecordHit(HIT_TELEPORTSOURCE);
	if (/*iFacet->texture && */iFacet->sh.countTrans) RecordHitOnTexture(iFacet, sHandle->currentParticle.flightTime, true, 2.0, 2.0);
	if (/*iFacet->direction && */iFacet->sh.countDirection) RecordDirectionVector(iFacet, sHandle->currentParticle.flightTime);
	ProfileFacet(iFacet, sHandle->currentParticle.flightTime, true, 2.0, 2.0);
	LogHit(iFacet);
	if (iFacet->sh.anglemapParams.record) RecordAngleMap(iFacet);

	// Relaunch particle from new facet
	auto[inTheta, inPhi] = CartesianToPolar(sHandle->currentParticle.direction, iFacet->sh.nU, iFacet->sh.nV, iFacet->sh.N);
	PolarToCartesian(destination, inTheta, inPhi, false);
	// Move particle to teleport destination point
	double u = iFacet->colU;
	double v = iFacet->colV;
	sHandle->currentParticle.position = destination->sh.O + u * destination->sh.U + v * destination->sh.V;
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
				sHandle->currentParticle.position = destination->sh.O + u * destination->sh.U + v * destination->sh.V;
				RecordHit(HIT_DES);
			}
		}
		nbTry++;
	}

	sHandle->currentParticle.lastHitFacet = destination;

	//Count hits on teleport facets
	/*iFacet->sh.tmpCounter.hit.nbAbsEquiv++;
	destination->sh.tmpCounter.hit.nbDesorbed++;*/

	double ortVelocity = sHandle->currentParticle.velocity*abs(Dot(sHandle->currentParticle.direction, iFacet->sh.N));
	//We count a teleport as a local hit, but not as a global one since that would affect the MFP calculation
	/*iFacet->sh.tmpCounter.hit.nbMCHit++;
	iFacet->sh.tmpCounter.hit.sum_1_per_ort_velocity += 2.0 / ortVelocity;
	iFacet->sh.tmpCounter.hit.sum_v_ort += 2.0*(sHandle->wp.useMaxwellDistribution ? 1.0 : 1.1781)*ortVelocity;*/
	IncreaseFacetCounter(iFacet, sHandle->currentParticle.flightTime, 1, 0, 0, 2.0 / ortVelocity, 2.0*(sHandle->wp.useMaxwellDistribution ? 1.0 : 1.1781)*ortVelocity);
	iFacet->hitted = true;
	/*destination->sh.tmpCounter.hit.sum_1_per_ort_velocity += 2.0 / sHandle->currentParticle.velocity;
	destination->sh.tmpCounter.hit.sum_v_ort += sHandle->currentParticle.velocity*abs(DOT3(
	sHandle->currentParticle.direction.x, sHandle->currentParticle.direction.y, sHandle->currentParticle.direction.z,
	destination->sh.N.x, destination->sh.N.y, destination->sh.N.z));*/
}

// Perform nbStep simulation steps (a step is a bounce)

bool SimulationMCStep(size_t nbStep) {

	// Perform simulation steps
	for (size_t i = 0; i < nbStep; i++) {

		//Prepare output values
		auto[found, collidedFacet, d] = Intersect(sHandle, sHandle->currentParticle.position, sHandle->currentParticle.direction);

		if (found) {

			// Move particle to intersection point
			sHandle->currentParticle.position = sHandle->currentParticle.position + d * sHandle->currentParticle.direction;
			//sHandle->currentParticle.distanceTraveled += d;

			double lastFLightTime = sHandle->currentParticle.flightTime; //memorize for partial hits
			sHandle->currentParticle.flightTime += d / 100.0 / sHandle->currentParticle.velocity; //conversion from cm to m

			if ((!sHandle->wp.calcConstantFlow && (sHandle->currentParticle.flightTime > sHandle->wp.latestMoment))
				|| (sHandle->wp.enableDecay && (sHandle->currentParticle.expectedDecayMoment < sHandle->currentParticle.flightTime))) {
				//hit time over the measured period - we create a new particle
				//OR particle has decayed
				double remainderFlightPath = sHandle->currentParticle.velocity*100.0*
					Min(sHandle->wp.latestMoment - lastFLightTime, sHandle->currentParticle.expectedDecayMoment - lastFLightTime); //distance until the point in space where the particle decayed
				sHandle->tmpGlobalResult.distTraveled_total += remainderFlightPath * sHandle->currentParticle.oriRatio;
				RecordHit(HIT_LAST);
				//sHandle->distTraveledSinceUpdate += sHandle->currentParticle.distanceTraveled;
				if (!StartFromSource())
					// desorptionLimit reached
					return false;
			}
			else { //hit within measured time, particle still alive
				if (collidedFacet->sh.teleportDest != 0) { //Teleport
					IncreaseDistanceCounters(d * sHandle->currentParticle.oriRatio);
					PerformTeleport(collidedFacet);
				}
				/*else if ((GetOpacityAt(collidedFacet, sHandle->currentParticle.flightTime) < 1.0) && (rnd() > GetOpacityAt(collidedFacet, sHandle->currentParticle.flightTime))) {
					//Transparent pass
					sHandle->tmpGlobalResult.distTraveled_total += d;
					PerformTransparentPass(collidedFacet);
				}*/
				else { //Not teleport
					IncreaseDistanceCounters(d * sHandle->currentParticle.oriRatio);
					double stickingProbability = GetStickingAt(collidedFacet, sHandle->currentParticle.flightTime);
					if (!sHandle->ontheflyParams.lowFluxMode) { //Regular stick or bounce
						if (stickingProbability == 1.0 || ((stickingProbability > 0.0) && (rnd() < (stickingProbability)))) {
							//Absorbed
							RecordAbsorb(collidedFacet);
							//sHandle->distTraveledSinceUpdate += sHandle->currentParticle.distanceTraveled;
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
							double oriRatioBeforeCollision = sHandle->currentParticle.oriRatio; //Local copy
							sHandle->currentParticle.oriRatio *= (stickingProbability); //Sticking part
							RecordAbsorb(collidedFacet);
							sHandle->currentParticle.oriRatio = oriRatioBeforeCollision * (1.0 - stickingProbability); //Reflected part
						}
						else
							sHandle->currentParticle.oriRatio *= (1.0 - stickingProbability);
						if (sHandle->currentParticle.oriRatio > sHandle->ontheflyParams.lowFluxCutoff) {
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
			sHandle->tmpGlobalResult.nbLeakTotal++;
			RecordLeakPos();
			if (!StartFromSource())
				// desorptionLimit reached
				return false;
		}
	}
	return true;
}

void IncreaseDistanceCounters(double distanceIncrement)
{
	sHandle->tmpGlobalResult.distTraveled_total += distanceIncrement;
	sHandle->tmpGlobalResult.distTraveledTotal_fullHitsOnly += distanceIncrement;
	sHandle->currentParticle.distanceTraveled += distanceIncrement;
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
		if (sHandle->totalDesorbed >= sHandle->ontheflyParams.desorptionLimit / sHandle->ontheflyParams.nbProcess) {
			sHandle->currentParticle.lastHitFacet = NULL;
			return false;
		}
	}

	// Select source
	srcRnd = rnd() * sHandle->wp.totalDesorbedMolecules;

	while (!found && j < sHandle->sh.nbSuper) { //Go through superstructures
		i = 0;
		while (!found && i < sHandle->structures[j].facets.size()) { //Go through facets in a structure
			SubprocessFacet& f = sHandle->structures[j].facets[i];
			if (f.sh.desorbType != DES_NONE) { //there is some kind of outgassing
				if (f.sh.useOutgassingFile) { //Using SynRad-generated outgassing map
					if (f.sh.totalOutgassing > 0.0) {
						found = (srcRnd >= sumA) && (srcRnd < (sumA + sHandle->wp.latestMoment * f.sh.totalOutgassing / (1.38E-23*f.sh.temperature)));
						if (found) {
							//look for exact position in map
							double rndRemainder = (srcRnd - sumA) / sHandle->wp.latestMoment*(1.38E-23*f.sh.temperature); //remainder, should be less than f.sh.totalOutgassing
							/*double sumB = 0.0;
							for (w = 0; w < f.sh.outgassingMapWidth && !foundInMap; w++) {
								for (h = 0; h < f.sh.outgassingMapHeight && !foundInMap; h++) {
									double cellOutgassing = f.outgassingMap[h*f.sh.outgassingMapWidth + w];
									if (cellOutgassing > 0.0) {
										foundInMap = (rndRemainder >= sumB) && (rndRemainder < (sumB + cellOutgassing));
										if (foundInMap) mapPositionW = w; mapPositionH = h;
										sumB += cellOutgassing;
									}
								}
							}*/
							double lookupValue = rndRemainder;
							int outgLowerIndex = my_lower_bound(lookupValue, f.outgassingMap); //returns line number AFTER WHICH LINE lookup value resides in ( -1 .. size-2 )
							outgLowerIndex++;
							mapPositionH = (size_t)((double)outgLowerIndex / (double)f.sh.outgassingMapWidth);
							mapPositionW = (size_t)outgLowerIndex - mapPositionH * f.sh.outgassingMapWidth;
							foundInMap = true;
							/*if (!foundInMap) {
								SetErrorSub("Starting point not found in imported desorption map");
								return false;
							}*/
						}
						sumA += sHandle->wp.latestMoment * f.sh.totalOutgassing / (1.38E-23*f.sh.temperature);
					}
				} //end outgassing file block
				else { //constant or time-dependent outgassing
					double facetOutgassing =
						(f.sh.outgassing_paramId >= 0)
						? sHandle->IDs[f.sh.IDid].back().second / (1.38E-23*f.sh.temperature)
						: sHandle->wp.latestMoment*f.sh.outgassing / (1.38E-23*f.sh.temperature);
					found = (srcRnd >= sumA) && (srcRnd < (sumA + facetOutgassing));
					sumA += facetOutgassing;
				} //end constant or time-dependent outgassing block
			} //end 'there is some kind of outgassing'
			if (!found) i++;
			if (f.sh.is2sided) reverse = rnd() > 0.5;
			else reverse = false;
		}
		if (!found) j++;
	}
	if (!found) {
		SetErrorSub("No starting point, aborting");
		return false;
	}
	src = &(sHandle->structures[j].facets[i]);

	sHandle->currentParticle.lastHitFacet = src;
	//sHandle->currentParticle.distanceTraveled = 0.0;  //for mean free path calculations
	//sHandle->currentParticle.flightTime = sHandle->desorptionStartTime + (sHandle->desorptionStopTime - sHandle->desorptionStartTime)*rnd();
	sHandle->currentParticle.flightTime = GenerateDesorptionTime(src);
	if (sHandle->wp.useMaxwellDistribution) sHandle->currentParticle.velocity = GenerateRandomVelocity(src->sh.CDFid);
	else sHandle->currentParticle.velocity = 145.469*sqrt(src->sh.temperature / sHandle->wp.gasMass);  //sqrt(8*R/PI/1000)=145.47
	sHandle->currentParticle.oriRatio = 1.0;
	if (sHandle->wp.enableDecay) { //decaying gas
		sHandle->currentParticle.expectedDecayMoment = sHandle->currentParticle.flightTime + sHandle->wp.halfLife*1.44269*-log(rnd()); //1.44269=1/ln2
		//Exponential distribution PDF: probability of 't' life = 1/TAU*exp(-t/TAU) where TAU = half_life/ln2
		//Exponential distribution CDF: probability of life shorter than 't" = 1-exp(-t/TAU)
		//Equation: rnd()=1-exp(-t/TAU)
		//Solution: t=TAU*-log(1-rnd()) and 1-rnd()=rnd() therefore t=half_life/ln2*-log(rnd())
	}
	else {
		sHandle->currentParticle.expectedDecayMoment = 1e100; //never decay
	}
	//sHandle->temperature = src->sh.temperature; //Thermalize particle
	sHandle->currentParticle.nbBounces = 0;
	sHandle->currentParticle.distanceTraveled = 0;

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
			sHandle->currentParticle.position = src->sh.O + u * src->sh.U + v * src->sh.V;
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
			sHandle->currentParticle.position = src->sh.O + u * src->sh.U + v * src->sh.V;
			src->colU = u;
			src->colV = v;
		}
		else {
			src->colU = 0.5;
			src->colV = 0.5;
			sHandle->currentParticle.position = sHandle->structures[j].facets[i].sh.center;
		}

	}

	if (src->sh.isMoving && sHandle->wp.motionType) RecordHit(HIT_MOVING);
	else RecordHit(HIT_DES); //create blue hit point for created particle

	//See docs/theta_gen.png for further details on angular distribution generation
	switch (src->sh.desorbType) {
	case DES_UNIFORM:
		sHandle->currentParticle.direction = PolarToCartesian(src, acos(rnd()), rnd()*2.0*PI, reverse);
		break;
	case DES_NONE: //for file-based
	case DES_COSINE:
		sHandle->currentParticle.direction = PolarToCartesian(src, acos(sqrt(rnd())), rnd()*2.0*PI, reverse);
		break;
	case DES_COSINE_N:
		sHandle->currentParticle.direction = PolarToCartesian(src, acos(pow(rnd(), 1.0 / (src->sh.desorbTypeN + 1.0))), rnd()*2.0*PI, reverse);
		break;
	case DES_ANGLEMAP:
	{
		auto [theta, thetaLowerIndex, thetaOvershoot] = src->angleMap.GenerateThetaFromAngleMap(src->sh.anglemapParams);
		auto phi = src->angleMap.GeneratePhiFromAngleMap(thetaLowerIndex, thetaOvershoot, src->sh.anglemapParams);
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
		sHandle->currentParticle.direction = PolarToCartesian(src, PI - theta, phi, false); //angle map contains incident angle (between N and source dir) and theta is dir (between N and dest dir)
		_ASSERTE(sHandle->currentParticle.direction.x == sHandle->currentParticle.direction.x);

	}
	}

	// Current structure
	if (src->sh.superIdx == -1) {
		std::ostringstream out;
		out << "Facet " << (src->globalId + 1) << " is in all structures, it shouldn't desorb.";
		SetErrorSub(out.str().c_str());
		return false;
	}
	sHandle->currentParticle.structureId = src->sh.superIdx;
	sHandle->currentParticle.teleportedFrom = -1;

	// Count

	src->hitted = true;
	sHandle->totalDesorbed++;
	sHandle->tmpGlobalResult.globalHits.hit.nbDesorbed++;
	//sHandle->nbPHit = 0;

	if (src->sh.isMoving) {
		TreatMovingFacet();
	}

	double ortVelocity = sHandle->currentParticle.velocity*abs(Dot(sHandle->currentParticle.direction, src->sh.N));
	/*src->sh.tmpCounter.hit.nbDesorbed++;
	src->sh.tmpCounter.hit.sum_1_per_ort_velocity += 2.0 / ortVelocity; //was 2.0 / ortV
	src->sh.tmpCounter.hit.sum_v_ort += (sHandle->wp.useMaxwellDistribution ? 1.0 : 1.1781)*ortVelocity;*/
	IncreaseFacetCounter(src, sHandle->currentParticle.flightTime, 0, 1, 0, 2.0 / ortVelocity, (sHandle->wp.useMaxwellDistribution ? 1.0 : 1.1781)*ortVelocity);
	//Desorption doesn't contribute to angular profiles, nor to angle maps
	ProfileFacet(src, sHandle->currentParticle.flightTime, false, 2.0, 1.0); //was 2.0, 1.0
	LogHit(src);
	if (/*src->texture && */src->sh.countDes) RecordHitOnTexture(src, sHandle->currentParticle.flightTime, true, 2.0, 1.0); //was 2.0, 1.0
	//if (src->direction && src->sh.countDirection) RecordDirectionVector(src, sHandle->currentParticle.flightTime);

	// Reset volatile state
	if (sHandle->hasVolatile) {
		for (auto& s : sHandle->structures) {
			for (auto& f : s.facets) {
				f.ready = true;
			}
		}
	}

	found = false;
	return true;
}

std::tuple<double, int, double> Anglemap::GenerateThetaFromAngleMap(const AnglemapParams& anglemapParams)
{
	double lookupValue = rnd();
	int thetaLowerIndex = my_lower_bound(lookupValue, theta_CDF); //returns line number AFTER WHICH LINE lookup value resides in ( -1 .. size-2 )
	double theta, thetaOvershoot;

	if (thetaLowerIndex == -1) { //first half section
		thetaOvershoot = 0.5 + 0.5 * lookupValue / theta_CDF[0]; //between 0.5 and 1
		theta = GetTheta((double)thetaLowerIndex + 0.5 + thetaOvershoot, anglemapParams); //between 0 and the first section end
		return { theta, thetaLowerIndex, thetaOvershoot };
	}
	else if (thetaLowerIndex == (anglemapParams.thetaLowerRes + anglemapParams.thetaHigherRes - 1)) { //last half section //can this happen?
		thetaOvershoot = 0.5 * (lookupValue - theta_CDF[thetaLowerIndex])
			/ (1.0 - theta_CDF[thetaLowerIndex]); //between 0 and 0.5
		theta = GetTheta((double)thetaLowerIndex + 0.5 + thetaOvershoot, anglemapParams); //between 0 and the first section end
		return { theta, thetaLowerIndex, thetaOvershoot };
	}
	else { //regular section
		if (/*true || */phi_CDFsums[thetaLowerIndex] == phi_CDFsums[thetaLowerIndex + 1]) {
			//The pdf's slope is 0, linear interpolation
			thetaOvershoot = (lookupValue - theta_CDF[thetaLowerIndex]) / (theta_CDF[thetaLowerIndex + 1] - theta_CDF[thetaLowerIndex]);
			theta = GetTheta((double)thetaLowerIndex + 0.5 + thetaOvershoot, anglemapParams);
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
			double a = 0.5 * ((double)(phi_CDFsums[thetaLowerIndex + 1]) - (double)phi_CDFsums[thetaLowerIndex]) / (double)theta_CDFsum / Sqr(thetaStep); //pdf slope at lower index
			double dy = lookupValue - c;

			double dx = (-b + sqrt(Sqr(b) + 4 * a*dy)) / (2 * a); //Since b>=0 it's the + branch of the +- that is valid for us

			thetaOvershoot = dx / thetaStep;
			theta = GetTheta((double)thetaLowerIndex + 0.5 + thetaOvershoot, anglemapParams);
		}
	}
	_ASSERTE(theta == theta);
	return { theta, thetaLowerIndex, thetaOvershoot };
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

		double div, weigh;
		div = ((double)phi_CDFsums[thetaLowerIndex] * (1.0 - thetaOvershoot) + (double)phi_CDFsums[thetaLowerIndex + 1] * thetaOvershoot); // (w1*w3 + w2*w4)
		if (div > 0.0) {
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
		if (/*true ||*/ GetPhipdfValue(thetaIndex, phiLowerIndex, anglemapParams) == GetPhipdfValue(thetaIndex, phiLowerIndex + 1, anglemapParams)) {
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
			double b = GetPhipdfValue(thetaIndex, phiLowerIndex, anglemapParams) / GetPhiCDFSum(thetaIndex, anglemapParams) / phiStep; //pdf value at lower index
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
		return anglemapParams.thetaLimit + (PI / 2.0 - anglemapParams.thetaLimit) * (thetaIndex - (double)anglemapParams.thetaLowerRes) / (double)anglemapParams.thetaHigherRes;
	}
}

double Anglemap::GetPhi(const double & phiIndex, const AnglemapParams & anglemapParams)
//makes phiIndex circular and converts from index to -pi...pi
{
	double width = (double)anglemapParams.phiWidth;
	double correctedIndex = (phiIndex < width) ? phiIndex : phiIndex - width;
	return -PI + 2.0 * PI * correctedIndex / width;
}

double Anglemap::GetPhipdfValue(const double & thetaIndex, const int & phiLowerIndex, const AnglemapParams& anglemapParams)
//phiLowerIndex is circularized
{
	if (thetaIndex < 0.5) {
		return (double)pdf[IDX(phiLowerIndex, anglemapParams.phiWidth)];
	}
	else if (thetaIndex > (double)(anglemapParams.thetaLowerRes + anglemapParams.thetaHigherRes) - 0.5) {
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
	else if (thetaIndex > (double)(anglemapParams.thetaLowerRes + anglemapParams.thetaHigherRes) - 0.5) {
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
	else if (thetaIndex > (double)(anglemapParams.thetaLowerRes + anglemapParams.thetaHigherRes) - 0.5) {
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

	sHandle->tmpGlobalResult.globalHits.hit.nbMCHit++; //global
	sHandle->tmpGlobalResult.globalHits.hit.nbHitEquiv += sHandle->currentParticle.oriRatio;

	// Handle super structure link facet. Can be 
	if (iFacet->sh.superDest) {
		IncreaseFacetCounter(iFacet, sHandle->currentParticle.flightTime, 1, 0, 0, 0, 0);
		sHandle->currentParticle.structureId = iFacet->sh.superDest - 1;
		if (iFacet->sh.isMoving) { //A very special case where link facets can be used as transparent but moving facets
			RecordHit(HIT_MOVING);
			TreatMovingFacet();
		}
		else {
			// Count this hit as a transparent pass
			RecordHit(HIT_TRANS);
		}
		LogHit(iFacet);
		ProfileFacet(iFacet, sHandle->currentParticle.flightTime, true, 2.0, 2.0);
		if (iFacet->sh.anglemapParams.record) RecordAngleMap(iFacet);
		if (/*iFacet->texture &&*/ iFacet->sh.countTrans) RecordHitOnTexture(iFacet, sHandle->currentParticle.flightTime, true, 2.0, 2.0);
		if (/*iFacet->direction &&*/ iFacet->sh.countDirection) RecordDirectionVector(iFacet, sHandle->currentParticle.flightTime);

		return;

	}

	// Handle volatile facet
	if (iFacet->sh.isVolatile) {

		if (iFacet->ready) {
			IncreaseFacetCounter(iFacet, sHandle->currentParticle.flightTime, 0, 0, 1, 0, 0);
			iFacet->ready = false;
			LogHit(iFacet);
			ProfileFacet(iFacet, sHandle->currentParticle.flightTime, true, 2.0, 1.0);
			if (/*iFacet->texture && */iFacet->sh.countAbs) RecordHitOnTexture(iFacet, sHandle->currentParticle.flightTime, true, 2.0, 1.0);
			if (/*iFacet->direction && */iFacet->sh.countDirection) RecordDirectionVector(iFacet, sHandle->currentParticle.flightTime);
		}
		return;

	}

	if (iFacet->sh.is2sided) {
		// We may need to revert normal in case of 2 sided hit
		revert = Dot(sHandle->currentParticle.direction, iFacet->sh.N) > 0.0;
	}

	//Texture/Profile incoming hit


	//Register (orthogonal) velocity
	double ortVelocity = sHandle->currentParticle.velocity*abs(Dot(sHandle->currentParticle.direction, iFacet->sh.N));

	/*iFacet->sh.tmpCounter.hit.nbMCHit++; //hit facet
	iFacet->sh.tmpCounter.hit.sum_1_per_ort_velocity += 1.0 / ortVelocity;
	iFacet->sh.tmpCounter.hit.sum_v_ort += (sHandle->wp.useMaxwellDistribution ? 1.0 : 1.1781)*ortVelocity;*/

	IncreaseFacetCounter(iFacet, sHandle->currentParticle.flightTime, 1, 0, 0, 1.0 / ortVelocity, (sHandle->wp.useMaxwellDistribution ? 1.0 : 1.1781)*ortVelocity);
	sHandle->currentParticle.nbBounces++;
	if (/*iFacet->texture &&*/ iFacet->sh.countRefl) RecordHitOnTexture(iFacet, sHandle->currentParticle.flightTime, true, 1.0, 1.0);
	if (/*iFacet->direction &&*/ iFacet->sh.countDirection) RecordDirectionVector(iFacet, sHandle->currentParticle.flightTime);
	LogHit(iFacet);
	ProfileFacet(iFacet, sHandle->currentParticle.flightTime, true, 1.0, 1.0);
	if (iFacet->sh.anglemapParams.record) RecordAngleMap(iFacet);

	// Relaunch particle
	UpdateVelocity(iFacet);
	//Sojourn time
	if (iFacet->sh.enableSojournTime) {
		double A = exp(-iFacet->sh.sojournE / (8.31*iFacet->sh.temperature));
		sHandle->currentParticle.flightTime += -log(rnd()) / (A*iFacet->sh.sojournFreq);
	}

	if (iFacet->sh.reflection.diffusePart > 0.999999) { //Speedup branch for most common, diffuse case
		sHandle->currentParticle.direction = PolarToCartesian(iFacet, acos(sqrt(rnd())), rnd()*2.0*PI, revert);
	}
	else {
		double reflTypeRnd = rnd();
		if (reflTypeRnd < iFacet->sh.reflection.diffusePart)
		{
			//diffuse reflection
			//See docs/theta_gen.png for further details on angular distribution generation
			sHandle->currentParticle.direction = PolarToCartesian(iFacet, acos(sqrt(rnd())), rnd()*2.0*PI, revert);
		}
		else  if (reflTypeRnd < (iFacet->sh.reflection.diffusePart + iFacet->sh.reflection.specularPart))
		{
			//specular reflection
			auto [inTheta, inPhi] = CartesianToPolar(sHandle->currentParticle.direction, iFacet->sh.nU, iFacet->sh.nV, iFacet->sh.N);
			sHandle->currentParticle.direction = PolarToCartesian(iFacet, PI - inTheta, inPhi, false);

		}
		else {
			//Cos^N reflection
			sHandle->currentParticle.direction = PolarToCartesian(iFacet, acos(pow(rnd(), 1.0 / (iFacet->sh.reflection.cosineExponent + 1.0))), rnd()*2.0*PI, revert);
		}
	}

	if (iFacet->sh.isMoving) {
		TreatMovingFacet();
	}

	//Texture/Profile outgoing particle
	//Register outgoing velocity
	ortVelocity = sHandle->currentParticle.velocity*abs(Dot(sHandle->currentParticle.direction, iFacet->sh.N));

	/*iFacet->sh.tmpCounter.hit.sum_1_per_ort_velocity += 1.0 / ortVelocity;
	iFacet->sh.tmpCounter.hit.sum_v_ort += (sHandle->wp.useMaxwellDistribution ? 1.0 : 1.1781)*ortVelocity;*/
	IncreaseFacetCounter(iFacet, sHandle->currentParticle.flightTime, 0, 0, 0, 1.0 / ortVelocity, (sHandle->wp.useMaxwellDistribution ? 1.0 : 1.1781)*ortVelocity);
	if (/*iFacet->texture &&*/ iFacet->sh.countRefl) RecordHitOnTexture(iFacet, sHandle->currentParticle.flightTime, false, 1.0, 1.0); //count again for outward velocity
	ProfileFacet(iFacet, sHandle->currentParticle.flightTime, false, 1.0, 1.0);
	//no direction count on outgoing, neither angle map

	if (iFacet->sh.isMoving && sHandle->wp.motionType) RecordHit(HIT_MOVING);
	else RecordHit(HIT_REF);
	sHandle->currentParticle.lastHitFacet = iFacet;
	//sHandle->nbPHit++;
}

void PerformTransparentPass(SubprocessFacet *iFacet) { //disabled, caused finding hits with the same facet
	/*double directionFactor = abs(DOT3(
		sHandle->currentParticle.direction.x, sHandle->currentParticle.direction.y, sHandle->currentParticle.direction.z,
		iFacet->sh.N.x, iFacet->sh.N.y, iFacet->sh.N.z));
	iFacet->sh.tmpCounter.hit.nbMCHit++;
	iFacet->sh.tmpCounter.hit.sum_1_per_ort_velocity += 2.0 / (sHandle->currentParticle.velocity*directionFactor);
	iFacet->sh.tmpCounter.hit.sum_v_ort += 2.0*(sHandle->wp.useMaxwellDistribution ? 1.0 : 1.1781)*sHandle->currentParticle.velocity*directionFactor;
	iFacet->hitted = true;
	if (iFacet->texture && iFacet->sh.countTrans) RecordHitOnTexture(iFacet, sHandle->currentParticle.flightTime + iFacet->colDist / 100.0 / sHandle->currentParticle.velocity,
		true, 2.0, 2.0);
	if (iFacet->direction && iFacet->sh.countDirection) RecordDirectionVector(iFacet, sHandle->currentParticle.flightTime + iFacet->colDist / 100.0 / sHandle->currentParticle.velocity);
	ProfileFacet(iFacet, sHandle->currentParticle.flightTime + iFacet->colDist / 100.0 / sHandle->currentParticle.velocity,
		true, 2.0, 2.0);
	RecordHit(HIT_TRANS);
	sHandle->lastHit = iFacet;*/
}

void RecordAbsorb(SubprocessFacet *iFacet) {
	sHandle->tmpGlobalResult.globalHits.hit.nbMCHit++; //global	
	sHandle->tmpGlobalResult.globalHits.hit.nbHitEquiv += sHandle->currentParticle.oriRatio;
	sHandle->tmpGlobalResult.globalHits.hit.nbAbsEquiv += sHandle->currentParticle.oriRatio;

	RecordHistograms(iFacet);

	RecordHit(HIT_ABS);
	double ortVelocity = sHandle->currentParticle.velocity*abs(Dot(sHandle->currentParticle.direction, iFacet->sh.N));
	IncreaseFacetCounter(iFacet, sHandle->currentParticle.flightTime, 1, 0, 1, 2.0 / ortVelocity, (sHandle->wp.useMaxwellDistribution ? 1.0 : 1.1781)*ortVelocity);
	LogHit(iFacet);
	ProfileFacet(iFacet, sHandle->currentParticle.flightTime, true, 2.0, 1.0); //was 2.0, 1.0
	if (iFacet->sh.anglemapParams.record) RecordAngleMap(iFacet);
	if (/*iFacet->texture &&*/ iFacet->sh.countAbs) RecordHitOnTexture(iFacet, sHandle->currentParticle.flightTime, true, 2.0, 1.0); //was 2.0, 1.0
	if (/*iFacet->direction &&*/ iFacet->sh.countDirection) RecordDirectionVector(iFacet, sHandle->currentParticle.flightTime);
}

void RecordHistograms(SubprocessFacet * iFacet)
{
	//Record in global and facet histograms
	for (size_t m = 0; m <= sHandle->moments.size(); m++) {
		if (m == 0 || abs(sHandle->currentParticle.flightTime - sHandle->moments[m - 1]) < sHandle->wp.timeWindowSize / 2.0) {
			size_t binIndex;
			if (sHandle->wp.globalHistogramParams.recordBounce) {
				binIndex = Min(sHandle->currentParticle.nbBounces / sHandle->wp.globalHistogramParams.nbBounceBinsize, sHandle->wp.globalHistogramParams.GetBounceHistogramSize() - 1);
				sHandle->tmpGlobalHistograms[m].nbHitsHistogram[binIndex] += sHandle->currentParticle.oriRatio;
			}
			if (sHandle->wp.globalHistogramParams.recordDistance) {
				binIndex = Min(static_cast<size_t>(sHandle->currentParticle.distanceTraveled / sHandle->wp.globalHistogramParams.distanceBinsize), sHandle->wp.globalHistogramParams.GetDistanceHistogramSize() - 1);
				sHandle->tmpGlobalHistograms[m].distanceHistogram[binIndex] += sHandle->currentParticle.oriRatio;
			}
			if (sHandle->wp.globalHistogramParams.recordTime) {
				binIndex = Min(static_cast<size_t>(sHandle->currentParticle.flightTime / sHandle->wp.globalHistogramParams.timeBinsize), sHandle->wp.globalHistogramParams.GetTimeHistogramSize() - 1);
				sHandle->tmpGlobalHistograms[m].timeHistogram[binIndex] += sHandle->currentParticle.oriRatio;
			}
			if (iFacet->sh.facetHistogramParams.recordBounce) {
				binIndex = Min(sHandle->currentParticle.nbBounces / iFacet->sh.facetHistogramParams.nbBounceBinsize, iFacet->sh.facetHistogramParams.GetBounceHistogramSize() - 1);
				iFacet->tmpHistograms[m].nbHitsHistogram[binIndex] += sHandle->currentParticle.oriRatio;
			}
			if (iFacet->sh.facetHistogramParams.recordDistance) {
				binIndex = Min(static_cast<size_t>(sHandle->currentParticle.distanceTraveled / iFacet->sh.facetHistogramParams.distanceBinsize), iFacet->sh.facetHistogramParams.GetDistanceHistogramSize() - 1);
				iFacet->tmpHistograms[m].distanceHistogram[binIndex] += sHandle->currentParticle.oriRatio;
			}
			if (iFacet->sh.facetHistogramParams.recordTime) {
				binIndex = Min(static_cast<size_t>(sHandle->currentParticle.flightTime / iFacet->sh.facetHistogramParams.timeBinsize), iFacet->sh.facetHistogramParams.GetTimeHistogramSize() - 1);
				iFacet->tmpHistograms[m].timeHistogram[binIndex] += sHandle->currentParticle.oriRatio;
			}
		}
	}
}

void RecordHitOnTexture(SubprocessFacet *f, double time, bool countHit, double velocity_factor, double ortSpeedFactor) {

	size_t tu = (size_t)(f->colU * f->sh.texWidthD);
	size_t tv = (size_t)(f->colV * f->sh.texHeightD);
	size_t add = tu + tv * (f->sh.texWidth);
	double ortVelocity = (sHandle->wp.useMaxwellDistribution ? 1.0 : 1.1781)*sHandle->currentParticle.velocity*abs(Dot(sHandle->currentParticle.direction, f->sh.N)); //surface-orthogonal velocity component

	for (size_t m = 0; m <= sHandle->moments.size(); m++)
		if (m == 0 || abs(time - sHandle->moments[m - 1]) < sHandle->wp.timeWindowSize / 2.0) {
			if (countHit) f->texture[m][add].countEquiv += sHandle->currentParticle.oriRatio;
			f->texture[m][add].sum_1_per_ort_velocity += sHandle->currentParticle.oriRatio * velocity_factor / ortVelocity;
			f->texture[m][add].sum_v_ort_per_area += sHandle->currentParticle.oriRatio * ortSpeedFactor*ortVelocity*f->textureCellIncrements[add]; // sum ortho_velocity[m/s] / cell_area[cm2]
		}
}

void RecordDirectionVector(SubprocessFacet *f, double time) {
	size_t tu = (size_t)(f->colU * f->sh.texWidthD);
	size_t tv = (size_t)(f->colV * f->sh.texHeightD);
	size_t add = tu + tv * (f->sh.texWidth);

	for (size_t m = 0; m <= sHandle->moments.size(); m++) {
		if (m == 0 || abs(time - sHandle->moments[m - 1]) < sHandle->wp.timeWindowSize / 2.0) {
			f->direction[m][add].dir = f->direction[m][add].dir + sHandle->currentParticle.oriRatio * sHandle->currentParticle.direction * sHandle->currentParticle.velocity;
			f->direction[m][add].count++;
		}
	}

}

void ProfileFacet(SubprocessFacet *f, double time, bool countHit, double velocity_factor, double ortSpeedFactor) {

	size_t nbMoments = sHandle->moments.size();

	if (countHit && f->sh.profileType == PROFILE_ANGULAR) {
		double dot = Dot(f->sh.N, sHandle->currentParticle.direction);
		double theta = acos(abs(dot));     // Angle to normal (PI/2 => PI)
		size_t pos = (size_t)(theta / (PI / 2)*((double)PROFILE_SIZE)); // To Grad
		Saturate(pos, 0, PROFILE_SIZE - 1);
		for (size_t m = 0; m <= nbMoments; m++) {
			if (m == 0 || abs(time - sHandle->moments[m - 1]) < sHandle->wp.timeWindowSize / 2.0) {
				f->profile[m][pos].countEquiv += sHandle->currentParticle.oriRatio;
			}
		}
	}
	else if (f->sh.profileType == PROFILE_U || f->sh.profileType == PROFILE_V) {
		size_t pos = (size_t)((f->sh.profileType == PROFILE_U ? f->colU : f->colV)*(double)PROFILE_SIZE);
		if (pos >= 0 && pos < PROFILE_SIZE) {
			for (size_t m = 0; m <= nbMoments; m++) {
				if (m == 0 || abs(time - sHandle->moments[m - 1]) < sHandle->wp.timeWindowSize / 2.0) {
					if (countHit) f->profile[m][pos].countEquiv += sHandle->currentParticle.oriRatio;
					double ortVelocity = sHandle->currentParticle.velocity*abs(Dot(f->sh.N, sHandle->currentParticle.direction));
					f->profile[m][pos].sum_1_per_ort_velocity += sHandle->currentParticle.oriRatio * velocity_factor / ortVelocity;
					f->profile[m][pos].sum_v_ort += sHandle->currentParticle.oriRatio * ortSpeedFactor*(sHandle->wp.useMaxwellDistribution ? 1.0 : 1.1781)*ortVelocity;
				}
			}
		}
	}
	else if (countHit && (f->sh.profileType == PROFILE_VELOCITY || f->sh.profileType == PROFILE_ORT_VELOCITY || f->sh.profileType == PROFILE_TAN_VELOCITY)) {
		double dot;
		if (f->sh.profileType == PROFILE_VELOCITY) {
			dot = 1.0;
		}
		else if (f->sh.profileType == PROFILE_ORT_VELOCITY) {
			dot = abs(Dot(f->sh.N, sHandle->currentParticle.direction));  //cos(theta) as "dot" value
		}
		else { //Tangential
			dot = sqrt(1 - Sqr(abs(Dot(f->sh.N, sHandle->currentParticle.direction))));  //tangential
		}
		size_t pos = (size_t)(dot*sHandle->currentParticle.velocity / f->sh.maxSpeed*(double)PROFILE_SIZE); //"dot" default value is 1.0
		if (pos >= 0 && pos < PROFILE_SIZE) {
			for (size_t m = 0; m <= nbMoments; m++) {
				if (m == 0 || abs(time - sHandle->moments[m - 1]) < sHandle->wp.timeWindowSize / 2.0) {
					f->profile[m][pos].countEquiv += sHandle->currentParticle.oriRatio;
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
		std::tie(log.hitTheta, log.hitPhi) = CartesianToPolar(sHandle->currentParticle.direction, f->sh.nU, f->sh.nV, f->sh.N);
		log.oriRatio = sHandle->currentParticle.oriRatio;
		log.particleDecayMoment = sHandle->currentParticle.expectedDecayMoment;
		log.time = sHandle->currentParticle.flightTime;
		log.velocity = sHandle->currentParticle.velocity;
		sHandle->tmpParticleLog.push_back(log);
	}
}

void RecordAngleMap(SubprocessFacet* collidedFacet) {
	auto[inTheta, inPhi] = CartesianToPolar(sHandle->currentParticle.direction, collidedFacet->sh.nU, collidedFacet->sh.nV, collidedFacet->sh.N);
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
		if (sHandle->wp.useMaxwellDistribution) sHandle->currentParticle.velocity = GenerateRandomVelocity(collidedFacet->sh.CDFid);
		else sHandle->currentParticle.velocity = 145.469*sqrt(collidedFacet->sh.temperature / sHandle->wp.gasMass);
	}
	else {
		double oldSpeed2 = pow(sHandle->currentParticle.velocity, 2);
		double newSpeed2;
		if (sHandle->wp.useMaxwellDistribution) newSpeed2 = pow(GenerateRandomVelocity(collidedFacet->sh.CDFid), 2);
		else newSpeed2 = /*145.469*/ 29369.939*(collidedFacet->sh.temperature / sHandle->wp.gasMass);
		//sqrt(29369)=171.3766= sqrt(8*R*1000/PI)*3PI/8, that is, the constant part of the v_avg=sqrt(8RT/PI/m/0.001)) found in literature, multiplied by
		//the corrective factor of 3PI/8 that accounts for moving from volumetric speed distribution to wall collision speed distribution
		sHandle->currentParticle.velocity = sqrt(oldSpeed2 + (newSpeed2 - oldSpeed2)*collidedFacet->sh.accomodationFactor);
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
		return rnd()*sHandle->wp.latestMoment; //continous desorption between 0 and latestMoment
	}
}

double GetStickingAt(SubprocessFacet *f, double time) {
	if (f->sh.sticking_paramId == -1) //constant sticking
		return f->sh.sticking;
	else return sHandle->parameters[f->sh.sticking_paramId].InterpolateY(time, false);
}

double GetOpacityAt(SubprocessFacet *f, double time) {
	if (f->sh.opacity_paramId == -1) //constant sticking
		return f->sh.opacity;
	else return sHandle->parameters[f->sh.opacity_paramId].InterpolateY(time, false);
}

void TreatMovingFacet() {
	Vector3d localVelocityToAdd;
	if (sHandle->wp.motionType == 1) {
		localVelocityToAdd = sHandle->wp.motionVector2;
	}
	else if (sHandle->wp.motionType == 2) {
		Vector3d distanceVector = 0.01*(sHandle->currentParticle.position - sHandle->wp.motionVector1); //distance from base, with cm->m conversion
		localVelocityToAdd = CrossProduct(sHandle->wp.motionVector2, distanceVector);
	}
	Vector3d oldVelocity, newVelocity;
	oldVelocity = sHandle->currentParticle.direction*sHandle->currentParticle.velocity;
	newVelocity = oldVelocity + localVelocityToAdd;
	sHandle->currentParticle.direction = newVelocity.Normalized();
	sHandle->currentParticle.velocity = newVelocity.Norme();
}

void IncreaseFacetCounter(SubprocessFacet *f, double time, size_t hit, size_t desorb, size_t absorb, double sum_1_per_v, double sum_v_ort) {
	size_t nbMoments = sHandle->moments.size();
	for (size_t m = 0; m <= nbMoments; m++) {
		if (m == 0 || abs(time - sHandle->moments[m - 1]) < sHandle->wp.timeWindowSize / 2.0) {
			f->tmpCounter[m].hit.nbMCHit += hit;
			double hitEquiv = static_cast<double>(hit)*sHandle->currentParticle.oriRatio;
			f->tmpCounter[m].hit.nbHitEquiv += hitEquiv;
			f->tmpCounter[m].hit.nbDesorbed += desorb;
			f->tmpCounter[m].hit.nbAbsEquiv += static_cast<double>(absorb)*sHandle->currentParticle.oriRatio;
			f->tmpCounter[m].hit.sum_1_per_ort_velocity += sHandle->currentParticle.oriRatio * sum_1_per_v;
			f->tmpCounter[m].hit.sum_v_ort += sHandle->currentParticle.oriRatio * sum_v_ort;
			f->tmpCounter[m].hit.sum_1_per_velocity += (hitEquiv + static_cast<double>(desorb)) / sHandle->currentParticle.velocity;
		}
	}
}

void SubprocessFacet::ResetCounter() {
	std::fill(tmpCounter.begin(), tmpCounter.end(), FacetHitBuffer());
}

void SubprocessFacet::ResizeCounter(size_t nbMoments) {
	tmpCounter = std::vector<FacetHitBuffer>(nbMoments + 1); //Includes 0-init
	tmpCounter.shrink_to_fit();
	tmpHistograms = std::vector<FacetHistogramBuffer>(nbMoments + 1);
	tmpHistograms.shrink_to_fit();
}

void SubprocessFacet::RegisterTransparentPass()
{
	double directionFactor = abs(Dot(sHandle->currentParticle.direction, this->sh.N));
	IncreaseFacetCounter(this, sHandle->currentParticle.flightTime + this->colDist / 100.0 / sHandle->currentParticle.velocity, 1, 0, 0, 2.0 / (sHandle->currentParticle.velocity*directionFactor), 2.0*(sHandle->wp.useMaxwellDistribution ? 1.0 : 1.1781)*sHandle->currentParticle.velocity*directionFactor);

	this->hitted = true;
	if (/*this->texture &&*/ this->sh.countTrans) {
		RecordHitOnTexture(this, sHandle->currentParticle.flightTime + this->colDist / 100.0 / sHandle->currentParticle.velocity,
			true, 2.0, 2.0);
	}
	if (/*this->direction &&*/ this->sh.countDirection) {
		RecordDirectionVector(this, sHandle->currentParticle.flightTime + this->colDist / 100.0 / sHandle->currentParticle.velocity);
	}
	LogHit(this);
	ProfileFacet(this, sHandle->currentParticle.flightTime + this->colDist / 100.0 / sHandle->currentParticle.velocity,
		true, 2.0, 2.0);
	if (this->sh.anglemapParams.record) RecordAngleMap(this);
}
