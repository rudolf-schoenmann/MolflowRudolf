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
#ifdef WIN
#define NOMINMAX
//#include <windows.h> // For GetTickCount()
#include <Process.h> // For _getpid()
#else
//#include <time.h>
//#include <sys/time.h>
#endif

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "Simulation.h"
#include "IntersectAABB_shared.h"
#include "Random.h"
#include <sstream>

// Global handles

SubprocessFacet** THitCache;
SIMULATION *sHandle;

// Timing stuff

#ifdef WIN
bool usePerfCounter;         // Performance counter usage
LARGE_INTEGER perfTickStart; // First tick
double perfTicksPerSec;      // Performance counter (number of tick per second)
#endif
DWORD tickStart;

void InitSimulation() {

	// Global handle allocation
	sHandle = new SIMULATION();
	THitCache = new SubprocessFacet*[MAX_THIT]; // Transparent hit cache

#ifdef WIN
	{
		LARGE_INTEGER qwTicksPerSec;
		usePerfCounter = QueryPerformanceFrequency(&qwTicksPerSec);
		if (usePerfCounter) {
			QueryPerformanceCounter(&perfTickStart);
			perfTicksPerSec = (double)qwTicksPerSec.QuadPart;
		}
		tickStart = GetTickCount();
	}
#else
	tickStart = (DWORD)time(NULL);
#endif

}

void ClearSimulation() {

	int i, j;

	// Free old stuff
	sHandle->CDFs = std::vector<std::vector<std::pair<double, double>>>(); //clear CDF distributions

	SAFE_FREE(sHandle->vertices3);
	for (j = 0; j < sHandle->nbSuper; j++) {
		for (i = 0; i < sHandle->str[j].nbFacet; i++) {
			SubprocessFacet *f = sHandle->str[j].facets[i];
			if (f) {
				SAFE_FREE(f->indices);
				SAFE_FREE(f->vertices2);
				//SAFE_FREE(f->fullElem);
				SAFE_FREE(f->inc);
				SAFE_FREE(f->largeEnough);
				for (size_t m = 0; m < (sHandle->nbMoments + 1); m++) {
					if (f->texture) {
						SAFE_FREE(f->texture[m]);
					}
					if (f->profile) {
						SAFE_FREE(f->profile[m]);
					}
					if (f->direction) {
						SAFE_FREE(f->direction[m]);
					}
					//if (f->velocityHistogram) SAFE_FREE(f->velocityHistogram);
				}
				SAFE_FREE(f->outgassingMap);
				SAFE_FREE(f->angleMap.pdf);
				SAFE_FREE(f->angleMap.phi_CDFs);
				SAFE_FREE(f->angleMap.phi_CDFsums);
				SAFE_FREE(f->angleMap.theta_CDF);
				SAFE_FREE(f->texture);
				SAFE_FREE(f->profile);
				SAFE_FREE(f->direction);
				//SAFE_FREE(f->velocityHistogram);
				f->counter.clear();
				f->counter.shrink_to_fit();
				delete(f); f = NULL;
			}

		}
		SAFE_FREE(sHandle->str[j].facets);
		if (sHandle->str[j].aabbTree) {
			DestroyAABB(sHandle->str[j].aabbTree->left);
			DestroyAABB(sHandle->str[j].aabbTree->right);
			free(sHandle->str[j].aabbTree);
			sHandle->str[j].aabbTree = NULL;
		}
	}

	sHandle->CDFs.clear(); sHandle->CDFs.shrink_to_fit();
	sHandle->IDs.clear();sHandle->IDs.shrink_to_fit();
	sHandle->parameters.clear();sHandle->parameters.shrink_to_fit();
	sHandle->temperatures.clear();
	sHandle->moments.clear();
	sHandle->desorptionParameterIDs.clear();
	sHandle->tmpParticleLog.clear();sHandle->tmpParticleLog.shrink_to_fit();
	

	ClearACMatrix();
	memset(sHandle, 0, sizeof(SIMULATION)); //Will it corrupt vectors?

}

DWORD RevertBit(DWORD dw) {
	DWORD dwIn = dw;
	DWORD dwOut = 0;
	DWORD dWMask = 1;
	int i;
	for (i = 0; i < 32; i++) {
		if (dwIn & 0x80000000UL) dwOut |= dWMask;
		dwIn = dwIn << 1;
		dWMask = dWMask << 1;
	}
	return dwOut;
}

DWORD GetSeed() {

	/*#ifdef WIN
	DWORD r;
	_asm {
	rdtsc
	mov r,eax
	}
	return RevertBit(r ^ (DWORD)(_getpid()*65519));
	#else*/
	int processId;
#ifdef  WIN
	processId = _getpid();
#else
	processId = ::getpid();
#endif //  WIN


	return (DWORD)((int)(GetTick()*1000.0)*_getpid());
	//#endif

}

/*double Norme(Vector3d *v) { //already defined
	return sqrt(v->x*v->x + v->y*v->y + v->z*v->z);
}*/

bool LoadSimulation(Dataport *loader) {

	size_t i, j, idx;
	BYTE *buffer;
	BYTE *globalBuff;
	BYTE *bufferStart;
	GeomProperties *shGeom;
	Vector3d *shVert;
	double t1, t0;
	DWORD seed;
	char err[128];

	t0 = GetTick();

	sHandle->loadOK = false;

	SetState(PROCESS_STARTING, "Clearing previous simulation");
	ClearSimulation();
	
	/* //Mutex not necessary: by the time the COMMAND_LOAD is issued the interface releases the handle, concurrent reading is safe and it's only destroyed by the interface when all processes are ready loading
	   //Result: faster, parallel loading
	// Connect the dataport
	SetState(PROCESS_STARTING, "Waiting for 'loader' dataport access...");
	if (!AccessDataportTimed(loader,15000)) {
		SetErrorSub("Failed to connect to DP");
		return false;
	}
	*/

	SetState(PROCESS_STARTING, "Loading simulation");

	buffer = (BYTE *)loader->buff;
	bufferStart = buffer; //memorize start for later

	// Load new geom from the dataport

	shGeom = (GeomProperties *)buffer;
	if (shGeom->nbSuper > MAX_STRUCT) {
		//ReleaseDataport(loader);
		SetErrorSub("Too many structures");
		return false;
	}
	if (shGeom->nbSuper <= 0) {
		//ReleaseDataport(loader);
		SetErrorSub("No structures");

		return false;
	}

	sHandle->nbVertex = shGeom->nbVertex;
	sHandle->nbSuper = shGeom->nbSuper;
	sHandle->totalFacet = shGeom->nbFacet;
	sHandle->nbMoments = shGeom->nbMoments;
	sHandle->latestMoment = shGeom->latestMoment;
	sHandle->totalDesorbedMolecules = shGeom->totalDesorbedMolecules;
	sHandle->finalOutgassingRate = shGeom->finalOutgassingRate;
	sHandle->gasMass = shGeom->gasMass;
	sHandle->enableDecay = shGeom->enableDecay;
	sHandle->halfLife = shGeom->halfLife;
	sHandle->timeWindowSize = shGeom->timeWindowSize;
	sHandle->useMaxwellDistribution = shGeom->useMaxwellDistribution;
	sHandle->calcConstantFlow = shGeom->calcConstantFlow;
	sHandle->motionType = shGeom->motionType;
	sHandle->motionVector1 = shGeom->motionVector1;
	sHandle->motionVector2 = shGeom->motionVector2;
	// Prepare super structure (allocate memory for facets)
	buffer += sizeof(GeomProperties);
	sHandle->ontheflyParams = READBUFFER(OntheflySimulationParams);
	if (sHandle->ontheflyParams.enableLogging) sHandle->tmpParticleLog.reserve(sHandle->ontheflyParams.logLimit / sHandle->ontheflyParams.nbProcess);
	buffer += sizeof(Vector3d)*sHandle->nbVertex;
	for (i = 0; i < sHandle->totalFacet; i++) {
		FacetProperties *shFacet = (FacetProperties *)buffer;
		sHandle->str[shFacet->superIdx].nbFacet++;
		buffer += sizeof(FacetProperties) + shFacet->nbIndex*(sizeof(size_t) + sizeof(Vector2d));
		if (shFacet->useOutgassingFile) buffer += sizeof(double)*shFacet->outgassingMapWidth*shFacet->outgassingMapHeight;
		if (shFacet->anglemapParams.hasRecorded) buffer += sizeof(size_t)*shFacet->anglemapParams.phiWidth*(shFacet->anglemapParams.thetaLowerRes+shFacet->anglemapParams.thetaHigherRes);
		if (shFacet->isTextured) buffer += sizeof(double)*shFacet->texWidth*shFacet->texHeight;
	}
	for (i = 0; i < sHandle->nbSuper; i++) {
		int nbF = sHandle->str[i].nbFacet;
		if (nbF == 0) {
			/*ReleaseDataport(loader);
			sprintf(err,"Structure #%d has no facets!",i+1);
			SetErrorSub(err);
			return false;*/
		}
		else {

			sHandle->str[i].facets = (SubprocessFacet **)malloc(nbF * sizeof(SubprocessFacet *));
			memset(sHandle->str[i].facets, 0, nbF * sizeof(SubprocessFacet *));
			//sHandle->str[i].nbFacet = 0;
		}
		sHandle->str[i].nbFacet = 0;
	}

	globalBuff = buffer; //after facets and outgassing maps, with inc values

	//buffer = (BYTE *)loader->buff;
	buffer = bufferStart; //start from beginning again

	// Name
	memcpy(sHandle->name, shGeom->name, 64);

	// Vertices

	sHandle->vertices3 = (Vector3d *)malloc(sHandle->nbVertex * sizeof(Vector3d));
	if (!sHandle->vertices3) {
		SetErrorSub("Not enough memory to load vertices");
		return false;
	}
	buffer += sizeof(GeomProperties);
	buffer += sizeof(OntheflySimulationParams);
	shVert = (Vector3d *)(buffer);
	memcpy(sHandle->vertices3, shVert, sHandle->nbVertex * sizeof(Vector3d));
	buffer += sizeof(Vector3d)*sHandle->nbVertex;

	// Facets
	for (i = 0; i < sHandle->totalFacet; i++) {

		FacetProperties *shFacet = (FacetProperties *)buffer;
		SubprocessFacet *f = new SubprocessFacet;
		if (!f) {
			SetErrorSub("Not enough memory to load facets");
			return false;
		}
		memset(f, 0, sizeof(SubprocessFacet));
		memcpy(&(f->sh), shFacet, sizeof(FacetProperties));
		f->ResizeCounter(sHandle->nbMoments); //Initialize counter

		//f->sh.maxSpeed = 4.0 * sqrt(2.0*8.31*f->sh.temperature/0.001/sHandle->gasMass);

		sHandle->hasVolatile |= f->sh.isVolatile;
		sHandle->hasDirection |= f->sh.countDirection;

		idx = f->sh.superIdx;
		sHandle->str[idx].facets[sHandle->str[idx].nbFacet] = f;
		sHandle->str[idx].facets[sHandle->str[idx].nbFacet]->globalId = i;
		sHandle->str[idx].nbFacet++;

		if (f->sh.superDest || f->sh.isVolatile) {
			// Link or volatile facet, overides facet settings
			// Must be full opaque and 0 sticking
			// (see SimulationMC.c::PerformBounce)
			//f->sh.isOpaque = true;
			f->sh.opacity = 1.0;
			f->sh.opacity_paramId = -1;
			f->sh.sticking = 0.0;
			f->sh.sticking_paramId = -1;
			if (((f->sh.superDest - 1) >= sHandle->nbSuper || f->sh.superDest < 0)) {
				// Geometry error
				ClearSimulation();
				//ReleaseDataport(loader);
				sprintf(err, "Invalid structure (wrong link on F#%zd)", i + 1);
				SetErrorSub(err);
				return false;
			}
		}

		// Reset counter in local memory
		//memset(&(f->sh.counter), 0, sizeof(FacetHitBuffer));
		f->indices = (size_t *)malloc(f->sh.nbIndex * sizeof(size_t));
		buffer += sizeof(FacetProperties);
		memcpy(f->indices, buffer, f->sh.nbIndex * sizeof(size_t));
		buffer += f->sh.nbIndex * sizeof(size_t);
		f->vertices2 = (Vector2d *)malloc(f->sh.nbIndex * sizeof(Vector2d));
		if (!f->vertices2) {
			SetErrorSub("Not enough memory to load vertices");
			return false;
		}
		memcpy(f->vertices2, buffer, f->sh.nbIndex * sizeof(Vector2d));
		buffer += f->sh.nbIndex * sizeof(Vector2d);
		//Outgassing map
		if (f->sh.useOutgassingFile) {
			size_t nbE = f->sh.outgassingMapWidth*f->sh.outgassingMapHeight;
			f->outgassingMap = (double*)malloc(sizeof(double)*nbE);
			if (!f->outgassingMap) {
				SetErrorSub("Not enough memory to load outgassing map");
				return false;
			}
			memcpy(f->outgassingMap, buffer, sizeof(double)*nbE);
			buffer += sizeof(double)*nbE;
			//Precalc actual outgassing map width and height for faster generation:
			f->outgassingMapWidthD = f->sh.U.Norme() * f->sh.outgassingFileRatio;
			f->outgassingMapHeightD = f->sh.V.Norme() * f->sh.outgassingFileRatio;
			
			for (size_t i = 1; i < nbE; i++) {
				f->outgassingMap[i] += f->outgassingMap[i - 1]; //Convert p.d.f. to cumulative distr. f.
			}
		}
		//Incident angle map
		if (f->sh.anglemapParams.hasRecorded) {
			//Record or use to generate

			f->angleMapSize = f->sh.anglemapParams.phiWidth * (f->sh.anglemapParams.thetaLowerRes + f->sh.anglemapParams.thetaHigherRes) * sizeof(size_t);
			f->angleMap.pdf = (size_t*)malloc(f->angleMapSize);
			if (!f->angleMap.pdf) {
				SetErrorSub("Not enough memory to load incident angle map (values)");
				return false;
			}
			//Copy values if "USE" mode, copy 0 values if "RECORD" mode
			memcpy(f->angleMap.pdf, buffer, f->angleMapSize);
			buffer += f->angleMapSize;

			if (f->sh.desorbType == DES_ANGLEMAP) {
				//Construct CDFs				
				f->angleMap.phi_CDFsums = (size_t*)malloc(sizeof(size_t) * (f->sh.anglemapParams.thetaLowerRes + f->sh.anglemapParams.thetaHigherRes));
				if (!f->angleMap.phi_CDFsums) {
					SetErrorSub("Not enough memory to load incident angle map (phi CDF line sums)");
					return false;
				}
				f->angleMap.theta_CDF = (double*)malloc(sizeof(double) * (f->sh.anglemapParams.thetaLowerRes + f->sh.anglemapParams.thetaHigherRes));
				if (!f->angleMap.theta_CDF) {
					SetErrorSub("Not enough memory to load incident angle map (line sums, CDF)");
					return false;
				}
				f->angleMap.phi_CDFs = (double*)malloc(sizeof(double) * f->sh.anglemapParams.phiWidth * (f->sh.anglemapParams.thetaLowerRes + f->sh.anglemapParams.thetaHigherRes));
				if (!f->angleMap.phi_CDFs) {
					SetErrorSub("Not enough memory to load incident angle map (CDF)");
					return false;
				}
				
				//First pass: determine sums
				f->angleMap.theta_CDFsum = 0;
				memset(f->angleMap.phi_CDFsums, 0, sizeof(size_t) * (f->sh.anglemapParams.thetaLowerRes + f->sh.anglemapParams.thetaHigherRes));
				for (size_t thetaIndex = 0; thetaIndex < (f->sh.anglemapParams.thetaLowerRes + f->sh.anglemapParams.thetaHigherRes); thetaIndex++) {
					for (size_t phiIndex = 0; phiIndex < f->sh.anglemapParams.phiWidth; phiIndex++) {
						f->angleMap.phi_CDFsums[thetaIndex] +=	f->angleMap.pdf[thetaIndex*f->sh.anglemapParams.phiWidth + phiIndex];	
					}
					f->angleMap.theta_CDFsum += f->angleMap.phi_CDFsums[thetaIndex];
				}
				if (!f->angleMap.theta_CDFsum) {
					std::stringstream err; err << "Facet " << f->globalId + 1 << " has all-zero recorded angle map.";
					SetErrorSub(err.str().c_str());
					return false;
				}

				//Second pass: write CDFs
				double thetaNormalizingFactor = 1.0 / (double)f->angleMap.theta_CDFsum;
				for (size_t thetaIndex = 0; thetaIndex < (f->sh.anglemapParams.thetaLowerRes + f->sh.anglemapParams.thetaHigherRes); thetaIndex++) {	
					if (f->angleMap.theta_CDFsum == 0) { //no hits in this line, generate CDF of uniform distr.
						f->angleMap.theta_CDF[thetaIndex] = (0.5 + (double)thetaIndex) / (double)(f->sh.anglemapParams.thetaLowerRes + f->sh.anglemapParams.thetaHigherRes);
					}
					else {
						if (thetaIndex == 0) {
							//First CDF value, covers half of first segment
							f->angleMap.theta_CDF[thetaIndex] = 0.5 * (double)f->angleMap.phi_CDFsums[0] * thetaNormalizingFactor;
						}
						else {
							//value covering second half of last segment and first of current segment
							f->angleMap.theta_CDF[thetaIndex] = f->angleMap.theta_CDF[thetaIndex - 1] + (double)(f->angleMap.phi_CDFsums[thetaIndex - 1] + f->angleMap.phi_CDFsums[thetaIndex])*0.5*thetaNormalizingFactor;
						}
					}
					double phiNormalizingFactor = 1.0 / (double)f->angleMap.phi_CDFsums[thetaIndex];
					for (size_t phiIndex = 0; phiIndex < f->sh.anglemapParams.phiWidth; phiIndex++) {
						size_t index = f->sh.anglemapParams.phiWidth * thetaIndex + phiIndex;
						if (f->angleMap.phi_CDFsums[thetaIndex] == 0) { //no hits in this line, create CDF of uniform distr.
								f->angleMap.phi_CDFs[index] = (0.5 + (double)phiIndex) / (double)f->sh.anglemapParams.phiWidth;
							}
						else {
							if (phiIndex == 0) {
								//First CDF value, covers half of first segment
								f->angleMap.phi_CDFs[index] = 0.5 * (double)f->angleMap.pdf[f->sh.anglemapParams.phiWidth * thetaIndex] * phiNormalizingFactor;
							}
							else {
								//value covering second half of last segment and first of current segment
								f->angleMap.phi_CDFs[index] = f->angleMap.phi_CDFs[f->sh.anglemapParams.phiWidth * thetaIndex + phiIndex - 1] + (double)(f->angleMap.pdf[f->sh.anglemapParams.phiWidth * thetaIndex + phiIndex - 1] + f->angleMap.pdf[f->sh.anglemapParams.phiWidth * thetaIndex + phiIndex])*0.5*phiNormalizingFactor;
							}
						}
					}
				}

			}
			else {
				//Record
				f->angleMap.phi_CDFs = NULL;
				f->angleMap.theta_CDF = NULL;
				f->angleMap.phi_CDFsums = NULL;
			}
		}
		else {
			f->angleMapSize = 0;
			f->angleMap.pdf = NULL;
		}

		//Textures
		if (f->sh.isTextured) {
			size_t nbE = f->sh.texWidth*f->sh.texHeight;
			f->textureSize = nbE * sizeof(TextureCell);

			if ((f->texture = (TextureCell **)malloc(sizeof(TextureCell *)* (1 + sHandle->nbMoments))) == NULL) {
				//ReleaseDataport(loader);
				SetErrorSub("Couldn't allocate memory (time moments container, textures)");
				return false;
			}
			memset(f->texture, 0, sizeof(TextureCell *)* (1 + sHandle->nbMoments)); //set all pointers to NULL so we can use SAFE_FREE later
			for (size_t m = 0; m < (1 + sHandle->nbMoments); m++) {
				if ((f->texture[m] = (TextureCell *)malloc(f->textureSize)) == NULL) { //steady-state plus one for each moment
					//ReleaseDataport(loader);
					SetErrorSub("Couldn't allocate memory (textures)");
					return false;
				}
				memset(f->texture[m], 0, f->textureSize);
			}

			//Load inc values (1/area)
				f->inc = (double *)malloc(nbE * sizeof(double));
				f->largeEnough = (bool *)malloc(sizeof(bool)*nbE);
				//f->fullElem = (bool *)malloc(sizeof(bool)*nbE);
				if (!(f->inc && f->largeEnough /*&& f->fullElem*/)) {
					SetErrorSub("Not enough memory to load");
					return false;
				}
				f->fullSizeInc = 1E30;
				for (j = 0; j < nbE; j++) {
					double incVal = READBUFFER(double);
					f->inc[j] = incVal;
					if ((f->inc[j] > 0.0) && (f->inc[j] < f->fullSizeInc)) f->fullSizeInc = f->inc[j];
				}
				for (j = 0; j < nbE; j++) { //second pass, filter out very small cells
					f->largeEnough[j] = (f->inc[j] < ((5.0f)*f->fullSizeInc));
				}
				sHandle->textTotalSize += f->textureSize*(1 + sHandle->nbMoments);

				f->iw = 1.0 / (double)f->sh.texWidthD;
				f->ih = 1.0 / (double)f->sh.texHeightD;
				f->rw = f->sh.U.Norme() * f->iw;
				f->rh = f->sh.V.Norme() * f->ih;
		}
		else f->textureSize = 0;

		//Profiles
		if (f->sh.isProfile) {
			f->profileSize = PROFILE_SIZE * sizeof(ProfileSlice);
			/*f->profile = (llong *)malloc(f->profileSize);
			memset(f->profile,0,f->profileSize);*/
			if ((f->profile = (ProfileSlice **)malloc(sizeof(ProfileSlice *)* (1 + sHandle->nbMoments))) == NULL) {
				//ReleaseDataport(loader);
				SetErrorSub("Couldn't allocate memory (time moments container, profiles)");
				return false;
			}
			memset(f->profile, 0, sizeof(ProfileSlice *)* (1 + sHandle->nbMoments)); //set all pointers to NULL so we can use SAFE_FREE later
			for (size_t m = 0; m < (1 + sHandle->nbMoments); m++) {
				if ((f->profile[m] = (ProfileSlice *)malloc(f->profileSize)) == NULL) { //steady-state plus one for each moment
					//ReleaseDataport(loader);
					SetErrorSub("Couldn't allocate memory (profiles)");
					return false;
				}
				memset(f->profile[m], 0, f->profileSize);
			}
			sHandle->profTotalSize += f->profileSize*(1 + sHandle->nbMoments);
		}
		else f->profileSize = 0;

		//Direction
		if (f->sh.countDirection) {
			f->directionSize = f->sh.texWidth*f->sh.texHeight * sizeof(DirectionCell);
			/*f->direction = (DirectionCell *)malloc(f->directionSize);
			memset(f->direction,0,f->directionSize);*/
			if ((f->direction = (DirectionCell **)malloc(sizeof(DirectionCell *)* (1 + sHandle->nbMoments))) == NULL) {
				//ReleaseDataport(loader);
				SetErrorSub("Couldn't allocate memory (time moments container, direction vectors)");
				return false;
			}
			memset(f->direction, 0, sizeof(DirectionCell *)* (1 + sHandle->nbMoments)); //set all pointers to NULL so we can use SAFE_FREE later
			for (size_t m = 0; m < (1 + sHandle->nbMoments); m++) {
				if ((f->direction[m] = (DirectionCell *)malloc(f->directionSize)) == NULL) { //steady-state plus one for each moment
					//ReleaseDataport(loader);
					SetErrorSub("Couldn't allocate memory (direction vectors)");
					return false;
				}
				memset(f->direction[m], 0, f->directionSize);
			}
			sHandle->dirTotalSize += f->directionSize*(1 + sHandle->nbMoments);
		}
		else f->directionSize = 0;

	}

	
	buffer = globalBuff;

	//CDFs
	size_t size1 = READBUFFER(size_t);
	sHandle->CDFs.reserve(size1);
	for (size_t i = 0; i < size1; i++) {
		std::vector<std::pair<double, double>> newCDF;
		size_t size2 = READBUFFER(size_t);
		newCDF.reserve(size2);
		for (size_t j = 0; j < size2; j++) {
			double valueX = READBUFFER(double);
			double valueY = READBUFFER(double);
			newCDF.push_back(std::make_pair(valueX, valueY));
		}
		sHandle->CDFs.push_back(newCDF);
	}

	//IDs
	size1 = READBUFFER(size_t);
	sHandle->IDs.reserve(size1);
	for (size_t i = 0; i < size1; i++) {
		std::vector<std::pair<double, double>> newID;
		size_t size2 = READBUFFER(size_t);
		newID.reserve(size2);
		for (size_t j = 0; j < size2; j++) {
			double valueX = READBUFFER(double);
			double valueY = READBUFFER(double);
			newID.push_back(std::make_pair(valueX, valueY));
		}
		sHandle->IDs.push_back(newID);
	}

	//Parameters
	size1 = READBUFFER(size_t);
	sHandle->parameters.reserve(size1);
	for (size_t i = 0; i < size1; i++) {
		Parameter newParam = Parameter();
		std::vector<std::pair<double, double>> newValues;
		size_t size2 = READBUFFER(size_t);
		newValues.reserve(size2);
		for (size_t j = 0; j < size2; j++) {
			double valueX = READBUFFER(double);
			double valueY = READBUFFER(double);
			newValues.push_back(std::make_pair(valueX, valueY));
		}
		newParam.SetValues(newValues, false);
		sHandle->parameters.push_back(newParam);
	}

	//Temperatures
	size1 = READBUFFER(size_t);
	sHandle->temperatures.reserve(size1);
	for (size_t i = 0; i < size1; i++) {
		double valueX = READBUFFER(double);
		sHandle->temperatures.push_back(valueX);
	}

	//Time moments
	//size1=READBUFFER(size_t);
	sHandle->moments.reserve(sHandle->nbMoments); //nbMoments already passed
	for (size_t i = 0; i < sHandle->nbMoments; i++) {
		double valueX = READBUFFER(double);
		sHandle->moments.push_back(valueX);
	}

	//Desorption parameter IDs
	size1 = READBUFFER(size_t);
	sHandle->desorptionParameterIDs.reserve(size1);
	for (size_t i = 0; i < size1; i++) {
		size_t valueX = READBUFFER(size_t);
		sHandle->desorptionParameterIDs.push_back(valueX);
	}

	/*//DEBUG

		Parameter *param1=new Parameter();
		param1->AddValue(0.5,0);
		param1->AddValue(0.5001,0.001);
		param1->AddValue(0.6,0.001);
		param1->AddValue(0.8,0);
		sHandle->parameters.push_back(*param1);

		sHandle->str[0].facets[0]->sh.sticking_paramId=0; //set facet 1 as time-dependent desorption*/
	
	//ReleaseDataport(loader); //Commented out as AccessDataport removed

	// Build all AABBTrees
	size_t maxDepth=0;
	for (i = 0; i < sHandle->nbSuper; i++)
		sHandle->str[i].aabbTree = BuildAABBTree(sHandle->str[i].facets, sHandle->str[i].nbFacet, 0,maxDepth);

	// Initialise simulation

	seed = GetSeed();
	rseed(seed);
	sHandle->loadOK = true;
	t1 = GetTick();
	printf("  Load %s successful\n", sHandle->name);
	printf("  Geometry: %zd vertex %zd facets\n", sHandle->nbVertex, sHandle->totalFacet);

	printf("  Geom size: %zd bytes\n", (size_t)(buffer - bufferStart));
	printf("  Number of stucture: %zd\n", sHandle->nbSuper);
	printf("  Global Hit: %zd bytes\n", sizeof(GlobalHitBuffer));
	printf("  Facet Hit : %zd bytes\n", sHandle->totalFacet * sizeof(FacetHitBuffer));
	printf("  Texture   : %zd bytes\n", sHandle->textTotalSize);
	printf("  Profile   : %zd bytes\n", sHandle->profTotalSize);
	printf("  Direction : %zd bytes\n", sHandle->dirTotalSize);

	printf("  Total     : %zd bytes\n", GetHitsSize());
	printf("  Seed: %lu\n", seed);
	printf("  Loading time: %.3f ms\n", (t1 - t0)*1000.0);
	return true;

}

bool UpdateOntheflySimuParams(Dataport *loader) {
	// Connect the dataport
	

	if (!AccessDataportTimed(loader, 2000)) {
		SetErrorSub("Failed to connect to loader DP");
		return false;
	}
	BYTE* buffer = (BYTE *)loader->buff;

	sHandle->ontheflyParams = READBUFFER(OntheflySimulationParams);
	ReleaseDataport(loader);

	return true;
}

void UpdateHits(Dataport *dpHit, Dataport* dpLog,int prIdx, DWORD timeout) {
	switch (sHandle->sMode) {
	case MC_MODE:
	{
		UpdateMCHits(dpHit, prIdx, sHandle->nbMoments, timeout);
		if (dpLog) UpdateLog(dpLog, timeout);
	}
		break;
	case AC_MODE:

		UpdateACHits(dpHit, prIdx, timeout);
		break;
	}

}

size_t GetHitsSize() {
	return sHandle->textTotalSize + sHandle->profTotalSize + sHandle->dirTotalSize + sHandle->totalFacet * sizeof(FacetHitBuffer) + sizeof(GlobalHitBuffer);
}

void ResetTmpCounters() {
	SetState(NULL, "Resetting local cache...", false, true);

	memset(&sHandle->tmpGlobalCount, 0, sizeof(FacetHitBuffer));

	sHandle->distTraveledSinceUpdate_total = 0.0;
	sHandle->distTraveledSinceUpdate_fullHitsOnly = 0.0;
	sHandle->nbLeakSinceUpdate = 0;
	sHandle->hitCacheSize = 0;
	sHandle->leakCacheSize = 0;

	for (int j = 0; j < sHandle->nbSuper; j++) {
		for (int i = 0; i < sHandle->str[j].nbFacet; i++) {
			SubprocessFacet *f = sHandle->str[j].facets[i];
			f->ResetCounter();
			f->hitted = false;

			if (f->texture) {
				for (size_t m = 0; m < (sHandle->nbMoments + 1); m++) {
					memset(f->texture[m], 0, f->textureSize);
				}
			}

			
			if (f->profile) {
				for (size_t m = 0; m < (sHandle->nbMoments + 1); m++) {
					memset(f->profile[m], 0, f->profileSize);
				}
			}

			
			if (f->direction) {
				for (size_t m = 0; m < (sHandle->nbMoments + 1); m++) {
					memset(f->direction[m], 0, f->directionSize);
				}
			}

			if (f->sh.anglemapParams.record) {
				memset(f->angleMap.pdf, 0, f->angleMapSize);
			}
		}
	}

}

void ResetSimulation() {
	sHandle->lastHitFacet = NULL;
	sHandle->totalDesorbed = 0;
	ResetTmpCounters();
	sHandle->tmpParticleLog.clear();
	if (sHandle->acDensity) memset(sHandle->acDensity, 0, sHandle->nbAC * sizeof(ACFLOAT));

}

bool StartSimulation(size_t sMode) {
	sHandle->sMode = sMode;
	switch (sMode) {
	case MC_MODE:
		if (!sHandle->lastHitFacet) StartFromSource();
		return (sHandle->lastHitFacet != NULL);
	case AC_MODE:
		if (sHandle->prgAC != 100) {
			SetErrorSub("AC matrix not calculated");
			return false;
		}
		else {
			sHandle->stepPerSec = 0.0;
			return true;
		}
	}

	SetErrorSub("Unknown simulation mode");
	return false;
}

void RecordHit(const int &type) {
	if (sHandle->hitCacheSize < HITCACHESIZE) {
		sHandle->hitCache[sHandle->hitCacheSize].pos = sHandle->pPos;
		sHandle->hitCache[sHandle->hitCacheSize].type = type;
		sHandle->hitCacheSize++;
	}
}

void RecordLeakPos() {
	// Source region check performed when calling this routine 
	// Record leak for debugging
	RecordHit(HIT_REF);
	RecordHit(HIT_LAST);
	if (sHandle->leakCacheSize < LEAKCACHESIZE) {
		sHandle->leakCache[sHandle->leakCacheSize].pos = sHandle->pPos;
		sHandle->leakCache[sHandle->leakCacheSize].dir = sHandle->pDir;
		sHandle->leakCacheSize++;
	}
}

bool SimulationRun() {

	// 1s step
	double t0, t1;
	int    nbStep = 1;
	bool   goOn;

	if (sHandle->stepPerSec == 0.0) {
		switch (sHandle->sMode) {
		case MC_MODE:

			nbStep = 250;
			break;
		case AC_MODE:
			nbStep = 1;
			break;
		}

	}

	if (sHandle->stepPerSec != 0.0)
		nbStep = (int)(sHandle->stepPerSec + 0.5);
	if (nbStep < 1) nbStep = 1;
	t0 = GetTick();
	switch (sHandle->sMode) {
	case MC_MODE:

		goOn = SimulationMCStep(nbStep);
		break;
	case AC_MODE:
		goOn = SimulationACStep(nbStep);
		break;
	}

	t1 = GetTick();
	sHandle->stepPerSec = (double)(nbStep) / (t1 - t0);
#ifdef _DEBUG
	printf("Running: stepPerSec = %f\n", sHandle->stepPerSec);
#endif

	return !goOn;

}

double GetTick() {

	// Number of sec since the application startup

#ifdef WIN

	if (usePerfCounter) {

		LARGE_INTEGER t, dt;
		QueryPerformanceCounter(&t);
		dt.QuadPart = t.QuadPart - perfTickStart.QuadPart;
		return (double)(dt.QuadPart) / perfTicksPerSec;

	}
	else {

		return (double)((GetTickCount() - tickStart) / 1000.0);

	}

#else

	if (tickStart < 0)
		tickStart = time(NULL);

	struct timeval tv;
	gettimeofday(&tv, NULL);
	return ((double)(tv.tv_sec - tickStart)*1000.0 + (double)tv.tv_usec / 1000.0);

#endif

	}

int GetIDId(int paramId) {

	int i;
	for (i = 0; i < (int)sHandle->desorptionParameterIDs.size() && (paramId != sHandle->desorptionParameterIDs[i]); i++); //check if we already had this parameter Id
	if (i >= (int)sHandle->desorptionParameterIDs.size()) i = -1; //not found
	return i;
}

