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
#pragma once

#include "MolflowTypes.h"
#include "Buffer_shared.h" //Facetproperties
#include "SMP.h"
#include <vector>
#include "Vector.h"
#include "Parameter.h"
#include <tuple>

class Anglemap {
public:
	std::vector<size_t>   pdf;		  // Incident angle distribution, phi and theta, not normalized. Used either for recording or for 2nd order interpolation
	std::vector<double>   phi_CDFs;    // A table containing phi distributions for each theta, starting from 0 for every line (1 line = 1 theta value). For speed we keep it in one memory block, 1 pointer
	std::vector<size_t>   phi_CDFsums; // since CDF runs only to the middle of the last segment, for each theta a line sum is stored here. Also a pdf for theta
	std::vector<double>   theta_CDF;	  // Theta CDF, not normalized. nth value is the CDF at the end of region n (beginning of first section is always 0)
	size_t   theta_CDFsum; // since theta CDF only runs till the middle of the last segment, the map sum is here

	double GetTheta(const double& thetaIndex, const AnglemapParams& anglemapParams);
	double GetPhi(const double& phiIndex, const AnglemapParams& anglemapParams);
	double GetPhipdfValue(const double & thetaIndex, const int & phiIndex, const AnglemapParams & anglemapParams);
	double GetPhiCDFValue(const double& thetaIndex, const int& phiIndex, const AnglemapParams& anglemapParams);
	double GetPhiCDFSum(const double & thetaIndex, const AnglemapParams & anglemapParams);
	std::tuple<double, int, double> GenerateThetaFromAngleMap(const AnglemapParams& anglemapParams);
	double GeneratePhiFromAngleMap(const int& thetaLowerIndex, const double& thetaOvershoot, const AnglemapParams& anglemapParams);
};

// Local facet structure
class SubprocessFacet {
public:
	FacetProperties sh;

	std::vector<size_t>      indices;          // Indices (Reference to geometry vertex)
	std::vector<Vector2d> vertices2;        // Vertices (2D plane space, UV coordinates)
	std::vector<std::vector<TextureCell>>     texture;            // Texture hit recording (taking area, temperature, mass into account), 1+nbMoments
	std::vector<double>   textureCellIncrements;              // Texure increment
	std::vector<bool>     largeEnough;      // cells that are NOT too small for autoscaling
	double   fullSizeInc;       // Texture increment of a full texture element
	std::vector<std::vector<DirectionCell>>     direction;       // Direction field recording (average), 1+nbMoments
	//bool     *fullElem;         // Direction field recording (only on full element)
	std::vector<std::vector<ProfileSlice>> profile;         // Distribution and hit recording
	std::vector<double>   outgassingMap; // Cumulative outgassing map when desorption is based on imported file
	double outgassingMapWidthD; //actual outgassing file map width
	double outgassingMapHeightD; //actual outgassing file map height
	Anglemap angleMap;

	// Temporary var (used in Intersect for collision)
	double colDist;
	double colU;
	double colV;
	double rw;
	double rh;
	double iw;
	double ih;

	// Temporary var (used in FillHit for hit recording)
	bool   hitted;
	bool   ready;         // Volatile state
	size_t    textureSize;   // Texture size (in bytes)
	size_t    profileSize;   // profile size (in bytes)
	size_t    directionSize; // direction field size (in bytes)
	size_t    angleMapSize;  // incidentangle map size (in bytes)

	/*int CDFid; //Which probability distribution it belongs to (one CDF per temperature)
	int IDid;  //If time-dependent desorption, which is its ID*/
	size_t globalId; //Global index (to identify when superstructures are present)

	// Facet hit counters
	std::vector<FacetHitBuffer> tmpCounter; //1+nbMoment
	std::vector<FacetHistogramBuffer> tmpHistograms; //1+nbMoment
	void  ResetCounter();
	void	ResizeCounter(size_t nbMoments);
	bool  InitializeOnLoad(const size_t& globalId);

	void InitializeHistogram();

	bool InitializeDirectionTexture();

	bool InitializeProfile();

	bool InitializeTexture();

	bool InitializeAngleMap();

	void InitializeOutgassingMap();

	bool InitializeLinkAndVolatile(const size_t & id);

	void RegisterTransparentPass(); //Allows one shared Intersect routine between MolFlow and Synrad
};

// Local simulation structure

class AABBNODE;

class SuperStructure {
public:
	SuperStructure();
	~SuperStructure();
	std::vector<SubprocessFacet>  facets;   // Facet handles
	AABBNODE* aabbTree; // Structure AABB tree
};

class CurrentParticleStatus {
public:
	Vector3d position;    // Position
	Vector3d direction;    // Direction
	double oriRatio;
	size_t   nbBounces; // Number of hit (current particle) since desorption
	double   distanceTraveled;
	double   velocity;
	double   flightTime;
	double   expectedDecayMoment; //for radioactive gases
	size_t   structureId;        // Current structure
	int      teleportedFrom;   // We memorize where the particle came from: we can teleport back
	SubprocessFacet *lastHitFacet;     // Last hitted facet
	std::vector<SubprocessFacet*> transparentHitBuffer; //Storing this buffer simulation-wide is cheaper than recreating it at every Intersect() call
};

class Simulation {
public:

	Simulation();
	GlobalHitBuffer tmpGlobalResult; //Global results since last UpdateMCHits
	std::vector<FacetHistogramBuffer> tmpGlobalHistograms; //Recorded histogram since last UpdateMCHits, 1+nbMoment copies
	std::vector<ParticleLoggerItem> tmpParticleLog; //Recorded particle log since last UpdateMCHits

	llong totalDesorbed;           // Total number of desorptions (for this process, not reset on UpdateMCHits)

	std::vector<std::vector<std::pair<double, double>>> CDFs; //cumulative distribution function for each temperature
	std::vector<std::vector<std::pair<double, double>>> IDs; //integrated distribution function for each time-dependent desorption type
	std::vector<double> temperatures; //keeping track of all temperatures that have a CDF already generated
	std::vector<double> moments;      //time values (seconds) when a simulation state is measured
	std::vector<size_t> desorptionParameterIDs; //time-dependent parameters which are used as desorptions, therefore need to be integrated
	std::vector<Parameter> parameters; //Time-dependent parameters 

	// Geometry
	GeomProperties sh;
	WorkerParams wp;
	OntheflySimulationParams ontheflyParams;

	std::vector<Vector3d>   vertices3;        // Vertices
	std::vector<SuperStructure> structures; //They contain the facets  

	double stepPerSec;  // Avg number of step per sec
	//Facet size counters
	size_t textTotalSize;  // Texture total size
	size_t profTotalSize;  // Profile total size
	size_t dirTotalSize;   // Direction field total size
	size_t angleMapTotalSize;
	size_t histogramTotalSize;
	bool loadOK;        // Load OK flag
	bool lastHitUpdateOK;  // Last hit update timeout
	bool lastLogUpdateOK; // Last log update timeout
	bool hasVolatile;   // Contains volatile facet
	double calcACTime;  // AC matrix calculation time

	// Particle coordinates (MC)
	CurrentParticleStatus currentParticle;


	// Angular coefficient (opaque facets)
	size_t     nbAC;
	ACFLOAT *acMatrix;
	ACFLOAT *acDensity;
	ACFLOAT *acDesorb;
	ACFLOAT *acAbsorb;
	ACFLOAT *acRho;
	ACFLOAT *acArea;
	double  *acLines;
	size_t     prgAC;

	// Angular coefficient (transparent facets)
	size_t     nbACT;
	ACFLOAT *acTMatrix;
	ACFLOAT *acTDensity;
	ACFLOAT *acTArea;
	double  *acTLines;

#ifdef JACOBI_ITERATION
	ACFLOAT *acDensityTmp;
#endif

};
// -- Methods ---------------------------------------------------

void RecordHitOnTexture(SubprocessFacet *f, double time, bool countHit, double velocity_factor, double ortSpeedFactor);
void RecordDirectionVector(SubprocessFacet *f, double time);
void ProfileFacet(SubprocessFacet *f, double time, bool countHit, double velocity_factor, double ortSpeedFactor);
void LogHit(SubprocessFacet *f);
void RecordAngleMap(SubprocessFacet* collidedFacet);
void InitSimulation();
void ClearSimulation();
void SetState(size_t state, const char *status, bool changeState = true, bool changeStatus = true);
void SetErrorSub(const char *msg);
void ClearACMatrix();
bool LoadSimulation(Dataport *loader);
bool UpdateOntheflySimuParams(Dataport *loader);
bool StartSimulation(size_t sMode);
void ResetSimulation();
bool SimulationRun();
bool SimulationMCStep(size_t nbStep);
void IncreaseDistanceCounters(double d);
bool SimulationACStep(int nbStep);
void RecordHit(const int& type);
void RecordLeakPos();
bool StartFromSource();
void PerformBounce(SubprocessFacet *iFacet);
void RecordAbsorb(SubprocessFacet *iFacet);
void RecordHistograms(SubprocessFacet * iFacet);
void PerformTeleport(SubprocessFacet *iFacet);
void PerformTransparentPass(SubprocessFacet *iFacet);
void UpdateHits(Dataport *dpHit, Dataport *dpLog, int prIdx, DWORD timeout);
void UpdateLog(Dataport *dpLog, DWORD timeout);
void UpdateMCHits(Dataport *dpHit, int prIdx, size_t nbMoments, DWORD timeout);
void UpdateACHits(Dataport *dpHit, int prIdx, DWORD timeout);
void ResetTmpCounters();

double GetTick();
size_t   GetHitsSize();
bool ComputeACMatrix(SHELEM_OLD *mesh);

int GetIDId(int paramId);

void   UpdateVelocity(SubprocessFacet *collidedFacet);
double GenerateRandomVelocity(int CDFId);
double GenerateDesorptionTime(SubprocessFacet* src);
double GetStickingAt(SubprocessFacet *src, double time);
double GetOpacityAt(SubprocessFacet *src, double time);
void   IncreaseFacetCounter(SubprocessFacet *f, double time, size_t hit, size_t desorb, size_t absorb, double sum_1_per_v, double sum_v_ort);
void   TreatMovingFacet();