#pragma once
#include "GLApp/GlTypes.h"
//#include "Buffer_shared.h"

// Desorption type
#define DES_NONE    0   // No desorption
#define DES_UNIFORM 1   // Uniform
#define DES_COSINE  2   // cos(theta)
#define DES_COSINE_N 3 // cos(theta)^N
#define DES_ANGLEMAP 4 //imported file

// (Old) Reflection types
#define REF_DIFFUSE 0   // Diffuse (cosine law)
#define REF_MIRROR  1   // Mirror
#define REF_UNIFORM 2   // Uniform (for testing)

// Profile type
#define REC_NONE       0  // No recording
#define REC_PRESSUREU  1  // Pressure and density profile (U direction)
#define REC_PRESSUREV  2  // Pressure and density profile (V direction)
#define REC_ANGULAR    3  // Angular profile
#define REC_VELOCITY   4 //Velocity distribution
#define REC_ORT_VELOCITY 5 //Orthogonal velocity component

// Hit type
#define HIT_DES   1
#define HIT_ABS   2
#define HIT_REF   3
#define HIT_TRANS 4
#define HIT_TELEPORT 5
#define HIT_MOVING 6
#define HIT_LAST 7

#define MC_MODE 0         // Monte Carlo simulation mode
#define AC_MODE 1         // Angular coefficient simulation mode

typedef float ACFLOAT;

// Density/Hit field stuff
#define HITMAX 1E38


typedef struct {
	llong count;
	double sum_v_ort;
	double sum_1_per_ort_velocity;
} APROFILE;

typedef struct {
	llong count;
	double sum_v_ort_per_area;
	double sum_1_per_ort_velocity;
} AHIT;

//Texture limit types
typedef struct {
	double all;
	double moments_only;
} TEXTURE_MOMENT_TYPE;

typedef struct {
	TEXTURE_MOMENT_TYPE min;
	TEXTURE_MOMENT_TYPE max;
} TEXTURE_MIN_MAX;

typedef struct {
	TEXTURE_MIN_MAX manual;
	TEXTURE_MIN_MAX autoscale;
} TEXTURE_SCALE_TYPE;

class AnglemapParams {
public:
	bool   record; // Record incident angle 2-dim distribution
	bool hasRecorded;
	size_t phiWidth; //resolution between -PI and +PI
	double thetaLimit; //angle map can have a different resolution under and over the limit. Must be between 0 and PI/2
	size_t thetaLowerRes; //resolution between 0 and angleMapThetaLimit
	size_t thetaHigherRes; //resolution between angleMapThetaLimit and PI/2
};

class Reflection {
public:
	double diffusePart;
	double specularPart;
	double cosineExponent; //Cos^N part: 1-diffuse-specular
};

//Just for AC matrix calculation in Molflow, old mesh structure:
typedef struct {

	float   area;     // Area of element
	float   uCenter;  // Center coordinates
	float   vCenter;  // Center coordinates
	int     elemId;   // Element index (MESH array)
	bool    full;     // Element is full

} SHELEM_OLD;

