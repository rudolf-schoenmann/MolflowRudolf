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
#pragma once

#include <crtdbg.h> //To debug heap corruptions in memory

#include "Interface.h"
class Worker;
class ImportDesorption;
class TimeSettings;
class Movement;
class FacetAdvParams;
class FacetDetails;
class Viewer3DSettings;
class TextureScaling;
class GlobalSettings;
class HistoryWin;
class ProfilePlotter;
class PressureEvolution;
class TimewisePlotter;
class TexturePlotter;
class OutgassingMap;
class MomentsEditor;
class ParameterEditor;

class Error;

//const double carbondiameter = 2 * 76E-12;//deprecated since we now consider H2O as contaminant
const double carbondiameter = 276E-12; // 96 pm O-H bond length + 180 pm hydrogen brigde bond between H and next O
//This is not very accurated, espacially considering a quadratic lattice of H2O molecules on the surface.
//Improvement necessary!
const double kb = 1.38E-23;
const double tau = 1E-13;

class MolFlow : public Interface
{
public:
    MolFlow();
	
	//Public textfields so we can disable them from "Advanced facet parameters":
	GLTextField   *facetFlow;
	GLTextField   *facetFlowArea;
	
    
    void LoadFile(char *fName=NULL);
	void InsertGeometry(bool newStr,char *fName=NULL);
	void InsertGeometryBuffer(char *shortName);
	void SaveFile();
    
	//void ImportDesorption_DES(); //Deprecated
	void ExportProfiles();
	void ExportAngleMaps();
	void ImportAngleMaps();
	void CopyAngleMapToClipboard();
	void ClearAngleMapsOnSelection();
    void ClearFacetParams();
	
    void ApplyFacetParams();
	void UpdateFacetParams(bool updateSelection = false);
    void StartStopSimulation();
    void SaveConfig();
    void LoadConfig();
    
	void PlaceComponents();
    void UpdateFacetHits(bool allRows=false);
	void QuickPipe();
	float GetAppTime();
	void ClearParameters();
	void UpdatePlotters();

	//Flow/sticking coeff. conversion
	void calcFlow();
	void calcSticking();
	void calcStickingnew();
	void calcCoverage();
	void calcCovering();
	double calcDesorption();
	void calcDesorptionRate();
	double calcNmono();
	double calcdNsurf();
	
	bool EvaluateVariable(VLIST *v);

	//char* appTitle;

	GLButton      *texturePlotterShortcut;
	GLButton      *profilePlotterShortcut;
    //GLButton      *statusSimu;
    
	
    GLTextField   *facetSticking;
	
    GLCombo       *facetDesType;
	GLTextField   *facetDesTypeN;
    GLCombo       *facetRecType;
	GLLabel       *facetUseDesFileLabel;
	GLLabel       *modeLabel;
	
	GLLabel       *facetPumpingLabel;
	GLTextField   *facetPumping;	
    GLLabel       *facetSLabel;

	GLLabel       *facetDesRateLabel;
	GLTextField   *facetDesRate;

	GLLabel       *facetcoverageLabel;
	GLTextField   *facetcoverage;
	GLLabel       *facetcoveringLabel;
	GLTextField   *facetcovering;
	
    
	GLLabel       *facetTempLabel;
	GLTextField   *facetTemperature;
    GLLabel       *facetDLabel;
    GLLabel       *facetReLabel;
    GLToggle       *facetFILabel;
	GLToggle      *facetFIAreaLabel;
    //GLLabel       *facetMLabel;
	GLButton      *compACBtn;
	GLButton      *singleACBtn;

	GLButton      *profilePlotterBtn;
	GLButton      *texturePlotterBtn;
	GLButton      *textureScalingBtn;
	GLButton      *globalSettingsBtn;
	GLButton      *historyBtn;

	GLTitledPanel *inputPanel;
	GLTitledPanel *outputPanel;

    //Dialog
	ImportDesorption *importDesorption;
	TimeSettings     *timeSettings;
	Movement         *movement;
    FacetAdvParams   *facetAdvParams;
    FacetDetails     *facetDetails;
    Viewer3DSettings *viewer3DSettings;
    TextureScaling  *textureScaling;
	GlobalSettings	 *globalSettings;
	HistoryWin		*history;
    ProfilePlotter   *profilePlotter;
	PressureEvolution *pressureEvolution;
	TimewisePlotter  *timewisePlotter;
    TexturePlotter   *texturePlotter;
	OutgassingMap    *outgassingMap;
	MomentsEditor    *momentsEditor;
	ParameterEditor  *parameterEditor;
	char *nbF;

    // Testing
    //int     nbSt;
    //void LogProfile();
    void BuildPipe(double ratio,int steps=0);
	void EmptyGeometry();
	void CrashHandler(Error *e);
	void ExportHitBufferToFile(); //new function to export hit buffer for simulation on Linux HPC, added by Rudi.
	void ExportLoadBufferToFile(); //new function to export load buffer for simulation on Linux HPC, added by Rudi.
	void ImportHitBuffer(char *fName); //new function to import hit buffer of simulation on Linux HPC, added by Rudi.
	void ImportLoadBuffer(char *fName); //new function to import hit buffer of simulation on Linux HPC, added by Rudi.

protected:
	void LoadParameterCatalog();
    int  OneTimeSceneInit();
    int  RestoreDeviceObjects();
	int  InvalidateDeviceObjects();
    int  FrameMove();
    void ProcessMessage(GLComponent *src,int message);
};