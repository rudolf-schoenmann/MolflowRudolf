/*
Program:     ContaminationFlow
Description: Monte Carlo simulator for satellite contanimation studies
Authors:     Rudolf Sch�nmann / Hoai My Van
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
#ifndef _GLOBALSETTINGSH_
#define _GLOBALSETTINGSH_

#include "GLApp/GLWindow.h"
class GLButton;
class GLTextField;
class GLLabel;
class GLToggle;
class GLTitledPanel;

class Worker;
class GLList;

class GlobalSettings : public GLWindow {

public:

  // Construction
  GlobalSettings(Worker *w);
  void UpdateOutgassing();
  void SMPUpdate();

  // Implementation
  void ProcessMessage(GLComponent *src,int message);
  void Update();

private:

	void RestartProc();
	  Worker      *worker;
  GLList      *processList;
  GLButton    *restartButton;
  GLButton    *maxButton;
  GLTextField *nbProcText;
  GLTextField *autoSaveText;
 

  int lastUpdate;
  //float lastCPUTime[MAX_PROCESS];
  //float lastCPULoad[MAX_PROCESS];

  GLToggle      *chkAntiAliasing;
  GLToggle      *chkWhiteBg;
  GLToggle		*leftHandedToggle;
  GLToggle      *chkSimuOnly;
  GLToggle      *chkCheckForUpdates;
  GLToggle      *chkAutoUpdateFormulas;
  GLToggle      *chkCompressSavedFiles;
  GLButton    *applyButton;
  GLButton    *recalcButton;

  GLTextField *outgassingText;
  GLTextField *influxText;
  GLTextField *gasMassText;
  GLToggle    *enableDecay;
  GLTextField *halfLifeText;
  GLToggle	*lowFluxToggle;
  GLButton    *lowFluxInfo;
  GLTextField *cutoffText;
  GLTextField* gasDiameterText;
};

#endif /* _GLOBALSETTINGSH_ */
