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
#ifndef _FACETDETAILSH_
#define _FACETDETAILSH_

#include "GLApp/GLWindow.h"
#include "GLApp/GLList.h"
#include "GLApp/GLButton.h"
#include "GLApp/GLToggle.h"
#include "GLApp/GLTitledPanel.h"
#include "Worker.h"

class Facet;

#define NB_FDCOLUMN 27

class FacetDetails : public GLWindow {

public:

  // Construction
  FacetDetails();

  // Component method
  void Display(Worker *w);
  void Update();

  // Implementation
  void ProcessMessage(GLComponent *src,int message);
  void SetBounds(int x,int y,int w,int h);

private:

  char *GetCountStr(Facet *f);
  void UpdateTable();
  char *FormatCell(size_t idx,Facet *f,size_t mode);
  void PlaceComponents();

  Worker      *worker;
  GLList      *facetListD;

  GLTitledPanel *sPanel;          // Toggle panel
  GLToggle      *show[NB_FDCOLUMN];
  size_t            shown[NB_FDCOLUMN];

  GLButton    *checkAllButton;
  GLButton    *uncheckAllButton;
  GLButton    *dismissButton;
  GLButton	  *updateButton;

};

#endif /* _FACETDETAILSH_ */
