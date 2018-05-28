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
#include "TexturePlotter.h"
#include "GLApp/GLToolkit.h"
#include "GLApp/GLMessageBox.h"
#include "GLApp/GLFileBox.h"
#include "GLApp/GLLabel.h"
#include "GLApp/GLButton.h"
#include "GLApp/GLList.h"
#include "GLApp/GLCombo.h"
#include "GLApp\GLToggle.h"
#include "Geometry_shared.h"
#include "Facet_shared.h"
#ifdef MOLFLOW
#include "MolFlow.h"
#endif

#ifdef SYNRAD
#include "SynRad.h"
#endif

#ifdef MOLFLOW
extern MolFlow *mApp;
#endif

#ifdef SYNRAD
extern SynRad*mApp;
#endif

static const char *fileFilters = "Text files\0*.txt";
static const int   nbFilter = sizeof(fileFilters) / (2 * sizeof(char *));

TexturePlotter::TexturePlotter() :GLWindow() {

	int wD = 500;
	int hD = 300;
	lastUpdate = 0.0f;
	strcpy(currentDir, ".");

	SetTitle("Texture plotter");
	SetResizable(true);
	SetIconfiable(true);
	SetMinimumSize(wD, hD);

	mapList = new GLList(0);
	mapList->SetColumnLabelVisible(true);
	mapList->SetRowLabelVisible(true);
	mapList->SetAutoColumnLabel(true);
	mapList->SetAutoRowLabel(true);
	mapList->SetRowLabelMargin(20);
	mapList->SetGrid(true);
	mapList->SetSelectionMode(BOX_CELL);
	mapList->SetCornerLabel("\202\\\201");
	Add(mapList);

	viewLabel = new GLLabel("View:");
	Add(viewLabel);
	viewCombo = new GLCombo(0);
	viewCombo->SetSize(9);
	viewCombo->SetValueAt(0, "Cell area (cm\262)");
	viewCombo->SetValueAt(1, "# of MC hits");
	viewCombo->SetValueAt(2, "Impingement rate [1/m\262/sec]");
	viewCombo->SetValueAt(3, "Particle density [1/m3]");
	viewCombo->SetValueAt(4, "Gas density [kg/m3]");
	viewCombo->SetValueAt(5, "Pressure [mbar]");
	viewCombo->SetValueAt(6, "Avg.speed estimate[m/s]");
	viewCombo->SetValueAt(7, "Incident velocity vector[m/s]");
	viewCombo->SetValueAt(8, "# of velocity vectors");

	viewCombo->SetSelectedIndex(5); //Pressure by default
	Add(viewCombo);

	saveButton = new GLButton(0, "Save");
	Add(saveButton);

	sizeButton = new GLButton(0, "Autosize");
	Add(sizeButton);

	maxButton = new GLButton(0, "Find Max.");
	Add(maxButton);

	cancelButton = new GLButton(0, "Dismiss");
	Add(cancelButton);

	autoSizeOnUpdate = new GLToggle(0, "Autosize on every update (disable for smooth scrolling)");
	autoSizeOnUpdate->SetState(true);
	Add(autoSizeOnUpdate);

	// Center dialog
	int wS, hS;
	GLToolkit::GetScreenSize(&wS, &hS);
	int xD = (wS - wD) / 2;
	int yD = (hS - hD) / 2;
	SetBounds(xD, yD, wD, hD);

	RestoreDeviceObjects();

	worker = NULL;

}

void TexturePlotter::PlaceComponents() {

	mapList->SetBounds(5, 5, width - 15, height - 80);
	saveButton->SetBounds(10, height - 70, 70, 19);
	sizeButton->SetBounds(10, height - 45, 70, 19);
	autoSizeOnUpdate->SetBounds(90, height - 45, 120, 19);
	maxButton->SetBounds(90, height - 70, 70, 19);
	viewLabel->SetBounds(320, height - 70, 30, 19);
	viewCombo->SetBounds(350, height - 70, 130, 19);
	cancelButton->SetBounds(width - 90, height - 45, 80, 19);

}

void TexturePlotter::SetBounds(int x, int y, int w, int h) {

	GLWindow::SetBounds(x, y, w, h);
	PlaceComponents();

}

void TexturePlotter::GetSelected() {

	if (!worker) return;

	Geometry *geom = worker->GetGeometry();
	selFacet = NULL;
	int i = 0;
	size_t nb = geom->GetNbFacet();
	while (!selFacet && i < nb) {
		if (geom->GetFacet(i)->selected) selFacet = geom->GetFacet(i);
		if (!selFacet) i++;
	}

	char tmp[32];
	sprintf(tmp, "Texture plotter [Facet #%d]", i + 1);
	SetTitle(tmp);

}

void TexturePlotter::Update(float appTime, bool force) {

	if (!IsVisible()) return;

	if (force) {
		UpdateTable();
		lastUpdate = appTime;
		return;
	}

	if ((appTime - lastUpdate > 1.0f)) {
		if (worker->isRunning) UpdateTable();
		lastUpdate = appTime;
	}

}

void TexturePlotter::UpdateTable() {
	size_t nbMoments = mApp->worker.moments.size();
	size_t facetHitsSize = (1 + nbMoments) * sizeof(FacetHitBuffer);
	maxValue = 0.0f;
	//double scale;
	GetSelected();
	if (!selFacet || !selFacet->cellPropertiesIds) {
		mapList->Clear();
		return;
	}

	//SHELEM *mesh = selFacet->mesh;
	if (selFacet->cellPropertiesIds) {

		char tmp[256];
		size_t w = selFacet->sh.texWidth;
		size_t h = selFacet->sh.texHeight;
		mapList->SetSize(w, h);
		mapList->SetAllColumnAlign(ALIGN_CENTER);

		int mode = viewCombo->GetSelectedIndex();

		switch (mode) {

		case 0: {// Cell area
			for (size_t i = 0; i < w; i++) {
				for (size_t j = 0; j < h; j++) {
					float val = selFacet->GetMeshArea(i + j*w);
					sprintf(tmp, "%g", val);
					if (val > maxValue) {
						maxValue = val;
						maxX = i; maxY = j;
					}
					mapList->SetValueAt(i, j, tmp);
				}
			}
			break; }

		case 1: {// MC Hits

					 // Lock during update
			BYTE *buffer = worker->GetHits();
			try {
				if (buffer) {
					GlobalHitBuffer *shGHit = (GlobalHitBuffer *)buffer;
					size_t profSize = (selFacet->sh.isProfile) ? (PROFILE_SIZE * sizeof(ProfileSlice)*(1 + nbMoments)) : 0;
					TextureCell *texture = (TextureCell *)((BYTE *)buffer + (selFacet->sh.hitOffset + facetHitsSize + profSize + mApp->worker.displayedMoment*w*h * sizeof(TextureCell)));
					for (size_t i = 0; i < w; i++) {
						for (size_t j = 0; j < h; j++) {
							//int tSize = selFacet->sh.texWidth*selFacet->sh.texHeight;

							double val = texture[i + j*w].countEquiv;
							if (val > maxValue) {
								maxValue = (double)val;
								maxX = i; maxY = j;
							}
							sprintf(tmp, "%g", val);
							mapList->SetValueAt(i, j, tmp);
						}
					}
					worker->ReleaseHits();
				}
			}
			catch (...) {
				worker->ReleaseHits();

			}

			break; }

		case 2: {// Impingement rate

					 // Lock during update
			BYTE *buffer = worker->GetHits();
			try {
				if (buffer) {
					GlobalHitBuffer *shGHit = (GlobalHitBuffer *)buffer;
					size_t profSize = (selFacet->sh.isProfile) ? (PROFILE_SIZE * sizeof(ProfileSlice)*(1 + nbMoments)) : 0;
					TextureCell *texture = (TextureCell *)((BYTE *)buffer + (selFacet->sh.hitOffset + facetHitsSize + profSize + mApp->worker.displayedMoment*w*h * sizeof(TextureCell)));
					double dCoef =1E4; //1E4: conversion m2->cm2
					/*if (shGHit->sMode == MC_MODE) dCoef *= ((mApp->worker.displayedMoment == 0) ? 1.0 : ((worker->desorptionStopTime - worker->desorptionStartTime)
						/ worker->timeWindowSize));*/
					if (shGHit->sMode == MC_MODE) dCoef *= mApp->worker.GetMoleculesPerTP(worker->displayedMoment);
					for (size_t i = 0; i < w; i++) {
						for (size_t j = 0; j < h; j++) {
							double area = (selFacet->GetMeshArea(i + j*w,true)); if (area == 0.0) area = 1.0;
							double val = texture[i + j*w].countEquiv / area*dCoef;
							if (val > maxValue) {
								maxValue = val;
								maxX = i; maxY = j;
							}
							sprintf(tmp, "%g", val);
							mapList->SetValueAt(i, j, tmp);
						}
					}
					worker->ReleaseHits();
				}
			}
			catch (...) { //incorrect hits reference
				worker->ReleaseHits();
			}

			break; }

		case 3: {// Particle density [1/m3]

					 // Lock during update
			BYTE *buffer = worker->GetHits();
			try {
				if (buffer) {
					GlobalHitBuffer *shGHit = (GlobalHitBuffer *)buffer;
					size_t profSize = (selFacet->sh.isProfile) ? (PROFILE_SIZE * sizeof(ProfileSlice)*(1 + nbMoments)) : 0;
					TextureCell *texture = (TextureCell *)((BYTE *)buffer + (selFacet->sh.hitOffset + facetHitsSize + profSize + mApp->worker.displayedMoment*w*h * sizeof(TextureCell)));
					double dCoef =1E4 * selFacet->DensityCorrection();   //1E4 m2 -> cm2
					
					if (shGHit->sMode == MC_MODE) dCoef *= mApp->worker.GetMoleculesPerTP(worker->displayedMoment);
					for (size_t i = 0; i < w; i++) {
						for (size_t j = 0; j < h; j++) {

							/*double v_avg = 2.0*(double)texture[i + j*w].count / texture[i + j*w].sum_1_per_ort_velocity;
							double imp_rate = texture[i + j*w].count / (selFacet->mesh[i + j*w].area*(selFacet->sh.is2sided ? 2.0 : 1.0))*dCoef;
							double rho = 4.0*imp_rate / v_avg;*/
							double rho = texture[i + j*w].sum_1_per_ort_velocity / selFacet->GetMeshArea(i + j*w,true)*dCoef;
							if (rho > maxValue) {
								maxValue = rho;
								maxX = i; maxY = j;
							}

							sprintf(tmp, "%g", rho);
							mapList->SetValueAt(i, j, tmp);
						}
					}
					worker->ReleaseHits();
				}
			}
			catch (...) { //incorrect hits reference
				worker->ReleaseHits();
			}

			break; }

		case 4: {// Gas density [kg/m3]

					 // Lock during update
			BYTE *buffer = worker->GetHits();
			try {
				if (buffer) {
					GlobalHitBuffer *shGHit = (GlobalHitBuffer *)buffer;
					size_t profSize = (selFacet->sh.isProfile) ? (PROFILE_SIZE * sizeof(ProfileSlice)*(1 + nbMoments)) : 0;
					TextureCell *texture = (TextureCell *)((BYTE *)buffer + (selFacet->sh.hitOffset + facetHitsSize + profSize + mApp->worker.displayedMoment*w*h * sizeof(TextureCell)));
					//float dCoef = (float)totalOutgassing / 8.31 * gasMass / 100 * MAGIC_CORRECTION_FACTOR;
					double dCoef = 1E4 * selFacet->DensityCorrection();

					if (shGHit->sMode == MC_MODE) dCoef *= worker->GetMoleculesPerTP(worker->displayedMoment);
					for (size_t i = 0; i < w; i++) {
						for (size_t j = 0; j < h; j++) {

							/*double v_avg = 2.0*(double)texture[i + j*w].count / texture[i + j*w].sum_1_per_ort_velocity;
							double imp_rate = texture[i + j*w].count / (selFacet->mesh[i + j*w].area*(selFacet->sh.is2sided ? 2.0 : 1.0))*dCoef;
							double rho = 4.0*imp_rate / v_avg;*/
							double rho = texture[i + j*w].sum_1_per_ort_velocity / selFacet->GetMeshArea(i + j*w,true)*dCoef;
							double rho_mass = rho*worker->gasMass / 1000.0 / 6E23;
							if (rho_mass > maxValue) {
								maxValue = rho_mass;
								maxX = i; maxY = j;
							}

							sprintf(tmp, "%g", rho_mass);
							mapList->SetValueAt(i, j, tmp);
						}
					}
					worker->ReleaseHits();
				}
			}
			catch (...) { //incorrect hits reference
				worker->ReleaseHits();
			}

			break; }

		case 5: {// Pressure

					 // Lock during update
			BYTE *buffer = worker->GetHits();
			try {
				if (buffer) {
					GlobalHitBuffer *shGHit = (GlobalHitBuffer *)buffer;
					size_t profSize = (selFacet->sh.isProfile) ? (PROFILE_SIZE * sizeof(ProfileSlice)*(1 + nbMoments)) : 0;
					TextureCell *texture = (TextureCell *)((BYTE *)buffer + (selFacet->sh.hitOffset + facetHitsSize + profSize + mApp->worker.displayedMoment*w*h * sizeof(TextureCell)));
					double dCoef = 1E4 * (worker->gasMass / 1000 / 6E23) * 0.0100;  //1E4 is conversion from m2 to cm2; 0.01 is Pa->mbar
					
					if (shGHit->sMode == MC_MODE) dCoef *= worker->GetMoleculesPerTP(worker->displayedMoment);
					for (size_t i = 0; i < w; i++) {
						for (size_t j = 0; j < h; j++) {

							double p = texture[i + j*w].sum_v_ort_per_area*dCoef;
							if (p > maxValue) {
								maxValue = p;
								maxX = i; maxY = j;
							}

							sprintf(tmp, "%g", p);

							mapList->SetValueAt(i, j, tmp);
						}
					}
					worker->ReleaseHits();
				}
			}
			catch (...) { //incorrect hits reference
				worker->ReleaseHits();
			}

			break; }

		case 6: {// Average gas velocity [m/s]

					// Lock during update
			BYTE *buffer = worker->GetHits();
			try {
				if (buffer) {
					GlobalHitBuffer *shGHit = (GlobalHitBuffer *)buffer;
					size_t profSize = (selFacet->sh.isProfile) ? (PROFILE_SIZE * sizeof(ProfileSlice)*(1 + nbMoments)) : 0;
					TextureCell *texture = (TextureCell *)((BYTE *)buffer + (selFacet->sh.hitOffset + facetHitsSize + profSize + mApp->worker.displayedMoment*w*h * sizeof(TextureCell)));
					for (size_t i = 0; i < w; i++) {
						for (size_t j = 0; j < h; j++) {
							size_t tSize = selFacet->sh.texWidth*selFacet->sh.texHeight;
							double val = 4.0*texture[i + j*w].countEquiv / texture[i + j*w].sum_1_per_ort_velocity;
							if (val > maxValue) {
								maxValue = val;
								maxX = i; maxY = j;
							}
							sprintf(tmp, "%g", val);
							mapList->SetValueAt(i, j, tmp);
						}
					}
					worker->ReleaseHits();
				}
			}
			catch (...) {
				worker->ReleaseHits();

			}

			break; }

		case 7: {// Gas velocity vector
			for (size_t i = 0; i < w; i++) {
				for (size_t j = 0; j < h; j++) {
					if (selFacet->dirCache) {
						sprintf(tmp, "%g,%g,%g",
							selFacet->dirCache[i + j*w].dir.x / (double)selFacet->dirCache[i + j*w].count,
							selFacet->dirCache[i + j*w].dir.y / (double)selFacet->dirCache[i + j*w].count,
							selFacet->dirCache[i + j*w].dir.z / (double)selFacet->dirCache[i + j*w].count);
						double vsum = (selFacet->dirCache[i + j*w].dir.x*selFacet->dirCache[i + j*w].dir.x +
							selFacet->dirCache[i + j*w].dir.y*selFacet->dirCache[i + j*w].dir.y +
							selFacet->dirCache[i + j*w].dir.z + selFacet->dirCache[i + j*w].dir.z);
						if (vsum > maxValue) {
							maxValue = vsum;
							maxX = i; maxY = j;
						}
					}
					else {
						sprintf(tmp, "Direction not recorded");
					}
					mapList->SetValueAt(i, j, tmp);
				}
			}
			break; }

		case 8: {// # of velocity vectors
			for (size_t i = 0; i < w; i++) {
				for (size_t j = 0; j < h; j++) {
					if (selFacet->dirCache) {
						llong val = selFacet->dirCache[i + j*w].count;
						if (val > maxValue) {
							maxValue = (double)val;
							maxX = i; maxY = j;
						}
						sprintf(tmp, "%I64d", val);
					}
					else {
						sprintf(tmp, "Direction not recorded");
					}
					mapList->SetValueAt(i, j, tmp);
				}
			}
			break; }
		}

	}
	if (autoSizeOnUpdate->GetState()) mapList->AutoSizeColumn();
}

void TexturePlotter::Display(Worker *w) {

	worker = w;
	UpdateTable();
	SetVisible(true);

}

void TexturePlotter::Close() {
	worker = NULL;
	if (selFacet) selFacet->UnselectElem();
	mapList->Clear();
}

void TexturePlotter::SaveFile() {

	if (!selFacet) return;

	FILENAME *fn = GLFileBox::SaveFile(currentDir, NULL, "Save File", fileFilters, nbFilter);

	if (fn) {

		size_t u, v, wu, wv;
		if (!mapList->GetSelectionBox(&u, &v, &wu, &wv)) {
			u = 0;
			v = 0;
			wu = mapList->GetNbRow();
			wv = mapList->GetNbColumn();
		}

		// Save tab separated text
		FILE *f = fopen(fn->fullName, "w");

		if (f == NULL) {
			char errMsg[512];
			sprintf(errMsg, "Cannot open file\nFile:%s", fn->fullName);
			GLMessageBox::Display(errMsg, "Error", GLDLG_OK, GLDLG_ICONERROR);
			return;
		}

		for (size_t i = u; i < u + wu; i++) {
			for (size_t j = v; j < v + wv; j++) {
				char *str = mapList->GetValueAt(j, i);
				if (str) fprintf(f, "%s", str);
				if (j < v + wv - 1)
					fprintf(f, "\t");
			}
			fprintf(f, "\r\n");
		}
		fclose(f);

	}

}

void TexturePlotter::ProcessMessage(GLComponent *src, int message) {

	switch (message) {

	case MSG_CLOSE:
		Close();
		break;

	case MSG_BUTTON:
		if (src == cancelButton) {
			Close();
			GLWindow::ProcessMessage(NULL, MSG_CLOSE);
		}
		else if (src == sizeButton) {
			mapList->AutoSizeColumn();
		}
		else if (src == saveButton) {
			SaveFile();
		}
		else if (src == maxButton) {
			size_t u, v, wu, wv;
			mapList->SetSelectedCell(maxX, maxY);
			if (mapList->GetSelectionBox(&v, &u, &wv, &wu))
				selFacet->SelectElem(u, v, wu, wv);
		}
		break;

	case MSG_LIST:
		if (src == mapList) {
			size_t u, v, wu, wv;
			if (mapList->GetSelectionBox(&v, &u, &wv, &wu))
				selFacet->SelectElem(u, v, wu, wv);
		}
		break;

	case MSG_COMBO:
		if (src == viewCombo) {
			UpdateTable();
			maxButton->SetEnabled(true);
			//maxButton->SetEnabled(viewCombo->GetSelectedIndex()!=2);
		}
		break;

	}

	GLWindow::ProcessMessage(src, message);
}
