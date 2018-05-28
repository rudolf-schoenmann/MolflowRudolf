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
#include "PressureEvolution.h"
#include "GLApp/GLToolkit.h"
#include "GLApp/GLMessageBox.h"
#include "GLApp/GLToggle.h"
#include "GLApp/MathTools.h"
#include "GLApp/GLList.h"
#include "GLApp/GLChart/GLChart.h"
#include "GLApp/GLLabel.h"
#include "GLApp/GLCombo.h"
#include "GLApp/GLButton.h"
#include "GLApp/GLParser.h"
#include "GLApp/GLTextField.h"
#include "Geometry_shared.h"
#include "Facet_shared.h"
#include <math.h>

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

extern const char*profType[];

PressureEvolution::PressureEvolution() :GLWindow() {

	int wD = 750;
	int hD = 425;

	SetTitle("Pressure evolution plotter");
	SetIconfiable(true);
	nbView = 0;
	worker = NULL;
	lastUpdate = 0.0f;

	nbColors = 8;
	colors[0] = new GLCColor(); colors[0]->r = 255; colors[0]->g = 000; colors[0]->b = 055; //red
	colors[1] = new GLCColor(); colors[1]->r = 000; colors[1]->g = 000; colors[1]->b = 255; //blue
	colors[2] = new GLCColor(); colors[2]->r = 000; colors[2]->g = 204; colors[2]->b = 051; //green
	colors[3] = new GLCColor(); colors[3]->r = 000; colors[3]->g = 000; colors[3]->b = 000; //black
	colors[4] = new GLCColor(); colors[4]->r = 255; colors[4]->g = 153; colors[4]->b = 051; //orange
	colors[5] = new GLCColor(); colors[5]->r = 153; colors[5]->g = 204; colors[5]->b = 255; //light blue
	colors[6] = new GLCColor(); colors[6]->r = 153; colors[6]->g = 000; colors[6]->b = 102; //violet
	colors[7] = new GLCColor(); colors[7]->r = 255; colors[7]->g = 230; colors[7]->b = 005; //yellow

	selectedSlice = 49;

	chart = new GLChart(0);
	chart->SetBorder(BORDER_BEVEL_IN);
	chart->GetY1Axis()->SetGridVisible(true);
	chart->GetXAxis()->SetGridVisible(true);
	chart->GetXAxis()->SetName("Time (s)");
	chart->GetY1Axis()->SetAutoScale(true);
	chart->GetY2Axis()->SetAutoScale(true);
	chart->GetY1Axis()->SetAnnotation(VALUE_ANNO);
	chart->GetXAxis()->SetAnnotation(VALUE_ANNO);
	Add(chart);

	dismissButton = new GLButton(0, "Dismiss");
	Add(dismissButton);

	selButton = new GLButton(0, "Show Facet");
	Add(selButton);

	addButton = new GLButton(0, "Add curve");
	Add(addButton);

	removeButton = new GLButton(0, "Remove curve");
	Add(removeButton);

	resetButton = new GLButton(0, "Remove all");
	Add(resetButton);

	profCombo = new GLCombo(0);
	profCombo->SetEditable(true);
	Add(profCombo);

	logXToggle = new GLToggle(0, "Log X");
	Add(logXToggle);

	logYToggle = new GLToggle(0, "Log Y");
	Add(logYToggle);

	normLabel = new GLLabel("Normalize");
	Add(normLabel);

	normCombo = new GLCombo(0);
	normCombo->SetEditable(true);
	normCombo->SetSize(5);
	normCombo->SetValueAt(0, "None (raw data)");
	normCombo->SetValueAt(1, "Pressure (mbar)");
	normCombo->SetValueAt(2, "Density (1/m3)");
	normCombo->SetValueAt(3, "Speed (m/s)");
	normCombo->SetValueAt(4, "Angle (deg)");
	normCombo->SetSelectedIndex(1);
	Add(normCombo);

	modeCombo = new GLCombo(0);
	modeCombo->SetEditable(true);
	modeCombo->SetSize(2);
	modeCombo->SetValueAt(0, "Average");
	modeCombo->SetValueAt(1, "Slice #");
	modeCombo->SetSelectedIndex(0);
	Add(modeCombo);

	selectedSliceText = new GLTextField(0, "50");
	selectedSliceText->SetEditable(false);
	Add(selectedSliceText);

	label1 = new GLLabel("Show evolution of:");
	Add(label1);

	correctForGas = new GLToggle(0, "Surface->Volume conversion");
	correctForGas->SetVisible(false);
	Add(correctForGas);

	formulaText = new GLTextField(0, "");
	formulaText->SetEditable(true);
	Add(formulaText);

	formulaBtn = new GLButton(0, "-> Plot");
	Add(formulaBtn);

	// Center dialog
	int wS, hS;
	GLToolkit::GetScreenSize(&wS, &hS);
	int xD = (wS - wD) / 2;
	int yD = (hS - hD) / 2;
	SetBounds(xD, yD, wD, hD);
	SetResizable(true);
	SetMinimumSize(wD, 220);

	RestoreDeviceObjects();

}

void PressureEvolution::SetBounds(int x, int y, int w, int h) {

	chart->SetBounds(7, 5, w - 15, h - 110);
	profCombo->SetBounds(7, h - 95, 117, 19);
	selButton->SetBounds(130, h - 95, 80, 19);
	addButton->SetBounds(215, h - 95, 80, 19);
	removeButton->SetBounds(300, h - 95, 80, 19);
	resetButton->SetBounds(385, h - 95, 80, 19);

	logXToggle->SetBounds(w - 110, h - 95, 40, 19);
	logYToggle->SetBounds(w - 55, h - 95, 40, 19);
	label1->SetBounds(w - 265, h - 70, 90, 19);
	normLabel->SetBounds(7, h - 68, 50, 19);
	normCombo->SetBounds(60, h - 70, 105, 19);
	correctForGas->SetBounds(240, h - 70, 80, 19);
	modeCombo->SetBounds(w - 170, h - 70, 105, 19);
	selectedSliceText->SetBounds(w - 60, h - 70, 50, 19);

	formulaText->SetBounds(7, h - 45, 525, 19);
	formulaBtn->SetBounds(537, h - 45, 80, 19);
	dismissButton->SetBounds(w - 100, h - 45, 90, 19);

	GLWindow::SetBounds(x, y, w, h);

}

void PressureEvolution::Refresh() {

	if (!worker) return;

	Geometry *geom = worker->GetGeometry();
	size_t nb = geom->GetNbFacet();
	size_t nbProf = 0;
	for (size_t i = 0; i < nb; i++)
		if (geom->GetFacet(i)->sh.isProfile) nbProf++;
	profCombo->Clear();
	if (nbProf) profCombo->SetSize(nbProf);
	nbProf = 0;
	for (size_t i = 0; i < nb; i++) {
		Facet *f = geom->GetFacet(i);
		if (f->sh.isProfile) {
			char tmp[128];
			sprintf(tmp, "F#%zd %s", i + 1, profType[f->sh.profileType]);
			profCombo->SetValueAt(nbProf, tmp, (int)i);
			nbProf++;
		}
	}
	profCombo->SetSelectedIndex(nbProf?0:-1);
	//Remove profiles that aren't present anymore
	for (size_t v = 0; v < nbView; v++)
		if (views[v]->userData1 >= geom->GetNbFacet() || !geom->GetFacet(views[v]->userData1)->sh.isProfile) {
			chart->GetY1Axis()->RemoveDataView(views[v]);
			SAFE_DELETE(views[v]);
			for (size_t j = v; j < nbView - 1; j++) views[j] = views[j + 1];
			nbView--;
		}
	refreshViews();

}

void PressureEvolution::Display(Worker *w) {

	worker = w;
	Refresh();
	SetVisible(true);

}

void PressureEvolution::Update(float appTime, bool force) {

	if (!IsVisible() || IsIconic()) return;

	if (force) {
		refreshViews();
		lastUpdate = appTime;
		return;
	}

	if ((appTime - lastUpdate > 1.0f || force) && nbView) {
		if (worker->isRunning) refreshViews();
		lastUpdate = appTime;
	}

}

void PressureEvolution::plot() {

	GLParser *parser = new GLParser();
	parser->SetExpression(formulaText->GetText());
	if (!parser->Parse()) {
		GLMessageBox::Display(parser->GetErrorMsg(), "Error", GLDLG_OK, GLDLG_ICONERROR);
		SAFE_DELETE(parser);
		return;
	}

	int nbVar = parser->GetNbVariable();
	if (nbVar == 0) {
		GLMessageBox::Display("Variable 'x' not found", "Error", GLDLG_OK, GLDLG_ICONERROR);
		SAFE_DELETE(parser);
		return;
	}
	if (nbVar > 1) {
		GLMessageBox::Display("Too many variables or unknown constant", "Error", GLDLG_OK, GLDLG_ICONERROR);
		SAFE_DELETE(parser);
		return;
	}
	VLIST *var = parser->GetVariableAt(0);
	if (_stricmp(var->name, "x") != 0) {
		GLMessageBox::Display("Variable 'x' not found", "Error", GLDLG_OK, GLDLG_ICONERROR);
		SAFE_DELETE(parser);
		return;
	}

	Geometry *geom = worker->GetGeometry();
	GLDataView *v;

	// Check that view is not already added
	bool found = false;
	int i = 0;
	while (i < nbView && !found) {
		found = (views[i]->userData1 == -1);
		if (!found) i++;
	}

	if (found) {
		v = views[i];
		v->SetName(formulaText->GetText());
		v->Reset();
	}
	else {
		if (nbView < 50) {
			v = new GLDataView();
			v->SetName(formulaText->GetText());
			v->userData1 = -1;
			chart->GetY1Axis()->AddDataView(v);
			views[nbView] = v;
			nbView++;
		}
	}

	// Plot
	for (int p = 0; p < 1000; p++) {
		double x = (double)p;
		double y;
		var->value = x;
		parser->Evaluate(&y);
		v->Add(x, y, false);
	}
	v->CommitChange();

	delete parser;

}

void PressureEvolution::refreshViews() {

	// Lock during update
	BYTE *buffer = worker->GetHits();
	int displayMode = normCombo->GetSelectedIndex();
	if (!buffer) return;

	Geometry *geom = worker->GetGeometry();
	GlobalHitBuffer *gHits = (GlobalHitBuffer *)buffer;
	double nbDes = (double)gHits->total.hit.nbDesorbed;
	double scaleY;
	size_t facetHitsSize = (1 + worker->moments.size()) * sizeof(FacetHitBuffer);
	for (int i = 0; i < nbView; i++) {

		GLDataView *v = views[i];
		if (v->userData1 >= 0 && v->userData1 < geom->GetNbFacet()) {
			Facet *f = geom->GetFacet(v->userData1);
			//FacetHitBuffer *fCount = (FacetHitBuffer *)(buffer + f->sh.hitOffset);
			//double fnbDes = (double)fCount->hit.nbDesorbed;
			//double fnbHit = (double)fCount->hit.nbMCHit;
			v->Reset();
			for (int m = 1; m <= Min((int)worker->moments.size(), 10000); m++) { //max 10000 points
				ProfileSlice *profilePtr = (ProfileSlice *)(buffer + f->sh.hitOffset + facetHitsSize + m*sizeof(ProfileSlice)*PROFILE_SIZE);

				switch (displayMode) {
				case 0: { //Raw data 
					double val = 0;
					if (modeCombo->GetSelectedIndex() == 1) //plot one slice
						val = profilePtr[selectedSlice].countEquiv;
					else //plot sum
						for (int j = 0; j < PROFILE_SIZE; j++)
							val += profilePtr[j].countEquiv;
					v->Add(worker->moments[m - 1], val, false);
					break;
				}
				case 1: {//Pressure
					scaleY = 1.0 / nbDes / (f->sh.area / (double)PROFILE_SIZE*1E-4)* worker->gasMass / 1000 / 6E23 * 0.0100; //0.01: Pa->mbar
					scaleY *= worker->totalDesorbedMolecules / worker->timeWindowSize;
					if (f->sh.is2sided) scaleY *= 0.5;
					double val = 0.0;
					if (modeCombo->GetSelectedIndex() == 1) //plot one slice
						val = profilePtr[selectedSlice].sum_v_ort*scaleY;
					else {//plot avg
						for (int j = 0; j < PROFILE_SIZE; j++)
							val += profilePtr[j].sum_v_ort;
						val *= scaleY / (double)PROFILE_SIZE;
					}
					v->Add(worker->moments[m - 1], val, false);
					break;
				}
				case 2: {//Particle density
					scaleY = 1.0 / nbDes / (f->GetArea() / (double)PROFILE_SIZE*1E-4);
					scaleY *= worker->totalDesorbedMolecules / worker->timeWindowSize;
					scaleY *= f->DensityCorrection();

					double val = 0.0;
					if (modeCombo->GetSelectedIndex() == 1) //plot one slice
						val = profilePtr[selectedSlice].sum_1_per_ort_velocity*scaleY;
					else {//plot avg
						for (int j = 0; j < PROFILE_SIZE; j++)
							val += profilePtr[j].sum_1_per_ort_velocity;
						val *= scaleY / (double)PROFILE_SIZE;
					}
					v->Add(worker->moments[m - 1], val, false);
					break;
				}
				case 3: {//Velocity
					double scaleX = f->sh.maxSpeed / (double)PROFILE_SIZE;
					double sum = 0.0;
					double weighedSum = 0.0;
					double val;
					for (int j = 0; j < PROFILE_SIZE; j++) {
						if (!correctForGas->GetState())
							val = (double)profilePtr[j].countEquiv; //weighed sum
						else
							val = (double)profilePtr[j].countEquiv / (((double)j + 0.5));
						sum += val;
						weighedSum += ((double)j+0.5)*val;
					}
					weighedSum *= scaleX / sum;
					v->Add(worker->moments[m - 1], weighedSum, false);
					break; }
				case 4: {//Angle (deg)
					double scaleX = 90.0 / (double)PROFILE_SIZE;
					double sum = 0.0;
					double weighedSum = 0.0;
					double val;
					for (int j = 0; j < PROFILE_SIZE; j++) {
						if (!correctForGas->GetState())
							val = (double)profilePtr[j].countEquiv; //weighed sum
						else
							val = (double)profilePtr[j].countEquiv / sin(((double)j + 0.5)*PI / 2.0 / (double)PROFILE_SIZE);
						sum += val;
						weighedSum += ((double)j+0.5)*val;
					}
					weighedSum *= scaleX / sum;
					v->Add(worker->moments[m - 1], weighedSum, false);
					break; }
				}
			}
			v->CommitChange();
		}
		else {
			if (v->userData1 == -2 && nbDes != 0.0) {

				// Volatile profile
				v->Reset();
				size_t nb = geom->GetNbFacet();
				for (size_t j = 0; j < nb; j++) {
					Facet *f = geom->GetFacet(j);
					if (f->sh.isVolatile) {
						FacetHitBuffer *fCount = (FacetHitBuffer *)(buffer + f->sh.hitOffset);
						double z = geom->GetVertex(f->indices[0])->z;
						v->Add(z, fCount->hit.nbAbsEquiv / nbDes, false);
					}
				}
				// Last
				Facet *f = geom->GetFacet(28);
				FacetHitBuffer *fCount = (FacetHitBuffer *)(buffer + f->sh.hitOffset);
				double fnbAbs = fCount->hit.nbAbsEquiv;
				v->Add(1000.0, fnbAbs / nbDes, false);
				v->CommitChange();

				//v->Reset();
				//for(int j=0;j<BOUNCEMAX && nbAbs;j++)
				//  v->Add((double)j,(double)gHits->wallHits[j]/nbAbs);

			}
		}

	}

	worker->ReleaseHits();

}

void PressureEvolution::addView(int facet) {

	char tmp[128];
	Geometry *geom = worker->GetGeometry();

	// Check that view is not already added
	bool found = false;
	int i = 0;
	while (i < nbView && !found) {
		found = (views[i]->userData1 == facet);
		if (!found) i++;
	}
	if (worker->moments.size() > 10000) {
		GLMessageBox::Display("Only the first 10000 moments will be plotted", "Error", GLDLG_OK, GLDLG_ICONWARNING);
	}
	if (found) {
		GLMessageBox::Display("Profile already plotted", "Error", GLDLG_OK, GLDLG_ICONERROR);
		return;
	}
	if (nbView < 50) {
		Facet *f = geom->GetFacet(facet);
		GLDataView *v = new GLDataView();
		//sprintf(tmp, "F#%d %s", facet + 1, profType[f->sh.profileType]);
		sprintf(tmp, "F#%d", facet + 1);
		v->SetName(tmp);
		v->SetViewType(TYPE_BAR);
		v->SetMarker(MARKER_DOT);
		v->SetColor(*colors[nbView%nbColors]);
		v->SetMarkerColor(*colors[nbView%nbColors]);
		v->userData1 = facet;
		chart->GetY1Axis()->AddDataView(v);
		views[nbView] = v;
		nbView++;
	}

}

void PressureEvolution::remView(int facet) {

	Geometry *geom = worker->GetGeometry();

	bool found = false;
	int i = 0;
	while (i < nbView && !found) {
		found = (views[i]->userData1 == facet);
		if (!found) i++;
	}
	if (!found) {
		GLMessageBox::Display("Profile not plotted", "Error", GLDLG_OK, GLDLG_ICONERROR);
		return;
	}
	chart->GetY1Axis()->RemoveDataView(views[i]);
	SAFE_DELETE(views[i]);
	for (int j = i; j < nbView - 1; j++) views[j] = views[j + 1];
	nbView--;

}

void PressureEvolution::Reset() {

	chart->GetY1Axis()->ClearDataView();
	for (int i = 0; i < nbView; i++) SAFE_DELETE(views[i]);
	nbView = 0;

}

void PressureEvolution::ProcessMessage(GLComponent *src, int message) {
	Geometry *geom = worker->GetGeometry();
	switch (message) {
	case MSG_BUTTON:
		if (src == dismissButton) {
			SetVisible(false);
		}
		else if (src == selButton) {
			int idx = profCombo->GetSelectedIndex();
			if (idx >= 0) {
				geom->UnselectAll();
				geom->GetFacet(profCombo->GetUserValueAt(idx))->selected = true;
				geom->UpdateSelection();
				mApp->UpdateFacetParams(true);
				mApp->facetList->SetSelectedRow(profCombo->GetUserValueAt(idx));
				mApp->facetList->ScrollToVisible(profCombo->GetUserValueAt(idx), 1, true);
			}
		}
		else if (src == addButton) {
			int idx = profCombo->GetSelectedIndex();
			if (idx >= 0) addView(profCombo->GetUserValueAt(idx));
			refreshViews();
		}
		else if (src == removeButton) {
			int idx = profCombo->GetSelectedIndex();
			if (idx >= 0) remView(profCombo->GetUserValueAt(idx));
			refreshViews();
		}
		else if (src == resetButton) {
			Reset();
		}
		else if (src == formulaBtn) {
			plot();
		}
		break;
	case MSG_COMBO:
		if (src == normCombo) {
			int normMode = normCombo->GetSelectedIndex();
			correctForGas->SetVisible(normMode == 3 || normMode == 4);
			refreshViews();
			if (normCombo->GetSelectedIndex() == 0) {//pressure selected
				modeCombo->SetValueAt(0, "Sum");
				if (modeCombo->GetSelectedIndex() == 0) modeCombo->SetSelectedIndex(0); //re-select to update text
			}
			else { //something else selected
				modeCombo->SetValueAt(0, "Average");
				if (modeCombo->GetSelectedIndex() == 0) modeCombo->SetSelectedIndex(0); //re-select to update text
			}
			modeCombo->SetEditable(!(normMode == 3 || normMode == 4));
		}
		else if (src == modeCombo) {
			selectedSliceText->SetEditable(modeCombo->GetSelectedIndex() == 1);
		}
		refreshViews();
		break;
	case MSG_TEXT: //enter pressed
		if (src == selectedSliceText) {
			int sliceNum;
			if (!selectedSliceText->GetNumberInt(&sliceNum)) {
				GLMessageBox::Display("Invalid slice number", "Error", GLDLG_OK, GLDLG_ICONERROR);
				return;
			}
			if (!(sliceNum > 0 && sliceNum <= 100)) {
				GLMessageBox::Display("Choose a slice between 1 and 100", "Error", GLDLG_OK, GLDLG_ICONERROR);
				return;
			}
			//valid input
			selectedSlice = sliceNum - 1;
			refreshViews();
		}
		break;
	case MSG_TOGGLE:
		if (src == logXToggle) {
			chart->GetXAxis()->SetScale(logXToggle->GetState());
		}
		else if (src == logYToggle) {
			chart->GetY1Axis()->SetScale(logYToggle->GetState());
		}
		else if (src == correctForGas) {
			refreshViews();
		}
		break;
	}

	GLWindow::ProcessMessage(src, message);

}

