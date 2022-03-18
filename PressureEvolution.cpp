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
#include "GLApp/GLTitledPanel.h"
#include "GLApp/GLFileBox.h"
#include "Geometry_shared.h"
#include <math.h>

#include <fstream>
#include <iomanip>

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

PressureEvolution::PressureEvolution(Worker *w) :GLWindow() {

	int wD = 750;
	int hD = 800;//375;

	SetTitle("Pressure evolution plotter");
	SetIconfiable(true);
	lastUpdate = 0.0f;

	colors = {
		GLColor(255,000,055), //red
		GLColor(000,000,255), //blue
		GLColor(000,204,051), //green
		GLColor(000,000,000), //black
		GLColor(255,153,051), //orange
		GLColor(153,204,255), //light blue
		GLColor(153,000,102), //violet
		GLColor(255,230,005)  //yellow
	};

	panelGraph = new GLTitledPanel("Graph");
	Add(panelGraph);

	panelTable = new GLTitledPanel("Table");
	Add(panelTable);

	chart = new GLChart(0);
	chart->SetBorder(BORDER_BEVEL_IN);
	chart->GetY1Axis()->SetGridVisible(true);
	chart->GetXAxis()->SetGridVisible(true);
	chart->GetXAxis()->SetName("Time (s)");
	chart->GetY1Axis()->SetAutoScale(true);
	chart->GetY2Axis()->SetAutoScale(true);
	chart->GetY1Axis()->SetAnnotation(VALUE_ANNO);
	chart->GetXAxis()->SetAnnotation(VALUE_ANNO);
	panelGraph->Add(chart);
	historyList = new GLList(0);
	historyList->SetHScrollVisible(true);
	historyList->SetSize(0,0);
	historyList->SetColumnLabelVisible(true);
	historyList->SetSelectionMode(1);
	panelTable->Add(historyList);
	
	exportButton = new GLButton(0, "Export Table");
	panelTable->Add(exportButton);
	refreshButton = new GLButton(0, "Refresh");
	panelTable->Add(refreshButton);
	onlyViewsButton = new GLButton(0, "Show only selected facets");
	panelTable->Add(onlyViewsButton);
	onlyViews = false;

	selButton = new GLButton(0, "<-Show Facet");
	Add(selButton);

	removeButton = new GLButton(0, "Remove");
	Add(removeButton);

	removeAllButton = new GLButton(0, "Remove all");
	Add(removeAllButton);

	addButton = new GLButton(0, "Add selected facet");
	Add(addButton);

	profCombo = new GLCombo(0);
	Add(profCombo);

	logXToggle = new GLToggle(0, "Log X");
	Add(logXToggle);

	logYToggle = new GLToggle(0, "Log Y");
	Add(logYToggle);

	normLabel = new GLLabel("Y scale:");
	Add(normLabel);

	yScaleCombo = new GLCombo(0);
	yScaleCombo->SetEditable(true);
	yScaleCombo->SetSize(5);
	yScaleCombo->SetValueAt(0, "MC Hits");
	yScaleCombo->SetValueAt(1, "Equiv. hits");
	yScaleCombo->SetValueAt(2, "Pressure (mbar)");
	yScaleCombo->SetValueAt(3, "Density (1/m3)");
	yScaleCombo->SetValueAt(4, "Imp.rate (1/s/m2)");
	yScaleCombo->SetSelectedIndex(2); //Pressure by default
	Add(yScaleCombo);

	// Center dialog
	int wS, hS;
	GLToolkit::GetScreenSize(&wS, &hS);
	int xD = (wS - wD) / 2;
	int yD = (hS - hD) / 2;
	SetBounds(xD, yD, wD, hD);
	SetResizable(true);
	SetMinimumSize(wD, 220);

	RestoreDeviceObjects();

	worker = w;

	Geometry *geom = worker->GetGeometry();
	nb_Facets = geom->GetNbFacet();
	Refresh();

}

void PressureEvolution::SetBounds(int x, int y, int w, int h) {

	int mid = (h - 60) / 2;
	int offset = 80;

	panelGraph->SetBounds(5, 10, w - 10, mid+offset);
	panelTable->SetBounds(5, mid+10+offset, w - 10, mid-offset);

	panelGraph->SetCompBounds(chart, 10, 20, w-30, mid-30+offset);

	panelTable->SetCompBounds(historyList, 10, 20, w-30, mid-50-offset);
	panelTable->SetCompBounds(exportButton, 10, mid-offset-23, 80, 19);
	panelTable->SetCompBounds(refreshButton, 100, mid - offset - 23, 60, 19);
	panelTable->SetCompBounds(onlyViewsButton, 170, mid - offset - 23, 150, 19);

	//historyList->SetColumnWidthForAll((w - 40) / (nb_Facets+1));

	profCombo->SetBounds(7, h - 45, 117, 19);
	selButton->SetBounds(130, h - 45, 80, 19);
	removeButton->SetBounds(215, h - 45, 60, 19);
	removeAllButton->SetBounds(280, h - 45, 70, 19);
	addButton->SetBounds(370, h - 45, 110, 19);

	logXToggle->SetBounds(w - 105, h - 45, 40, 19);
	logYToggle->SetBounds(w - 55, h - 45, 40, 19);
	normLabel->SetBounds(490, h - 42, 50, 19);
	yScaleCombo->SetBounds(535, h - 45, 105, 19);

	GLWindow::SetBounds(x, y, w, h);

}

void PressureEvolution::Refresh() {
	//Rebuilds combo and calls refreshviews

	Geometry *geom = worker->GetGeometry();

	//Remove views that aren't present anymore
	for (auto i = views.begin();i!=views.end();) {
		if ((*i)->userData1 >= geom->GetNbFacet()) {
			chart->GetY1Axis()->RemoveDataView(*i);
			views.erase(i);
		}
		else {
			i++;
		}
	}

	//Construct combo
	profCombo->SetSize(views.size());
	size_t nbProf = 0;
	for (auto& v : views) {
		profCombo->SetValueAt(nbProf++, v->GetName(), v->userData1);
	}
	profCombo->SetSelectedIndex(nbProf ? (int)nbProf-1 : -1);

	//Refresh chart
	refreshChart();
}



void PressureEvolution::Update(float appTime, bool force) {
	//Calls refreshChart if needed
	if (!IsVisible() || IsIconic()) return;

	if (force) {
		refreshChart();
		lastUpdate = appTime;
		return;
	}

	if ((appTime - lastUpdate > 1.0f || force) && views.size() > 0) {
		if (worker->isRunning) refreshChart();
		lastUpdate = appTime;
	}

}

void PressureEvolution::exportList() {
	FILENAME *fn = GLFileBox::SaveFile(mApp->currentDir, NULL, "Save File", fileTexFilters, 0);

	if (fn) {
		try {
			char tmp[512];
			char *fileName = fn->fullName;

			std::ofstream f;
			char *ext, *dir;
			dir = strrchr(fileName, '\\');
			ext = strrchr(fileName, '.');
			ext++;

			bool ok = true;
			if (FileUtils::Exist(fileName)) {
				sprintf(tmp, "Overwrite existing file ?\n%s", fileName);
				ok = (GLMessageBox::Display(tmp, "Question", GLDLG_OK | GLDLG_CANCEL, GLDLG_ICONWARNING) == GLDLG_OK);
			}
			
			if (ok) {
				f.open(fileName, std::ofstream::out | std::ios::trunc);
				if (!f) {
					char tmp[256];
					sprintf(tmp, "Cannot open file for writing %s", fileName);
					throw Error(tmp);
				}

				int numCol = historyList->GetNbColumn();
				int numWidth = 20;

				int yScaleMode = yScaleCombo->GetSelectedIndex();
				switch (yScaleMode) {
				case 0: { //MC Hits
					f <<"Table - MC Hits\n";

					break;
				}
				case 1: { //Equiv Hits
					f<<"Table - Equiv. Hits\n";

					break;
				}
				case 2: {//Pressure
					f<<"Table - Pressure (mbar)\n";

					break;
				}
				case 3: {//Particle density
					f<<"Table - Density (1/m3)\n";
					break;
				}
				case 4: {//Imp.rate
					f<<"Table - Imp.rate (1/s/m2)\n";
					break;
				}
				}

				for (int row = 0;row < historyList->GetNbRow();row++)
				{
					if (row == 0) {
						for (int col = 0; col < numCol;col++) {
							if (col == 0) {
								f << left << setw(10) << historyList->getColumnLabel(col) <<'\t';
							}
							else{
								f << left << setw(numWidth) << historyList->getColumnLabel(col);
								if (col == numCol - 1) { f << '\n'; }
								else { f << "\t"; }
							}
						}
					}
					for (int col = 0; col < numCol;col++) {
						if (col == 0) {
							f << left << setw(10) << historyList->GetValueAt(col, row) << '\t';
						}
						else {
							f << left << setw(numWidth) << historyList->GetValueAt(col, row);
							if (col == numCol - 1) { f << '\n'; }
							else { f << "\t"; }
						}
					}
				}

				f.close();
			}

		}
		catch (Error &e) {
			char errMsg[512];
			sprintf(errMsg, "%s\nFile:%s", e.GetMsg(), fn->fullName);
			GLMessageBox::Display(errMsg, "Error", GLDLG_OK, GLDLG_ICONERROR);
		}
	}
}

void PressureEvolution::updateList(BYTE *buffer, Facet *f, int colnum) {
	char tmp[256];
	int yScaleMode = yScaleCombo->GetSelectedIndex();
	GlobalHitBuffer *gHits = (GlobalHitBuffer *)buffer;
	double nbDes = (double)gHits->globalHits.nbDesorbed;
	double scaleY;
	
	switch (yScaleMode) {
	case 0: { //MC Hits
		panelTable->SetTitle("Table - MC Hits");
		for (size_t m = 1; m <= Min(worker->moments.size(), (size_t)10000); m++) { //max 10000 points
			FacetHitBuffer* facetHits = (FacetHitBuffer*)(buffer + f->sh.hitOffset + m * sizeof(FacetHitBuffer));


			sprintf(tmp, "%g", (double)facetHits->nbMCHit);
			historyList->SetValueAt(colnum + 1, m - 1, tmp);
		}
		break;
	}
	case 1: { //Equiv Hits
		panelTable->SetTitle("Table - Equiv. Hits");
		for (size_t m = 1; m <= Min(worker->moments.size(), (size_t)10000); m++) { //max 10000 points
			FacetHitBuffer* facetHits = (FacetHitBuffer*)(buffer + f->sh.hitOffset + m * sizeof(FacetHitBuffer));

			sprintf(tmp, "%g", facetHits->nbHitEquiv);
			historyList->SetValueAt(colnum + 1, m - 1, tmp);
		}
		break;
	}
	case 2: {//Pressure
		panelTable->SetTitle("Table - Pressure (mbar)");
		scaleY = 1.0 / nbDes / (f->sh.area * 1E-4)* worker->wp.gasMass / 1000 / 6E23 * 0.0100; //0.01: Pa->mbar
		scaleY *= worker->wp.totalDesorbedMolecules / worker->wp.timeWindowSize;
		if (f->sh.is2sided) scaleY *= 0.5;
		for (size_t m = 1; m <= Min(worker->moments.size(), (size_t)10000); m++) { //max 10000 points
			FacetHitBuffer* facetHits = (FacetHitBuffer*)(buffer + f->sh.hitOffset + m * sizeof(FacetHitBuffer));


			sprintf(tmp, "%g", facetHits->sum_v_ort*scaleY);
			historyList->SetValueAt(colnum + 1, m - 1, tmp);
		}
		break;
	}
	case 3: {//Particle density
		panelTable->SetTitle("Table - Density (1/m3)");
		scaleY = 1.0 / nbDes / (f->GetArea() * 1E-4);
		scaleY *= worker->wp.totalDesorbedMolecules / worker->wp.timeWindowSize;
		scaleY *= f->DensityCorrection();
		for (size_t m = 1; m <= Min(worker->moments.size(), (size_t)10000); m++) { //max 10000 points
			FacetHitBuffer* facetHits = (FacetHitBuffer*)(buffer + f->sh.hitOffset + m * sizeof(FacetHitBuffer));


			sprintf(tmp, "%g", facetHits->sum_1_per_ort_velocity*scaleY);
			historyList->SetValueAt(colnum + 1, m - 1, tmp);
		}
		break;
	}
	case 4: {//Imp.rate
		panelTable->SetTitle("Table - Imp.rate (1/s/m2)");
		scaleY = 1.0 / nbDes / (f->GetArea() * 1E-4);
		scaleY *= worker->wp.totalDesorbedMolecules / worker->wp.timeWindowSize;
		for (size_t m = 1; m <= Min(worker->moments.size(), (size_t)10000); m++) { //max 10000 points
			FacetHitBuffer* facetHits = (FacetHitBuffer*)(buffer + f->sh.hitOffset + m * sizeof(FacetHitBuffer));


			sprintf(tmp, "%g", facetHits->nbHitEquiv*scaleY);
			historyList->SetValueAt(colnum + 1, m - 1, tmp);
		}
		break;
	}
	}
}

void PressureEvolution::refreshTable() {

	// Lock during update
	BYTE *buffer = worker->GetHits();

	if (!buffer) return;

	Geometry *geom = worker->GetGeometry();


	size_t numMomentsNow = Min(worker->moments.size(), (size_t)10000);

	char tmp[256];
	if (onlyViews) {
		historyList->SetSize(views.size() + 1, numMomentsNow);
		historyList->SetColumnLabels();

		for (size_t m = 1; m <= Min(worker->moments.size(), (size_t)10000); m++) { //max 10000 points
			sprintf(tmp, "%g", (double)worker->moments[m - 1]);
			historyList->SetValueAt(0, m - 1, tmp);
		}

		int viewnum = 0;
		for (auto& v : views) {

			if (v->userData1 >= 0 && v->userData1 < geom->GetNbFacet()) {
				Facet *f = geom->GetFacet(v->userData1);
				updateList(buffer, f, viewnum);
				historyList->SetColumnLabel(viewnum+1, v->userData1+1);
			}
			viewnum += 1;
		}
	}
	else {
		nb_Facets = geom->GetNbFacet();
		historyList->SetSize(nb_Facets + 1, numMomentsNow);
		historyList->SetColumnLabels();
		//historyList->SetColumnWidthForAll((this->width - 40) / (nb_Facets + 1));


		for (size_t m = 1; m <= Min(worker->moments.size(), (size_t)10000); m++) { //max 10000 points
			sprintf(tmp, "%g", (double)worker->moments[m - 1]);
			historyList->SetValueAt(0, m - 1, tmp);
			for (auto& v : views) {

			}
		}

		for (int i = 0;i < nb_Facets;i++) {
			Facet *f = geom->GetFacet(i);
			updateList(buffer, f, i);

		}
	}
	historyList->SetColumnWidthForAll(90);
	historyList->SetColumnWidth(0, 60);
	worker->ReleaseHits();
}

void PressureEvolution::refreshChart() {
	//refreshes chart values
	refreshTable();

	// Lock during update
	BYTE *buffer = worker->GetHits();
	int yScaleMode = yScaleCombo->GetSelectedIndex();
	if (!buffer) return;

	Geometry *geom = worker->GetGeometry();
	GlobalHitBuffer *gHits = (GlobalHitBuffer *)buffer;
	double nbDes = (double)gHits->globalHits.nbDesorbed;
	double scaleY;
	size_t facetHitsSize = (1 + worker->moments.size()) * sizeof(FacetHitBuffer);

	for (auto& v : views) {

		if (v->userData1 >= 0 && v->userData1 < geom->GetNbFacet()) {
			Facet *f = geom->GetFacet(v->userData1);
			v->Reset();
			
			switch (yScaleMode) {
			case 0: { //MC Hits
				for (size_t m = 1; m <= Min(worker->moments.size(), (size_t)10000); m++) { //max 10000 points
					FacetHitBuffer* facetHits = (FacetHitBuffer*)(buffer + f->sh.hitOffset + m * sizeof(FacetHitBuffer));
					v->Add(worker->moments[m - 1], (double)facetHits->nbMCHit, false);
				}
				break;
			}
			case 1: { //Equiv Hits
				for (size_t m = 1; m <= Min(worker->moments.size(), (size_t)10000); m++) { //max 10000 points
					FacetHitBuffer* facetHits = (FacetHitBuffer*)(buffer + f->sh.hitOffset + m * sizeof(FacetHitBuffer));
					v->Add(worker->moments[m - 1], facetHits->nbHitEquiv, false);
				}
				break;
			}
			case 2: {//Pressure
				scaleY = 1.0 / nbDes / (f->sh.area * 1E-4)* worker->wp.gasMass / 1000 / 6E23 * 0.0100; //0.01: Pa->mbar
				scaleY *= worker->wp.totalDesorbedMolecules / worker->wp.timeWindowSize;
				if (f->sh.is2sided) scaleY *= 0.5;
				for (size_t m = 1; m <= Min(worker->moments.size(), (size_t)10000); m++) { //max 10000 points
					FacetHitBuffer* facetHits = (FacetHitBuffer*)(buffer + f->sh.hitOffset + m * sizeof(FacetHitBuffer));
					v->Add(worker->moments[m - 1], facetHits->sum_v_ort*scaleY, false);
				}
				break;
			}
			case 3: {//Particle density
				scaleY = 1.0 / nbDes / (f->GetArea() * 1E-4);
				scaleY *= worker->wp.totalDesorbedMolecules / worker->wp.timeWindowSize;
				scaleY *= f->DensityCorrection();
				for (size_t m = 1; m <= Min(worker->moments.size(), (size_t)10000); m++) { //max 10000 points
					FacetHitBuffer* facetHits = (FacetHitBuffer*)(buffer + f->sh.hitOffset + m * sizeof(FacetHitBuffer));
					v->Add(worker->moments[m - 1], facetHits->sum_1_per_ort_velocity*scaleY, false);
				}
				break;
			}
			case 4: {//Imp.rate
				scaleY = 1.0 / nbDes / (f->GetArea() * 1E-4);
				scaleY *= worker->wp.totalDesorbedMolecules / worker->wp.timeWindowSize;
				for (size_t m = 1; m <= Min(worker->moments.size(), (size_t)10000); m++) { //max 10000 points
					FacetHitBuffer* facetHits = (FacetHitBuffer*)(buffer + f->sh.hitOffset + m * sizeof(FacetHitBuffer));
					v->Add(worker->moments[m - 1], facetHits->nbHitEquiv*scaleY, false);
				}
				break;
			}
			}
		}
		v->CommitChange();
	}

	worker->ReleaseHits();
}

void PressureEvolution::addView(size_t facetId) {

	Geometry *geom = worker->GetGeometry();

	// Check that view is not already added
	{
		bool found = false;

		for (auto v = views.begin(); v != views.end() && !found;v++) {
			found = ((*v)->userData1 == facetId);
		}
		if (found) {
			GLMessageBox::Display("Facet already on chart", "Error", GLDLG_OK, GLDLG_ICONERROR);
			return;
		}
		if (worker->moments.size() > 10000) {
			GLMessageBox::Display("Only the first 10000 moments will be plotted", "Error", GLDLG_OK, GLDLG_ICONWARNING);
		}
	}

	auto v = new GLDataView();
	std::ostringstream tmp;
	tmp << "Facet " << facetId + 1;
	v->SetName(tmp.str().c_str());
	v->SetViewType(TYPE_BAR);
	v->SetMarker(MARKER_DOT);
	v->SetColor(colors[views.size() % colors.size()]);
	v->SetMarkerColor(colors[views.size() % colors.size()]);
	v->userData1 = (int)facetId;
	views.push_back(v);
	chart->GetY1Axis()->AddDataView(v);
	Refresh();
}

void PressureEvolution::remView(size_t viewId) {

	chart->GetY1Axis()->RemoveDataView(views[viewId]);
	views.erase(views.begin()+viewId);

}

void PressureEvolution::Reset() {

	chart->GetY1Axis()->ClearDataView();
	for (auto v : views)
		delete v;
	views.clear();
	Refresh();
}

void PressureEvolution::ProcessMessage(GLComponent *src, int message) {
	Geometry *geom = worker->GetGeometry();
	switch (message) {
	case MSG_BUTTON:
		if (src == selButton) {
			int idx = profCombo->GetSelectedIndex();
			if (idx >= 0) {
				size_t facetId = profCombo->GetUserValueAt(idx);
				if (facetId >= 0 && facetId < geom->GetNbFacet()) {
					geom->UnselectAll();
					geom->GetFacet(facetId)->selected = true;
					geom->UpdateSelection();
					mApp->UpdateFacetParams(true);
					mApp->facetList->SetSelectedRow((int)facetId);
					mApp->facetList->ScrollToVisible(facetId, 1, true);
				}
			}
		}
		else if (src == addButton) {
			auto selFacets = geom->GetSelectedFacets();
			if (selFacets.size() != 1) {
				GLMessageBox::Display("Select exactly one facet", "Add selected facet to chart", { "Sorry!" }, GLDLG_ICONERROR);
				return;
			}
			else {
				addView(selFacets[0]); //Includes chart refresh
			}
		}
		else if (src == removeButton) {
			int idx = profCombo->GetSelectedIndex();
			if (idx >= 0) remView(idx);
			Refresh();
		}
		else if (src == removeAllButton) {
			Reset();
		}
		else if (src == refreshButton) {
			Refresh();
		}
		else if (src == exportButton) {
			exportList();
		}
		else if (src == onlyViewsButton) {
			if (onlyViews == true) {
				onlyViews = false;
				onlyViewsButton->SetText("Show only selected facets");
			}
			else {
				onlyViews = true;
				onlyViewsButton->SetText("Show all facets");
			}
			refreshTable();
		}
		break;
	case MSG_COMBO:
		if (src == yScaleCombo) {
			refreshChart();
		}
		break;
	case MSG_TOGGLE:
		if (src == logXToggle) {
			chart->GetXAxis()->SetScale(logXToggle->GetState());
		}
		else if (src == logYToggle) {
			chart->GetY1Axis()->SetScale(logYToggle->GetState());
		}
		break;
	}

	GLWindow::ProcessMessage(src, message);

}

