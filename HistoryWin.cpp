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

//#include "GLApp/GLLabel.h"
#include "HistoryWin.h"
#include "GLApp/GLToolkit.h"
#include "GLApp/GLToggle.h"
#include "GLApp/GLMessageBox.h"
#include "GLApp/GLList.h"
#include "GLApp/GLInputBox.h"
#include "GLApp/GLFileBox.h"
#include "Facet_shared.h"
#include "Geometry_shared.h"
#include "GLApp/MathTools.h"
#include "GLApp/GLButton.h"
#include "GLApp/GLTextField.h"
#include "GLApp/GLLabel.h"
#include "GLApp/GLToggle.h"
#include "GLApp/GLTitledPanel.h"
#include "Buffer_shared.h"
#include "Interface.h"
#include <fstream>
#ifdef MOLFLOW
#include "MolFlow.h"
#endif

#ifdef SYNRAD
#include "SynRad.h"
#endif

#include <direct.h>

extern GLApplication *theApp;

#ifdef MOLFLOW
extern MolFlow *mApp;
#endif

#ifdef SYNRAD
extern SynRad*mApp;
#endif

HistoryWin::HistoryWin(Worker *w) :GLWindow() {
	// worker
	worker = w;

	// window size
	int wD = 580;
	int hD = 525;

	// UI Setup
	SetTitle("History");
	SetIconfiable(true);


	geom = worker->GetGeometry();
	nb_Facets = geom->GetNbFacet();
	GLTitledPanel *panel3 = new GLTitledPanel("History of covering");
	panel3->SetBounds(5, 10, wD - 10, hD - 35);
	Add(panel3);

	historyList = new GLList(0);
	historyList->SetHScrollVisible(true);
	historyList->SetSize(nb_Facets+1, 1);
	historyList->SetColumnLabelVisible(true);
	historyList->SetSelectionMode(1);
	panel3->Add(historyList);
	panel3->SetCompBounds(historyList,10, 50, wD - 30, hD - 100);

	histBtn = new GLButton(0, "Update Covering");
	histBtn->SetVisible(true);
	panel3->Add(histBtn);
	panel3->SetCompBounds(histBtn, 10, 20, 100, 18);

	delBtn = new GLButton(0, "Delete Entry");
	delBtn->SetVisible(true);
	panel3->Add(delBtn);
	panel3->SetCompBounds(delBtn, 120, 20, 100, 18);

	imBtn = new GLButton(0, "Import");
	imBtn->SetVisible(true);
	panel3->Add(imBtn);
	panel3->SetCompBounds(imBtn, wD-125-110, 20, 100, 18);

	exBtn = new GLButton(0, "Export");
	exBtn->SetVisible(true);
	panel3->Add(exBtn);
	panel3->SetCompBounds(exBtn, wD-125, 20, 100, 18);

	// variable initialization
	l_hist = 0;
	pointintime_list = std::vector< std::pair<float, std::vector<boost::multiprecision::uint128_t>> >();
	selectedRows = std::vector<size_t>();
	// set initial list
	if(nb_Facets!=0)
		UpdateList();
	else {
		historyList->SetValueAt(0, 0, "No Geometry");
		historyList->SetColumnLabelVisible(false);
		historyList->AutoSizeColumn();
		histBtn->SetVisible(false);
		delBtn->SetVisible(false);
		exBtn->SetVisible(false);
		imBtn->SetVisible(false);
	}

	
	// Center dialog, position for window?
	int wS, hS;
	GLToolkit::GetScreenSize(&wS, &hS);
	int xD = (wS - wD) / 2;
	int yD = (hS - hD) / 2;
	SetBounds(xD, yD, wD, hD);

	RestoreDeviceObjects();

	//for (size_t i = 0; i < MAX_PROCESS; i++) lastCPUTime[i] = -1.0f;
	//memset(lastCPULoad, 0, MAX_PROCESS*sizeof(float));
}

void HistoryWin::UpdateList() {
	float time = 0;
	if (worker->isRunning)
		time = float(1000.0)*(worker->simuTime + (mApp->m_fTime - worker->startTime));
	else
		time = worker->simuTime * (float)1000.0;

	std::vector<boost::multiprecision::uint128_t> covering;
	covering = std::vector<boost::multiprecision::uint128_t>();


	for (int i = 0;i < nb_Facets;i++) {
		Facet *f = geom->GetFacet(i);
		boost::multiprecision::uint128_t cov = f->facetHitCache.covering;
		covering.push_back(cov);
		//sprintf(tmp, "%g", cov);
		//historyList->SetValueAt(i+1,l_hist, tmp);
	}
	//sprintf(tmp, "%g", time);
	//historyList->SetValueAt(0, l_hist, tmp);
	pointintime_list.push_back(std::make_pair(time, covering));

	//historyList->AutoSizeColumn();
	l_hist += 1;
	
	UpdateUI();
}

void HistoryWin::UpdateUI() {
	char tmp[256];

	historyList->SetSize(nb_Facets + 1, l_hist);
	historyList->SetColumnLabels();
	for (int i = 0;i < l_hist;i++) {
		sprintf(tmp, "%g", pointintime_list[i].first);
		historyList->SetValueAt(0, i, tmp);
		for (int j = 0; j < nb_Facets; j++) {
			strcpy(tmp, (pointintime_list[i].second[j].str()).c_str());
			historyList->SetValueAt(j+1,i, tmp);
		}
	}

	historyList->AutoSizeColumn();
	historyList->SetAllColumnAlign(1);
	historyList->SetColumnAlign(0, 0);
}

void HistoryWin::deleteRows() {
	selectedRows = historyList->GetSelectedRows();

	if (selectedRows.size() != 0) {
		std::sort(selectedRows.begin(), selectedRows.end(), std::less<size_t>());

		if (selectedRows.front() == 0)
			selectedRows.erase(selectedRows.begin());
		if (selectedRows.size() != 0) {
			if (selectedRows.back() == l_hist - 1)
				selectedRows.pop_back();
		}
	}

	if (selectedRows.size() != 0) {
		char tmp[512];
		sprintf(tmp, "Deleting %zi row(s)?\nNote: First and last row cannot be deleted\n", selectedRows.size());
		bool ok = (GLMessageBox::Display(tmp, "Question", GLDLG_OK | GLDLG_CANCEL, GLDLG_ICONWARNING) == GLDLG_OK);

		if (!ok)
			return;

		for (size_t i = 0; i < selectedRows.size(); i++) {

			size_t row = selectedRows[i]-i;
			pointintime_list.erase(row + pointintime_list.begin());
			l_hist -= 1;
		}
		UpdateUI();
		historyList->ClearSelection();
	}
	else {
		GLMessageBox::Display("No rows selected \nNote: First and last row cannot be deleted\n", "Error", GLDLG_OK | GLDLG_CANCEL, GLDLG_ICONWARNING);
	}
}

void HistoryWin::ProcessMessage(GLComponent *src, int message) {
	switch (message) {
	case MSG_BUTTON:
		if (src == histBtn) {
			UpdateList();
		}

		else if (src == delBtn) {
			deleteRows();
		}

		else if (src == exBtn) {
			exportUI();
		}

		else if (src == imBtn) {
			importUI();
		}
	}
	GLWindow::ProcessMessage(src, message);
}

void HistoryWin::exportList(char *fileName) {
	char tmp[512];

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

		for (int i = 0;i < pointintime_list.size();i++)
		{
			f << pointintime_list[i].first;

			for (int j = 0; j < pointintime_list[i].second.size();j++)
			{
				f << '\t' << pointintime_list[i].second[j];
			}

			f << '\n';

		}

		f.close();
	}
}

void HistoryWin::exportUI() {
	FILENAME *fn = GLFileBox::SaveFile(mApp->currentDir, NULL, "Save File", "All files\0*.*\0", 0);
	
	if (fn) {
		try {
			exportList(fn->fullName);
		}
		catch (Error &e) {
			char errMsg[512];
			sprintf(errMsg, "%s\nFile:%s", e.GetMsg(), fn->fullName);
			GLMessageBox::Display(errMsg, "Error", GLDLG_OK, GLDLG_ICONERROR);
		}
	}
}

void HistoryWin::importUI() {

	//GLMessageBox::Display("Please choose 'All files' when selecting the import buffer in the next dialog. Buffer files do not have a file ending. They are not included in 'All Molflow supported files'.", "Info", GLDLG_OK, GLDLG_ICONINFO);
	char fullName[512];
	char shortName[512];
	strcpy(fullName, "");


	FILENAME *fn = GLFileBox::OpenFile(mApp->currentDir, NULL, "Open file", "All files\0*.*\0", 0);
	if (fn)
		strcpy(fullName, fn->fullName);

	if (strlen(fullName) == 0) {
		return;
	}
	char *lPart = strrchr(fullName, '\\');
	if (lPart) strcpy(shortName, lPart + 1);
	else strcpy(shortName, fullName);
	try {
		importList(fullName);
	}
	catch (Error &e) {
		char errMsg[512];
		sprintf(errMsg, "%s\nFile:%s", e.GetMsg(), shortName);
		GLMessageBox::Display(errMsg, "Error", GLDLG_OK, GLDLG_ICONERROR);
	}
}

void HistoryWin::importList(char *fileName) {
	char CWD[MAX_PATH];
	_getcwd(CWD, MAX_PATH);
	std::string ext = FileUtils::GetExtension(fileName);
	if (ext != "txt")
		throw Error("void HistoryWin::importList(): Wrong file extension.");
	char tmp[512];
	sprintf(tmp, "Overwrite complete History?\n");
	bool ok = (GLMessageBox::Display(tmp, "Question", GLDLG_OK | GLDLG_CANCEL, GLDLG_ICONWARNING) == GLDLG_OK);

	if (!ok)
		return;

	std::ifstream input(fileName, std::ifstream::in);
	if (input) {
		worker->Stop_Public();

		
		// read data
		pointintime_list.clear();
		l_hist = 0;
		//std::string read = "/home/van/history"+std::to_string(num)+".txt";
		std::string line;

		while (std::getline(input, line)) {
			std::vector<boost::multiprecision::uint128_t> currentstep;
			currentstep = std::vector<boost::multiprecision::uint128_t>();

			boost::multiprecision::uint128_t covering;
			float time;
			std::istringstream is(line);

			is >> time;
			while (!is.eof()) {
				is >> covering;
				//std::cout <<time <<'\t';
				currentstep.push_back(covering);

			}
			pointintime_list.push_back(std::make_pair(time, currentstep));

			l_hist += 1;

		}
		//Update(mApp->m_fTime);
		UpdateCoveringfromList();
		UpdateUI();
	}
	else{
		char tmp[256];
		sprintf(tmp, "Cannot open file for writing %s", fileName);
		throw Error(tmp);
	}
	input.close();
}

void HistoryWin::UpdateCoveringfromList() {
	for (int i = 0;i < nb_Facets;i++) {
		Facet *f = geom->GetFacet(i);
		f->facetHitCache.covering = pointintime_list.back().second[i];
	}
	worker->simuTime = pointintime_list.back().first / (float)1000.0;
}