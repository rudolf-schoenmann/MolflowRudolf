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
#ifndef _HISTORYWINH_
#define _HISTORYWINH_

#include "GLApp/GLWindow.h"
#include <vector>
#include <boost/multiprecision/cpp_int.hpp>

class GLButton;
class GLTextField;
class GLLabel;
class GLToggle;
class GLTitledPanel;

class Worker;
class Geometry;
class GLList;


class HistoryWin : public GLWindow {

public:
	// Construction
	HistoryWin(Worker *w);
	std::vector< std::pair<float, std::vector<boost::multiprecision::uint128_t>> > pointintime_list;
	size_t l_hist;

	// Implementation
	void ProcessMessage(GLComponent *src, int message);
	void UpdateList();
	void UpdateUI();
	void deleteRows();

private:
	std::vector<size_t> selectedRows;
	//void RestartProc();
	Worker      *worker;
	Geometry *geom;
	size_t nb_Facets;
	GLList      *historyList;
	GLButton *histBtn;
	GLButton *delBtn;
	GLButton *imBtn;
	GLButton *exBtn;

	void exportUI();
	void exportList(char *fileName);
	void importList(char *fileName);
	void importUI();
	void UpdateCoveringfromList();

};

#endif /* _HISTORYWINH_ */
