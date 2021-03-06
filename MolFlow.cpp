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
#include <math.h>
//#include <malloc.h>
#include "MolFlow.h"
#include "Facet_shared.h"
#include "MolflowGeometry.h"
#include "File.h"
#include "GLApp/GLMessageBox.h"
#include "GLApp/GLInputBox.h"
#include "GLApp/GLFileBox.h"
#include "GLApp/GLToolkit.h"
#include "GLApp/GLWindowManager.h"
#include "GLApp/MathTools.h"
#include "GLApp/GLMenubar.h"
#include "GLApp/GLButton.h"
#include "GLApp/GLLabel.h"
#include "GLApp/GLCombo.h"
#include "GLApp/GLTextField.h"
#include "RecoveryDialog.h"
#include "direct.h"
#include <vector>
#include <string>
#include <io.h>
#include <thread>
#include <numeric> //std::iota

#include "Interface.h"
//#include "AppUpdater.h"
#include "Worker.h"
#include "ImportDesorption.h"
#include "TimeSettings.h"
#include "Movement.h"
#include "FacetAdvParams.h"
#include "FacetDetails.h"
#include "Viewer3DSettings.h"
#include "TextureScaling.h"
#include "GlobalSettings.h"
#include "HistoryWin.h"
#include "ProfilePlotter.h"
#include "PressureEvolution.h"
#include "TimewisePlotter.h"
#include "TexturePlotter.h"
#include "OutgassingMap.h"
#include "MomentsEditor.h"
#include "FacetCoordinates.h"
#include "VertexCoordinates.h"
#include "ParameterEditor.h"
#include "SmartSelection.h"
#include "FormulaEditor.h"
#include "ParticleLogger.h"
#include "HistogramSettings.h"
#include "HistogramPlotter.h"
#include <iostream>

//Hard-coded identifiers, update these on new release
//---------------------------------------------------
std::string appName = "ContaminationFlow";
int appVersionId = 0; //Compared with available updates. Recompile Interface.cpp if changed
std::string appVersionName = "0";
//---------------------------------------------------

static const char *fileLFilters = "All MolFlow supported files\0*.txt;*.xml;*.zip;*.geo;*.geo7z;*.syn;*.syn7z;*.str;*.stl;*.ase\0"
"All files\0*.*\0";
static const char *fileInsFilters = "All insertable geometries\0*.txt;*.xml;*.zip;*.geo;*.geo7z;*.syn;*.syn7z;*.stl\0"
"All files\0*.*\0";
const char *fileSFilters = "MolFlow saveable files\0*.xml;*.zip;*.geo;*.geo7z;*.txt\0All files\0*.*\0";
static const char *fileDesFilters = "Desorption files\0*.des\0All files\0*.*\0";

int cSize = 4;
int offset_width = 0;
int   cWidth[] = { 30+offset_width, 56, 50, 50 };
char *cName[] = { "#", "Hits", "Des", "Abs" };

#ifdef _DEBUG
std::string appTitle = "ContaminationFlow debug version (Compiled " __DATE__ " " __TIME__ ")";
#else
std::string appTitle = "ContaminationFlow " + appVersionName + " (" __DATE__ ")";
#endif

std::vector<string> formulaPrefixes = { "A","D","H","MCH","P","DEN","Z","V","T","AR","a","d","h","mch","p","den","z","v","t","ar","," };
std::string formulaSyntax =
R"(MC Variables: An (Absorption on facet n), Dn (Desorption on facet n), Hn (Hit on facet n)
Pn (Pressure [mbar] on facet n), DENn (Density [1/m3] on facet n)
Zn (Imp. rate on facet n), Vn (avg. speed [m/s] on facet n), Tn (temp[K] of facet n)
SUMABS (total absorbed), SUMDES (total desorbed), SUMHIT (total hit)

SUM(H,3,8)    calculates the sum of hits on facets 3,4,... ...7,8. (Works with H,A,D,AR).
AVG(P,3,8)    calculates the average pressure (area-wise) on facets 3 to 8 (Works with P,D,Z)
SUM(H,S3)    calculates the sum of hits on selection group 3 (works with H,A,D,AR)
AVG(DEN,S2) calculates the average (area-wise) on facets belonging to sel. group 2
SUM(H,SEL)    calculates the sum of hits on the current selection. (Works with H,A,D,AR)
AVG(Z,SEL)    calculates the average impingement rate on the current selection
For the last two, might need to manually refresh formulas after you change the selection.

Area variables: ARn (Area of facet n), DESAR (total desorption area), ABSAR (total absorption area)

Final (constant) outgassing rate [mbar*l/s]: QCONST
Final (constant) outgassing rate [molecules/s]: QCONST_N
Total desorbed molecules until last moment: [molecules]: NTOT
Gas mass [g/mol]: sh.gasMass

Mean Pumping Path: MPP (average path of molecules in the system before absorption)
Mean Free Path:      MFP (average path of molecules between two wall hits)

Math functions: sin(), cos(), tan(), sinh(), cosh(), tanh(), asin(), acos(),
                     atan(), exp(), ln(), pow(x,y), log2(), log10(), inv(), sqrt(), abs()

Constants:  Kb (Boltzmann's constant), R (Gas constant), Na (Avogadro's number), PI
)";
int formulaSyntaxHeight = 380;

MolFlow *mApp;

//Menu elements, Molflow specific:
//#define MENU_FILE_IMPORTDES_DES 140
#define MENU_FILE_IMPORTDES_SYN 141

#define MENU_FILE_EXPORTTEXTURE_AREA 151
#define MENU_FILE_EXPORTTEXTURE_MCHITS 152
#define MENU_FILE_EXPORTTEXTURE_IMPINGEMENT 153
#define MENU_FILE_EXPORTTEXTURE_PART_DENSITY 154
#define MENU_FILE_EXPORTTEXTURE_GAS_DENSITY 155
#define MENU_FILE_EXPORTTEXTURE_PRESSURE 156
#define MENU_FILE_EXPORTTEXTURE_AVG_V 157
#define MENU_FILE_EXPORTTEXTURE_V_VECTOR 158
#define MENU_FILE_EXPORTTEXTURE_N_VECTORS 159

#define MENU_FILE_EXPORTTEXTURE_AREA_COORD 171
#define MENU_FILE_EXPORTTEXTURE_MCHITS_COORD  172
#define MENU_FILE_EXPORTTEXTURE_IMPINGEMENT_COORD  173
#define MENU_FILE_EXPORTTEXTURE_PART_DENSITY_COORD  174
#define MENU_FILE_EXPORTTEXTURE_GAS_DENSITY_COORD  175
#define MENU_FILE_EXPORTTEXTURE_PRESSURE_COORD  176
#define MENU_FILE_EXPORTTEXTURE_AVG_V_COORD  177
#define MENU_FILE_EXPORTTEXTURE_V_VECTOR_COORD  178
#define MENU_FILE_EXPORTTEXTURE_N_VECTORS_COORD  179
//#define MENU_FILE_EXPORTBUFFER 180
#define MENU_FILE_EXPORTBUFFER_LOAD 181
#define MENU_FILE_EXPORTBUFFER_HIT 182
//#define MENU_FILE_IMPORTBUFFER_HIT 183

#define MENU_TOOLS_MOVINGPARTS 410

#define MENU_FACET_MESH        360
#define MENU_SELECT_HASDESFILE 361
#define MENU_FACET_OUTGASSINGMAP 362

#define MENU_TIME_SETTINGS          900
#define MENU_TIMEWISE_PLOTTER       901
#define MENU_TIME_PRESSUREEVOLUTION 902
#define MENU_TIME_MOMENTS_EDITOR    903
#define MENU_TIME_PARAMETER_EDITOR  904

// Name: WinMain()
// Desc: Entry point to the program. Initializes everything, and goes into a
//       message-processing loop. Idle time is used to render the scene.

INT WINAPI WinMain(HINSTANCE hInst, HINSTANCE, LPSTR, INT)
{
	MolFlow *mApp = new MolFlow();

	if (!mApp->Create(1024, 800, false)) {
		char *logs = GLToolkit::GetLogs();
#ifdef WIN
		if (logs) MessageBox(NULL, logs, "Molflow [Fatal error]", MB_OK);
#else
		if (logs) {
			printf("Molflow [Fatal error]\n");
			printf(logs);
		}
#endif
		SAFE_FREE(logs);
		delete mApp;
		return -1;
	}
	try {
		mApp->Run();
	}
	catch (Error &e) {
		mApp->CrashHandler(&e);
	}
	delete mApp;
	return 0;
}

// Name: MolFlow()
// Desc: Application constructor. Sets default attributes for the app.

MolFlow::MolFlow()
{
	mApp = this; //to refer to the app as extern variable

	//Different Molflow implementation:
	facetAdvParams = NULL;
	facetDetails = NULL;
	viewer3DSettings = NULL;
	textureScaling = NULL;
	formulaEditor = NULL;
	globalSettings = NULL;
	history = NULL;
	profilePlotter = NULL;
	texturePlotter = NULL;

	//Molflow only:
	movement = NULL;
	timewisePlotter = NULL;
	pressureEvolution = NULL;
	outgassingMap = NULL;
	momentsEditor = NULL;
	parameterEditor = NULL;
	importDesorption = NULL;
	timeSettings = NULL;
	
}

// Name: OneTimeSceneInit()
// Desc: Called during initial app startup, this function performs all the
//       permanent initialization.

void MolFlow::LoadParameterCatalog()
{
	//Parameters from catalog
	//Find in parameter_catalog folder
	intptr_t file;
	_finddata_t filedata;
	file = _findfirst("parameter_catalog\\*.csv", &filedata);
	if (file != -1)
	{
		do
		{
			Parameter newParam;
			newParam.fromCatalog = true;
			newParam.isLogLog = true; //For now
			std::stringstream parName; parName /*<< "["*/ << filedata.name << " [catalog]";
			newParam.name = parName.str();

			std::vector<std::vector<string>> table;
			try {
				std::stringstream pathAndName; pathAndName << "parameter_catalog\\" << filedata.name;
				FileReader *f = new FileReader(pathAndName.str());
				table = worker.ImportCSV_string(f);
				SAFE_DELETE(f);
			}
			catch (Error &e) {
				char errMsg[512];
				sprintf(errMsg, "Failed to load CSV file.\n%s", e.GetMsg());
				GLMessageBox::Display(errMsg, "Error", GLDLG_OK, GLDLG_ICONERROR);
				continue;
			}
			//Parse
			for (size_t i = 0; i < table.size(); i++) {
				std::vector<std::string> row = table[i];
				if (row.size() != 2) {
					std::stringstream errMsg;
					errMsg << filedata.name << " Row " << i + 1 << "has " << row.size() << " values instead of 2.";
					GLMessageBox::Display(errMsg.str().c_str(), "Error", GLDLG_OK, GLDLG_ICONERROR);
					break;
				}
				else {
					double valueX, valueY;
					try {
						valueX = ::atof(row[0].c_str());
					}
					catch (std::exception err) {
						char tmp[256];
						sprintf(tmp, "Can't parse value \"%s\" in row %zd, first column:\n%s", row[0].c_str(), i + 1, err.what());
						GLMessageBox::Display(tmp, "Invalid parameter definition", GLDLG_OK, GLDLG_ICONWARNING);
						break;
					}
					try {
						valueY = ::atof(row[1].c_str());
					}
					catch (std::exception err) {
						char tmp[256];
						sprintf(tmp, "Can't parse value \"%s\" in row %zd, second column:\n%s", row[1].c_str(), i + 1, err.what());
						GLMessageBox::Display(tmp, "Invalid parameter definition", GLDLG_OK, GLDLG_ICONWARNING);
						break;
					}
					newParam.AddPair(valueX, valueY, true); //insert in correct position
				}
			}

			worker.parameters.push_back(newParam);
		} while (_findnext(file, &filedata) == 0);
	}
	_findclose(file);
}

int MolFlow::OneTimeSceneInit()
{
	/*
	//Enable memory check at EVERY malloc/free operation:
	int tmpFlag = _CrtSetDbgFlag( _CRTDBG_REPORT_FLAG );
	tmpFlag |= _CRTDBG_LEAK_CHECK_DF;
	_CrtSetDbgFlag( tmpFlag );
	*/

	OneTimeSceneInit_shared_pre();

	menu->GetSubMenu("File")->Add("Export selected textures");
	menu->GetSubMenu("File")->GetSubMenu("Export selected textures")->Add("Facet by facet");
	menu->GetSubMenu("File")->GetSubMenu("Export selected textures")->GetSubMenu("Facet by facet")->Add("Cell Area (cm\262)", MENU_FILE_EXPORTTEXTURE_AREA);
	menu->GetSubMenu("File")->GetSubMenu("Export selected textures")->GetSubMenu("Facet by facet")->Add("# of MC Hits", MENU_FILE_EXPORTTEXTURE_MCHITS);
	menu->GetSubMenu("File")->GetSubMenu("Export selected textures")->GetSubMenu("Facet by facet")->Add("Impingement rate (1/s/m\262)", MENU_FILE_EXPORTTEXTURE_IMPINGEMENT);
	menu->GetSubMenu("File")->GetSubMenu("Export selected textures")->GetSubMenu("Facet by facet")->Add("Particle density (1/m\263)", MENU_FILE_EXPORTTEXTURE_PART_DENSITY);
	menu->GetSubMenu("File")->GetSubMenu("Export selected textures")->GetSubMenu("Facet by facet")->Add("Gas density (kg/m\263)", MENU_FILE_EXPORTTEXTURE_GAS_DENSITY);
	menu->GetSubMenu("File")->GetSubMenu("Export selected textures")->GetSubMenu("Facet by facet")->Add("Pressure (mbar)", MENU_FILE_EXPORTTEXTURE_PRESSURE);
	menu->GetSubMenu("File")->GetSubMenu("Export selected textures")->GetSubMenu("Facet by facet")->Add("Avg. Velocity (m/s)", MENU_FILE_EXPORTTEXTURE_AVG_V);
	menu->GetSubMenu("File")->GetSubMenu("Export selected textures")->GetSubMenu("Facet by facet")->Add("Velocity vector (m/s)", MENU_FILE_EXPORTTEXTURE_V_VECTOR);
	menu->GetSubMenu("File")->GetSubMenu("Export selected textures")->GetSubMenu("Facet by facet")->Add("# of velocity vectors", MENU_FILE_EXPORTTEXTURE_N_VECTORS);
	
	menu->GetSubMenu("File")->GetSubMenu("Export selected textures")->Add("By X,Y,Z coordinates");
	menu->GetSubMenu("File")->GetSubMenu("Export selected textures")->GetSubMenu("By X,Y,Z coordinates")->Add("Cell Area (cm\262)", MENU_FILE_EXPORTTEXTURE_AREA_COORD);
	menu->GetSubMenu("File")->GetSubMenu("Export selected textures")->GetSubMenu("By X,Y,Z coordinates")->Add("# of MC Hits", MENU_FILE_EXPORTTEXTURE_MCHITS_COORD);
	menu->GetSubMenu("File")->GetSubMenu("Export selected textures")->GetSubMenu("By X,Y,Z coordinates")->Add("Impingement rate (1/s/m\262)", MENU_FILE_EXPORTTEXTURE_IMPINGEMENT_COORD);
	menu->GetSubMenu("File")->GetSubMenu("Export selected textures")->GetSubMenu("By X,Y,Z coordinates")->Add("Particle density (1/m\263)", MENU_FILE_EXPORTTEXTURE_PART_DENSITY_COORD);
	menu->GetSubMenu("File")->GetSubMenu("Export selected textures")->GetSubMenu("By X,Y,Z coordinates")->Add("Gas density (kg/m\263)", MENU_FILE_EXPORTTEXTURE_GAS_DENSITY_COORD);
	menu->GetSubMenu("File")->GetSubMenu("Export selected textures")->GetSubMenu("By X,Y,Z coordinates")->Add("Pressure (mbar)", MENU_FILE_EXPORTTEXTURE_PRESSURE_COORD);
	menu->GetSubMenu("File")->GetSubMenu("Export selected textures")->GetSubMenu("By X,Y,Z coordinates")->Add("Avg. Velocity (m/s)", MENU_FILE_EXPORTTEXTURE_AVG_V_COORD);
	menu->GetSubMenu("File")->GetSubMenu("Export selected textures")->GetSubMenu("By X,Y,Z coordinates")->Add("Velocity vector (m/s)", MENU_FILE_EXPORTTEXTURE_V_VECTOR_COORD);
	menu->GetSubMenu("File")->GetSubMenu("Export selected textures")->GetSubMenu("By X,Y,Z coordinates")->Add("# of velocity vectors", MENU_FILE_EXPORTTEXTURE_N_VECTORS_COORD);

	menu->GetSubMenu("File")->Add("Import desorption from SYN file",MENU_FILE_IMPORTDES_SYN);
	//menu->GetSubMenu("File")->GetSubMenu("Import desorption file")->Add("SYN file", );
	//menu->GetSubMenu("File")->GetSubMenu("Import desorption file")->Add("DES file (deprecated)", MENU_FILE_IMPORTDES_DES);
	menu->GetSubMenu("File")->Add(NULL); // Separator
	menu->GetSubMenu("File")->Add("Export buffer"); //, MENU_FILE_EXPORTBUFFER);
	menu->GetSubMenu("File")->GetSubMenu("Export buffer")->Add("Export load buffer", MENU_FILE_EXPORTBUFFER_LOAD);
	menu->GetSubMenu("File")->GetSubMenu("Export buffer")->Add("Export hit buffer", MENU_FILE_EXPORTBUFFER_HIT);
	menu->GetSubMenu("File")->Add("&Import buffer");
	menu->GetSubMenu("File")->GetSubMenu("Import buffer")->Add("&Import hit buffer", MENU_FILE_IMPORTBUFFER_HIT);
	menu->GetSubMenu("File")->GetSubMenu("Import buffer")->Add("&Import load buffer", MENU_FILE_IMPORTBUFFER_LOAD);

	menu->GetSubMenu("File")->Add(NULL); // Separator
	menu->GetSubMenu("File")->Add("E&xit", MENU_FILE_EXIT); //Moved here from OnetimeSceneinit_shared to assert it's the last menu item

	menu->GetSubMenu("Selection")->Add(NULL); // Separator
	menu->GetSubMenu("Selection")->Add("Select Desorption", MENU_FACET_SELECTDES);
	menu->GetSubMenu("Selection")->Add("Select Outgassing Map", MENU_SELECT_HASDESFILE);
	menu->GetSubMenu("Selection")->Add("Select Reflective", MENU_FACET_SELECTREFL);
	menu->GetSubMenu("Selection")->Add("Select volatile facets", MENU_FACET_SELECTVOL);

	menu->GetSubMenu("Tools")->Add(NULL);
	menu->GetSubMenu("Tools")->Add("Moving parts...", MENU_TOOLS_MOVINGPARTS);
	menu->GetSubMenu("Facet")->Add("Convert to outgassing map...", MENU_FACET_OUTGASSINGMAP);

	menu->Add("Time");
	menu->GetSubMenu("Time")->Add("Time settings...", MENU_TIME_SETTINGS, SDLK_i, ALT_MODIFIER);
	menu->GetSubMenu("Time")->Add("Edit moments...", MENU_TIME_MOMENTS_EDITOR);
	menu->GetSubMenu("Time")->Add("Edit parameters...", MENU_TIME_PARAMETER_EDITOR);
	menu->GetSubMenu("Time")->Add(NULL);
	menu->GetSubMenu("Time")->Add("Timewise plotter", MENU_TIMEWISE_PLOTTER);
	menu->GetSubMenu("Time")->Add("Pressure evolution", MENU_TIME_PRESSUREEVOLUTION);

	showFilter = new GLToggle(0, "Filtering");
	//togglePanel->Add(showFilter);

	viewerMoreButton = new GLButton(0, "<< View");
	togglePanel->Add(viewerMoreButton);

	shortcutPanel = new GLTitledPanel("Shortcuts");
	shortcutPanel->SetClosable(true);
	shortcutPanel->Close();
	Add(shortcutPanel);

	profilePlotterBtn = new GLButton(0, "Profile pl.");
	shortcutPanel->Add(profilePlotterBtn);

	texturePlotterBtn = new GLButton(0, "Texture pl.");
	shortcutPanel->Add(texturePlotterBtn);

	textureScalingBtn = new GLButton(0, "Tex.scaling");
	shortcutPanel->Add(textureScalingBtn);

	globalSettingsBtn = new GLButton(0, "<< Sim");
	simuPanel->Add(globalSettingsBtn);

	/*
	statusSimu = new GLButton(0,"...");
	simuPanel->Add(statusSimu);
	*/

	modeLabel = new GLLabel("Mode");
	//simuPanel->Add(modeLabel);

	modeCombo = new GLCombo(0);
	modeCombo->SetEditable(true);
	modeCombo->SetSize(2);
	modeCombo->SetValueAt(0, "Monte Carlo");
	modeCombo->SetValueAt(1, "Angular Coef");
	modeCombo->SetSelectedIndex(0);
	//simuPanel->Add(modeCombo);

	compACBtn = new GLButton(0, "Calc AC");
	compACBtn->SetEnabled(false);
	//simuPanel->Add(compACBtn);

	singleACBtn = new GLButton(0, "1");
	singleACBtn->SetEnabled(false);
	//simuPanel->Add(singleACBtn);

	inputPanel = new GLTitledPanel("Particles in");
	facetPanel->Add(inputPanel);

	facetDLabel = new GLLabel("Desorption");
	facetPanel->Add(facetDLabel);
	facetDesType = new GLCombo(0);
	facetDesType->SetSize(5);
	facetDesType->SetValueAt(0, "None");
	facetDesType->SetValueAt(1, "Uniform");
	facetDesType->SetValueAt(2, "Cosine");
	facetDesType->SetValueAt(3, "Cosine^N");
	facetDesType->SetValueAt(4, "Recorded");
	inputPanel->Add(facetDesType);

	facetDesTypeN = new GLTextField(0, NULL);
	facetDesTypeN->SetEditable(false);
	facetPanel->Add(facetDesTypeN);

	facetFILabel = new GLToggle(0, "Outgassing (mbar*l/s):");
	facetFILabel->SetEnabled(false);
	facetFILabel->SetState(true);
	inputPanel->Add(facetFILabel);
	facetFlow = new GLTextField(0, NULL);
	inputPanel->Add(facetFlow);

	facetFIAreaLabel = new GLToggle(1, "Outg/area(mbar*l/s/cm\262):");
	facetFIAreaLabel->SetEnabled(false);
	inputPanel->Add(facetFIAreaLabel);
	facetFlowArea = new GLTextField(0, NULL);
	inputPanel->Add(facetFlowArea);
	facetDesRateLabel = new GLLabel("Desorption rate [Pa m^3/s]:");
	inputPanel->Add(facetDesRateLabel);
	facetDesRate = new GLTextField(0, NULL);
	inputPanel->Add(facetDesRate);

	/*facetUseDesFileLabel = new GLLabel("Desorp. file");
	facetPanel->Add(facetUseDesFileLabel);
	facetUseDesFile = new GLCombo(0);
	facetUseDesFile->SetSize(1);
	facetUseDesFile->SetValueAt(0,"No desorption map");
	inputPanel->Add(facetUseDesFile);*/

	outputPanel = new GLTitledPanel("Contamination Level");
	facetPanel->Add(outputPanel);
	/*
	facetSLabel = new GLLabel("Sticking factor:");
	outputPanel->Add(facetSLabel);
	facetSticking = new GLTextField(0, NULL);
	outputPanel->Add(facetSticking);*/

	/*
	facetPumpingLabel = new GLLabel("Pumping Speed (l/s):");
	outputPanel->Add(facetPumpingLabel);
	facetPumping = new GLTextField(0, NULL);
	outputPanel->Add(facetPumping);*/
	
	facetcoverageLabel = new GLLabel("Coverage [C eq./monolayer]:");
	outputPanel->Add(facetcoverageLabel);
	facetcoverage = new GLTextField(0, NULL); //can probably be optimized in functions below
	outputPanel->Add(facetcoverage);

	facetcoveringLabel = new GLLabel("Covering [C eq. particles]:");
	outputPanel->Add(facetcoveringLabel);
	facetcovering = new GLTextField(0, NULL);
	outputPanel->Add(facetcovering);

	facetTempLabel = new GLLabel("Temperature (\260K):");
	facetPanel->Add(facetTempLabel);
	facetTemperature = new GLTextField(0, NULL);
	facetPanel->Add(facetTemperature);

	facetReLabel = new GLLabel("Profile:");
	facetPanel->Add(facetReLabel);
	facetRecType = new GLCombo(0);
	facetRecType->SetSize(7);
	facetRecType->SetValueAt(0, "None");
	facetRecType->SetValueAt(1, "Pressure, density (\201)");
	facetRecType->SetValueAt(2, "Pressure, density (\202)");
	facetRecType->SetValueAt(3, "Incident angle");
	facetRecType->SetValueAt(4, "Speed distribution");
	facetRecType->SetValueAt(5, "Orthogonal velocity");
	facetRecType->SetValueAt(6, "Tangential velocity");
	facetPanel->Add(facetRecType);

	facetAdvParamsBtn = new GLButton(0, "<< Adv");
	facetPanel->Add(facetAdvParamsBtn);

	facetList = new GLList(0);
	facetList->SetWorker(&worker);
	facetList->SetGrid(true);
	facetList->SetSelectionMode(MULTIPLE_ROW);
	facetList->SetSize(4, 0);
	facetList->SetColumnWidths((int*)cWidth);
	facetList->SetColumnLabels((char **)cName);
	facetList->SetColumnLabelVisible(true);
	facetList->Sortable = true;
	Add(facetList);

	historyBtn = new GLButton(0, "<< history");
	//historyBtn->SetVisible(false);
	simuPanel->Add(historyBtn);
	/*
	historyList = new GLList(0);
	historyList->SetWorker(&worker);
	historyList->SetGrid(true);
	historyList->SetSelectionMode(MULTIPLE_ROW);
	historyList->SetSize(1, 0);
	historyList->AutoSizeColumn();
	historyList->SetColumnLabels();
	historyList->SetColumnLabelVisible(true);
	historyList->Sortable = true;
	Add(historyList);*/



	facetAdvParams = new FacetAdvParams(&worker); //To use its UpdatefacetParams() routines

	LoadParameterCatalog();

	OneTimeSceneInit_shared_post();
	
	return GL_OK;
}

void MolFlow::PlaceComponents() {


	int sx = m_screenWidth - 205-offset_width;
	int sy = 30;

	

	Place3DViewer();

	geomNumber->SetBounds(sx, 3, 202+offset_width, 18);

	// Viewer settings ----------------------------------------
	togglePanel->SetBounds(sx, sy, 202+ offset_width, 112);

	togglePanel->SetCompBounds(showRule, 5, 20, 60, 18);
	togglePanel->SetCompBounds(showNormal, 70+int(offset_width/2), 20, 60, 18);
	togglePanel->SetCompBounds(showUV, 135 + offset_width, 20, 60, 18);

	togglePanel->SetCompBounds(showLine, 5, 42, 60, 18);
	togglePanel->SetCompBounds(showLeak, 70 + int(offset_width / 2), 42, 60, 18);
	togglePanel->SetCompBounds(showHit, 135 + offset_width, 42, 60, 18);

	togglePanel->SetCompBounds(showVolume, 5, 64, 60, 18);
	togglePanel->SetCompBounds(showTexture, 70 + int(offset_width / 2), 64, 60, 18);
	togglePanel->SetCompBounds(showFilter, 135 + offset_width, 64, 60, 18);

	togglePanel->SetCompBounds(viewerMoreButton, 5, 86, 55, 18);
	togglePanel->SetCompBounds(showVertex, 70 + int(offset_width / 2), 86, 60, 18);
	togglePanel->SetCompBounds(showIndex, 135 + offset_width, 86, 60, 18);

	sy += (togglePanel->GetHeight() + 5);

	// Selected facet -----------------------------------------
	int offset_output = 0;
	int offset_input = 25;
	int offset_facet = 150;
	int offset_text = 12;

	facetPanel->SetBounds(sx, sy, 202+ offset_width, 330 + offset_output + offset_input+offset_facet);

	facetPanel->SetCompBounds(inputPanel, 5, 16, 192+ offset_width, 90 + offset_input);

	int cursorY = 15;
	inputPanel->SetCompBounds(facetDLabel, 5, cursorY, 60, 18);
	inputPanel->SetCompBounds(facetDesType, 65 + int(offset_width / 2), cursorY, 80, 18);
	inputPanel->SetCompBounds(facetDesTypeN, 150 + offset_width, cursorY, 30, 18);

	inputPanel->SetCompBounds(facetFILabel, 5, cursorY += 25, 110, 18);
	inputPanel->SetCompBounds(facetFlow, 140 + offset_width, cursorY, 45, 18);

	inputPanel->SetCompBounds(facetFIAreaLabel, 5, cursorY += 25, 110, 18);
	inputPanel->SetCompBounds(facetFlowArea, 140 + offset_width, cursorY, 45, 18);

	inputPanel->SetCompBounds(facetDesRateLabel, 5, cursorY += 25, 110, 18);
	inputPanel->SetCompBounds(facetDesRate, 140 + offset_width, cursorY, 45, 18);

	//inputPanel->SetCompBounds(facetUseDesFileLabel,5,90,60,18);
	//inputPanel->SetCompBounds(facetUseDesFile,65,90,120,18);

	
	//facetPanel->SetCompBounds(outputPanel, 5, cursorY += 45, 192, 65);
	facetPanel->SetCompBounds(outputPanel, 5, cursorY += 45, 192 + offset_width, 65 + offset_output);

	/*outputPanel->SetCompBounds(facetSLabel, 7, cursorY = 15, 100, 18);
	outputPanel->SetCompBounds(facetSticking, 140, cursorY, 45, 18);*/
	cursorY = -10;
	//outputPanel->SetCompBounds(facetPumpingLabel, 7, cursorY += 25, 100, 18);
	//outputPanel->SetCompBounds(facetPumping, 140, cursorY, 45, 18);

	outputPanel->SetCompBounds(facetcoveringLabel, 7, cursorY += 25, 100, 18);
	outputPanel->SetCompBounds(facetcovering, 140 + offset_width, cursorY, 45, 18);

	outputPanel->SetCompBounds(facetcoverageLabel, 7, cursorY += 25, 100, 18);
	outputPanel->SetCompBounds(facetcoverage, 140 + offset_width, cursorY, 45, 18);

	facetPanel->SetCompBounds(facetSideLabel, 7, cursorY = 180 + offset_input+ offset_output, 50, 18);
	facetPanel->SetCompBounds(facetSideType, 65 , cursorY, 130 + offset_width, 18);

	facetPanel->SetCompBounds(facetTLabel, 7, cursorY += 25, 100, 18);
	facetPanel->SetCompBounds(facetOpacity, 110 + offset_text, cursorY, 82 + offset_width - offset_text, 18);

	facetPanel->SetCompBounds(facetTempLabel, 7, cursorY += 25, 100, 18);
	facetPanel->SetCompBounds(facetTemperature, 110 + offset_text, cursorY, 82 + offset_width - offset_text, 18);

	facetPanel->SetCompBounds(facetAreaLabel, 7, cursorY += 25, 100, 18);
	facetPanel->SetCompBounds(facetArea, 110 + offset_text, cursorY, 82 + offset_width - offset_text, 18);

	facetPanel->SetCompBounds(facetSurfaceFactorLabel, 7, cursorY += 25, 100, 18);
	facetPanel->SetCompBounds(facetSurfaceFactor, 110 + offset_text, cursorY, 82 + offset_width - offset_text, 18);

	facetPanel->SetCompBounds(facetDepthLabel, 7, cursorY += 25, 100, 18);
	facetPanel->SetCompBounds(facetDepth, 110 + offset_text, cursorY, 82 + offset_width - offset_text, 18);

	facetPanel->SetCompBounds(facetVolumeLabel, 7, cursorY += 25, 100, 18);
	facetPanel->SetCompBounds(facetVolume, 110 + offset_text, cursorY, 82 + offset_width - offset_text, 18);

	facetPanel->SetCompBounds(facetDiffusionLabel, 7, cursorY += 25, 100, 18);
	facetPanel->SetCompBounds(facetDiffusion, 110 + offset_text, cursorY, 82 + offset_width - offset_text, 18);

	facetPanel->SetCompBounds(facetConcentrationLabel, 7, cursorY += 25, 100, 18);
	facetPanel->SetCompBounds(facetConcentration, 110 + offset_text, cursorY, 82 + offset_width - offset_text, 18);

	facetPanel->SetCompBounds(facetMassLabel, 7, cursorY += 25, 100, 18);
	facetPanel->SetCompBounds(facetMass, 110 + offset_text, cursorY, 82 + offset_width - offset_text, 18);

	facetPanel->SetCompBounds(facetReLabel, 7, cursorY += 25, 60, 18);
	facetPanel->SetCompBounds(facetRecType, 65 , cursorY, 130 + offset_width, 18);

	facetPanel->SetCompBounds(facetAdvParamsBtn, 5, cursorY += 25, 48, 18);
	facetPanel->SetCompBounds(facetDetailsBtn, 56+int(offset_width / 3), cursorY, 45, 18);
	facetPanel->SetCompBounds(facetCoordBtn, 104+int(2*offset_width / 3), cursorY, 45, 18);
	facetPanel->SetCompBounds(facetApplyBtn, 153 + offset_width, cursorY, 44, 18);

	sy += facetPanel->GetHeight() + 5;

	shortcutPanel->SetBounds(sx, sy, 202 + offset_width, 40);
	shortcutPanel->SetCompBounds(profilePlotterBtn, 5 , 15, 60, 18);
	shortcutPanel->SetCompBounds(texturePlotterBtn, 70 + int(offset_width / 2), 15, 60, 18);
	shortcutPanel->SetCompBounds(textureScalingBtn, 135 + offset_width, 15, 60, 18);

	sy += shortcutPanel->GetHeight() + 5;

	// Simulation ---------------------------------------------
	simuPanel->SetBounds(sx, sy, 202 + offset_width, 169+50);

	simuPanel->SetCompBounds(globalSettingsBtn, 5, 20, 48, 19);
	simuPanel->SetCompBounds(startSimu, 58 + int(offset_width / 2), 20, 66, 19);
	simuPanel->SetCompBounds(resetSimu, 128 + offset_width, 20, 66, 19);
	//simuPanel->SetCompBounds(statusSimu,175,20,20,19);
	simuPanel->SetCompBounds(modeLabel, 5, 45, 30, 18);
	simuPanel->SetCompBounds(modeCombo, 40, 45, 85, 18);
	simuPanel->SetCompBounds(compACBtn, 130, 45, 65, 19);
	simuPanel->SetCompBounds(autoFrameMoveToggle, 5, 45, 65, 19);
	simuPanel->SetCompBounds(forceFrameMoveButton, 128 + offset_width, 45, 66, 19);

	//simuPanel->SetCompBounds(compACBtn,123,45,52,19);
	//simuPanel->SetCompBounds(singleACBtn,178,45,20,19);
	simuPanel->SetCompBounds(hitLabel, 5, 70, 30, 18);
	simuPanel->SetCompBounds(hitNumber, 40, 70, 155 + offset_width, 18);
	simuPanel->SetCompBounds(desLabel, 5, 95, 30, 18);
	simuPanel->SetCompBounds(desNumber, 40, 95, 155 + offset_width, 18);
	simuPanel->SetCompBounds(leakLabel, 5, 120, 30, 18);
	simuPanel->SetCompBounds(leakNumber, 40, 120, 155 + offset_width, 18);
	simuPanel->SetCompBounds(sTimeLabel, 5, 145, 30, 18);
	simuPanel->SetCompBounds(sTime, 40, 145, 155 + offset_width, 18);

	simuPanel->SetCompBounds(historyBtn, 5, 145+50, 50, 18);

	
	sy += (simuPanel->GetHeight() + 5);

	
	int lg = m_screenHeight -23 /*- (nbFormula * 25)*/;

	facetList->SetBounds(sx, sy, 202 + offset_width, (lg - sy));
	sy += (facetList->GetHeight() + 5);
	//historyList->SetBounds(sx, sy, 202, (lg - sy));

	/*
	for (int i = 0; i < nbFormula; i++) {
		formulas[i].name->SetBounds(sx, lg + 5, 95, 18);
		formulas[i].value->SetBounds(sx + 90, lg + 5, 87, 18);
		formulas[i].setBtn->SetBounds(sx + 182, lg + 5, 20, 18);
		lg += 25;
	}
	*/
}

// Name: ClearFacetParams()
// Desc: Reset selected facet parameters.

void MolFlow::ClearFacetParams() {
	facetPanel->SetTitle("Selected Facet (none)");
	/*facetSticking->Clear();
	facetSticking->SetEditable(false);*/
	facetFILabel->SetEnabled(false);
	facetFIAreaLabel->SetEnabled(false);
	facetFlow->Clear();
	facetFlow->SetEditable(false);
	facetFlowArea->Clear();
	facetFlowArea->SetEditable(false);
	facetDesRate->Clear();
	facetDesRate->SetEditable(false);
	facetArea->SetEditable(false);
	facetArea->Clear();
	facetSurfaceFactor->SetEditable(false);
	facetSurfaceFactor->Clear();
	facetDepth->SetEditable(false);
	facetDepth->Clear();
	facetVolume->SetEditable(false);
	facetVolume->Clear();
	facetDiffusion->SetEditable(false);
	facetDiffusion->Clear();
	facetConcentration->SetEditable(false);
	facetConcentration->Clear();
	facetMass->SetEditable(false);
	facetMass->Clear();
	/*facetPumping->SetEditable(false);
	facetPumping->Clear();*/
	facetcoverage->SetEditable(false);
	facetcoverage->Clear();
	facetcovering->SetEditable(false);
	facetcovering->Clear();
	facetOpacity->Clear();
	facetOpacity->SetEditable(false);
	facetTemperature->Clear();
	facetTemperature->SetEditable(false);
	facetSideType->SetSelectedValue("");
	facetSideType->SetEditable(false);
	facetDesType->SetSelectedValue("");
	facetDesType->SetEditable(false);
	facetDesTypeN->SetText("");
	facetDesTypeN->SetEditable(false);
	facetRecType->SetSelectedValue("");
	facetRecType->SetEditable(false);
}

// Name: ApplyFacetParams()
// Desc: Apply facet parameters.

void MolFlow::ApplyFacetParams() {

	Geometry *geom = worker.GetGeometry();
	size_t nbFacet = geom->GetNbFacet();
	/*
	// Sticking
	double sticking;
	bool stickingNotNumber;
	bool doSticking = false;
	if (facetSticking->GetNumber(&sticking)) {
		if (sticking<0.0 || sticking>1.0) {
			GLMessageBox::Display("Sticking must be in the range [0,1]", "Error", GLDLG_OK, GLDLG_ICONERROR);
			return;
		}
		doSticking = true;
		stickingNotNumber = false;
	}
	else {
		if (facetSticking->GetText() == "...") doSticking = false;
		//else {
		//	GLMessageBox::Display("Invalid sticking number","Error",GLDLG_OK,GLDLG_ICONERROR);
		//	UpdateFacetParams();
		//	return;
			doSticking = true;
			stickingNotNumber = true;
		}
	}*/

	// DesRate
	double DesRate;
	bool DesRateNotNumber;
	bool doDesRate = false;
	if (facetDesRate->GetNumber(&DesRate)) {
		if (DesRate <0.0 ) {
			GLMessageBox::Display("Desorption Rate must be positive", "Error", GLDLG_OK, GLDLG_ICONERROR);
			return;
		}
		doDesRate = true;
		DesRateNotNumber = false;
	}
	else {
		if (facetDesRate->GetText() == "...") doDesRate = false;
		else {/*
			GLMessageBox::Display("Invalid sticking number","Error",GLDLG_OK,GLDLG_ICONERROR);
			UpdateFacetParams();
			return;*/
			doDesRate = true;
			DesRateNotNumber = true;
		}
	}

	//facetArea
	double surfacefactor;
	double surfacefactorNotNumber;
	double dosurfacefactor = false;

	if (facetSurfaceFactor->GetNumber(&surfacefactor)) {
		if (surfacefactor < 1.0) {
			GLMessageBox::Display("Surfacefactor must be >= 1", "Error", GLDLG_OK, GLDLG_ICONERROR);
			return;
		}
		dosurfacefactor = true;
		surfacefactorNotNumber = false;
	}
	else {
		if (facetSurfaceFactor->GetText() == "...") dosurfacefactor = false;
		else {/*
			GLMessageBox::Display("Invalid sticking number","Error",GLDLG_OK,GLDLG_ICONERROR);
			UpdateFacetParams();
			return;*/
			dosurfacefactor = true;
			surfacefactorNotNumber = true;
		}
	}

	double depth;
	double depthNotNumber;
	double doDepth = false;

	if (facetDepth->GetNumber(&depth)) {
		if (depth < 0.0) {
			GLMessageBox::Display("Facet Depth must be positve", "Error", GLDLG_OK, GLDLG_ICONERROR);
			return;
		}
		doDepth = true;
		depthNotNumber = false;
	}
	else {
		if (facetDepth->GetText() == "...") doDepth = false;
		else {/*
			GLMessageBox::Display("Invalid sticking number","Error",GLDLG_OK,GLDLG_ICONERROR);
			UpdateFacetParams();
			return;*/
			doDepth = true;
			depthNotNumber = true;
		}
	}

	double diffusion;
	double diffusionNotNumber;
	double doDiffusion = false;

	if (facetDiffusion->GetNumber(&diffusion)) {
		if (diffusion < 0.0) {
			GLMessageBox::Display("Facet Diffusion must be positve", "Error", GLDLG_OK, GLDLG_ICONERROR);
			return;
		}
		doDiffusion = true;
		diffusionNotNumber = false;
	}
	else {
		if (facetDiffusion->GetText() == "...") doDiffusion = false;
		else {/*
			GLMessageBox::Display("Invalid sticking number","Error",GLDLG_OK,GLDLG_ICONERROR);
			UpdateFacetParams();
			return;*/
			doDiffusion = true;
			diffusionNotNumber = true;
		}
	}

	double concentration;
	double concentrationNotNumber;
	double doConcentration = false;

	if (facetConcentration->GetNumber(&concentration)) {
		if (concentration < 0.0) {
			GLMessageBox::Display("Facet Concentration must be positve", "Error", GLDLG_OK, GLDLG_ICONERROR);
			return;
		}
		doConcentration = true;
		concentrationNotNumber = false;
	}
	else {
		if (facetDiffusion->GetText() == "...") doConcentration = false;
		else {/*
			GLMessageBox::Display("Invalid sticking number","Error",GLDLG_OK,GLDLG_ICONERROR);
			UpdateFacetParams();
			return;*/
			doConcentration = true;
			concentrationNotNumber = true;
		}
	}

	// Coverage and Covering
	// coverage: one monolayer of contaminating particles means a coverage of 1,0
	// covering will cout each contaminating parcticle

	double coverage;
	llong covering;
	bool coverageNotNumber;
	bool docoverage = false;
	bool coveringNotNumber;
	bool docovering = false;

	if (facetcoverage->GetNumber(&coverage)) {
		if (coverage<0.0) {
			GLMessageBox::Display("Coverage must be positive", "Error", GLDLG_OK, GLDLG_ICONERROR);
			return;
		}
		//covering = llong(abs(coverage * calcNmono()/calcdNsurf()));
		docoverage = true;
		coverageNotNumber = false;
	}
	else {
		if (facetcoverage->GetText() == "...") docoverage = false;
		else {/*
			GLMessageBox::Display("Invalid sticking number","Error",GLDLG_OK,GLDLG_ICONERROR);
			UpdateFacetParams();
			return;*/
			docoverage = true;
			coverageNotNumber = true;
		}
	}

	if (facetcovering->GetNumberSizeT(&covering)) {
		if (covering < 0) {
			GLMessageBox::Display("Covering must be positive", "Error", GLDLG_OK, GLDLG_ICONERROR);
			return;
		}
		//covering = llong(abs(coverage * calcNmono() / calcdNsurf()));
		docovering = true;
		coveringNotNumber = false;
	}
	else {
		if (facetcovering->GetText() == "...") docovering = false;
		else {/*
			GLMessageBox::Display("Invalid sticking number","Error",GLDLG_OK,GLDLG_ICONERROR);
			UpdateFacetParams();
			return;*/
			docovering = true;
			coveringNotNumber = true;
		}
	}

	// opacity
	double opacity;
	bool doOpacity = false;
	bool opacityNotNumber;
	if (facetOpacity->GetNumber(&opacity)) {
		if (opacity<0.0 || opacity>1.0) {
			GLMessageBox::Display("Opacity must be in the range [0,1]", "Error", GLDLG_OK, GLDLG_ICONERROR);
			return;
		}
		doOpacity = true;
		opacityNotNumber = false;
	}
	else {
		if (facetOpacity->GetText() == "...") doOpacity = false;
		else {/*
			GLMessageBox::Display("Invalid opacity number","Error",GLDLG_OK,GLDLG_ICONERROR);
			UpdateFacetParams();
			return;*/
			doOpacity = true;
			opacityNotNumber = true;
		}
	}

	// temperature
	double temperature;
	bool doTemperature = false;
	if (facetTemperature->GetNumber(&temperature)) {
		if (temperature < 0.0) {
			GLMessageBox::Display("Temperature must be positive", "Error", GLDLG_OK, GLDLG_ICONERROR);
			return;
		}
		doTemperature = true;
	}
	else {
		if (facetTemperature->GetText() == "...") doTemperature = false;
		else {
			GLMessageBox::Display("Invalid temperature number", "Error", GLDLG_OK, GLDLG_ICONERROR);
			return;
		}
	}

	// Outgassing
	double outgassing = 0.0;
	bool doFlow = false;
	bool outgassingNotNumber;
	//Calculate outgassing
	if (facetFILabel->GetState() && facetFlow->GetText() != "..." && facetDesType->GetSelectedIndex() != 0
		&& facetDesType->GetSelectedValue() != "..." && facetFlow->IsEditable()) {  //We want outgassing
		if (facetFlow->GetNumber(&outgassing)) { //If we can parse the number
			if (!(outgassing > 0.0)) {
				GLMessageBox::Display("Outgassing must be positive", "Error", GLDLG_OK, GLDLG_ICONERROR);
				return;
			}
			doFlow = true;
			outgassingNotNumber = false;
		}
		else { //could not parse as number
			doFlow = true;
			outgassingNotNumber = true;
		}
	}

	// Outgassing per area
	double flowA = 0;
	bool doFlowA = false;
	//Calculate outgassing

	if (facetFIAreaLabel->GetState() && facetFlowArea->GetText()!= "..."
		&& facetDesType->GetSelectedIndex() != 0 && facetDesType->GetSelectedValue() != "..." && facetFlowArea->IsEditable()) { //We want outgassing per area
		if (facetFlowArea->GetNumber(&flowA)) { //Can be parsed as number
			if (!(flowA > 0.0)) {
				GLMessageBox::Display("Outgassing per area must be positive", "Error", GLDLG_OK, GLDLG_ICONERROR);
				return;
			}
			doFlowA = true;
		}
		else {
			GLMessageBox::Display("Invalid outgassing per area number", "Error", GLDLG_OK, GLDLG_ICONERROR);
			return;
		}
	}

	// Desorption type
	int desorbType = facetDesType->GetSelectedIndex();

	double desorbTypeN;
	bool doDesorbTypeN = false;
	if (desorbType == 3) {
		if (facetDesTypeN->GetNumber(&desorbTypeN)) {
			if (!(desorbTypeN > 0.0)) {
				GLMessageBox::Display("Desorption type exponent must be greater than 0.0", "Error", GLDLG_OK, GLDLG_ICONERROR);
				UpdateFacetParams();
				return;
			}
			doDesorbTypeN = true;
		}
		else {
			if (facetDesTypeN->GetText() == "...") doDesorbTypeN = false;
			else {
				GLMessageBox::Display("Invalid desorption type exponent", "Error", GLDLG_OK, GLDLG_ICONERROR);
				return;
			}
		}
	}

	// Record (profile) type
	int rType = facetRecType->GetSelectedIndex(); // -1 if "..."

	// 2sided
	int is2Sided = facetSideType->GetSelectedIndex();

	//Check complete, let's apply
	if (facetAdvParams && facetAdvParams->IsVisible()) {
		if (!facetAdvParams->Apply()) {
			return;
		}
	}
	if (!AskToReset()) return;
	
	changedSinceSave = true;

	// Update facets (local)
	for (int i = 0; i < nbFacet; i++) {
		Facet *f = geom->GetFacet(i);
		if (f->selected) {
			/*
			if (doSticking) {
				if (!stickingNotNumber) {
					f->sh.sticking = sticking;
					f->userSticking = "";
				}
				else {
					f->userSticking = facetSticking->GetText();
				}
			}*/
			if (dosurfacefactor) {
				if (!surfacefactorNotNumber) {
					f->sh.effectiveSurfaceFactor = surfacefactor;
					f->userSurfacefactor = "";
				}
				else {
					f->userSurfacefactor = facetSurfaceFactor->GetText();
				}
			}

			if (doDepth) {
				if (!depthNotNumber) {
					f->sh.facetDepth = depth;
					f->userSurfacefactor = "";
				}
				else {
					f->userDepth = facetDepth->GetText();
				}
			}

			if (doDiffusion) {
				if (!diffusionNotNumber) {
					f->sh.diffusionCoefficient = diffusion;
					f->userDiffusion = "";
				}
				else {
					f->userDiffusion = facetDiffusion->GetText();
				}
			}

			if (doConcentration) {
				if (!concentrationNotNumber) {
					f->sh.concentration = concentration;
					f->userConcentration = "";
				}
				else {
					f->userConcentration = facetConcentration->GetText();
				}
			}

			if (doDesRate) {
				if (!DesRateNotNumber) {
					f->sh.desorption = DesRate;
					f->userdesorption = "";
				}
				else {
					f->userdesorption = facetDesRate->GetText();
				}
			}

			if (docoverage) {
				if (!coverageNotNumber) {
					f->facetHitCache.hit.covering = covering;
					f->usercoverage = "";
				}
				else {
					f->usercoverage = facetcoverage->GetText();
				}
			}

			if (docovering) {
				if (!covering) {
					f->facetHitCache.hit.covering = covering;
					f->usercovering = "";
				}
				else {
					f->usercovering = facetcovering->GetText();
				}
			}

			if (doOpacity) {
				if (!opacityNotNumber) {
					f->sh.opacity = opacity;
					f->userOpacity = "";
				}
				else {
					f->userOpacity = facetOpacity->GetText();
				}
			}
			if (doTemperature) f->sh.temperature = temperature;
			if (doFlow) {
				if (!outgassingNotNumber) {
					f->sh.outgassing = outgassing*0.100; //0.1: mbar*l/s -> Pa*m3/s
					f->userOutgassing = "";
				}
				else {
					f->userOutgassing = facetFlow->GetText();
				}
			}
			if (doFlowA/* && !useMapA*/) f->sh.outgassing = flowA*f->GetArea()*0.100;
			if (desorbType >= 0) {
				if (desorbType == 0) f->sh.outgassing = 0.0;
				if (desorbType != 3) f->sh.desorbTypeN = 0.0;
				f->sh.desorbType = desorbType;
				if (doDesorbTypeN) f->sh.desorbTypeN = desorbTypeN;
			}

			if (rType >= 0) {
				f->sh.profileType = rType;
				//f->wp.isProfile = (rType!=PROFILE_NONE); //included below by f->UpdateFlags();
			}
			if (is2Sided >= 0) f->sh.is2sided = is2Sided;

			f->sh.maxSpeed = 4.0 * sqrt(2.0*8.31*f->sh.temperature / 0.001 / worker.wp.gasMass);
			f->UpdateFlags();
		}
	}

	// Mark "needsReload" to sync changes with workers on next simulation start
	try { worker.Reload(); }
	catch (Error &e) {
		GLMessageBox::Display((char *)e.GetMsg(), "Error", GLDLG_OK, GLDLG_ICONERROR);
		return;
	}
	worker.CalcTotalOutgassing();
	UpdateFacetParams();
	if (profilePlotter) profilePlotter->Refresh();
	if (pressureEvolution) pressureEvolution->Refresh();
	if (timewisePlotter) timewisePlotter->Refresh();
	//if (facetAdvParams) facetAdvParams->Refresh();
}

// Name: UpdateFacetParams()
// Desc: Update selected facet parameters.

void MolFlow::UpdateFacetParams(bool updateSelection) { //Calls facetAdvParams->Refresh()

	char tmp[256];

	// Update params
	Geometry *geom = worker.GetGeometry();
	// Get list of selected facet
	std::vector<size_t> selectedFacets = geom->GetSelectedFacets();
	size_t nbSel = selectedFacets.size();

	if (nbSel > 0) {

		Facet *f0;
		Facet *f;

		

		f0 = geom->GetFacet(selectedFacets[0]);

		double f0Area = f0->GetArea();
		double sumArea = f0Area; //sum facet area
		double sumVolume = f0Area * f0->sh.facetDepth * (f0->sh.is2sided ? 0.5 : 1); //sum facet area
		double sumMass = sumVolume * f0->sh.concentration;

		bool stickingE = true;
		bool opacityE = true;
		bool temperatureE = true;
		bool flowE = true;
		bool flowAreaE = true;
		bool desorbTypeE = true;
		bool desorbTypeNE = true;
		bool recordE = true;
		bool is2sidedE = true;
		bool coverageE = true;
		bool coveringE = true;
		bool desrateE = true;
		bool surfaceFactorE = true;
		bool facetDepthE = true;
		bool facetDiffusionE = true;
		bool facetConcentrationE = true;

 		for (size_t sel = 1; sel < selectedFacets.size();sel++) {
			f = geom->GetFacet(selectedFacets[sel]);
			double fArea = f->GetArea();
			double fVolume= fArea * f->sh.facetDepth * (f->sh.is2sided ? 0.5 : 1);
			//stickingE = stickingE && (f0->userSticking.compare(f->userSticking) == 0) && IsEqual(f0->sh.sticking, f->sh.sticking);
			opacityE = opacityE && (f0->userOpacity.compare(f->userOpacity) == 0) && IsEqual(f0->sh.opacity, f->sh.opacity);
			temperatureE = temperatureE && IsEqual(f0->sh.temperature, f->sh.temperature);
			flowE = flowE && f0->userOutgassing.compare(f->userOutgassing) == 0 && IsEqual(f0->sh.outgassing, f->sh.outgassing);
			flowAreaE = flowAreaE && IsEqual(f0->sh.outgassing / f0Area, f->sh.outgassing / fArea);
			is2sidedE = is2sidedE && (f0->sh.is2sided == f->sh.is2sided);
			desorbTypeE = desorbTypeE && (f0->sh.desorbType == f->sh.desorbType);
			desorbTypeNE = desorbTypeNE && IsEqual(f0->sh.desorbTypeN, f->sh.desorbTypeN);
			recordE = recordE && (f0->sh.profileType == f->sh.profileType);  //profiles
			sumArea += fArea;
			sumVolume += fVolume;
			sumMass += fVolume * f->sh.concentration;
			coverageE= coverageE && (f0->usercoverage.compare(f->usercoverage) == 0) && IsEqual(f0->facetHitCache.hit.covering, f->facetHitCache.hit.covering);
			coveringE = coveringE && (f0->usercovering.compare(f->usercovering) == 0) && IsEqual(f0->facetHitCache.hit.covering, f->facetHitCache.hit.covering);
			desrateE= desrateE && (f0->userdesorption.compare(f->userdesorption) == 0) && IsEqual(f0->sh.desorption, f->sh.desorption);
			surfaceFactorE = surfaceFactorE &&(f0->userSurfacefactor.compare(f->userSurfacefactor)==0)&& IsEqual(f0->sh.effectiveSurfaceFactor, f->sh.effectiveSurfaceFactor);
			facetDepthE = facetDepthE && (f0->userDepth.compare(f->userDepth) == 0) && IsEqual(f0->sh.facetDepth, f->sh.facetDepth);
			facetDiffusionE = facetDiffusionE && (f0->userDiffusion.compare(f->userDiffusion) == 0) && IsEqual(f0->sh.diffusionCoefficient, f->sh.diffusionCoefficient);
			facetConcentrationE = facetConcentrationE && (f0->userConcentration.compare(f->userConcentration) == 0) && IsEqual(f0->sh.concentration, f->sh.concentration);
		}

		if (nbSel == 1)
			sprintf(tmp, "Selected Facet (#%zd)", selectedFacets[0] + 1);
		else
			sprintf(tmp, "Selected Facet (%zd selected)", selectedFacets.size());

		// Old STR compatibility
		//if (stickingE && f0->wp.superDest) stickingE = false;

		facetPanel->SetTitle(tmp);
		if (selectedFacets.size() > 1) { facetAreaLabel->SetText("Sum Area (cm\262):"); facetVolumeLabel->SetText("Sum Volume (cm\263):"); facetMass->SetText("Mass (kg):"); }
		else { facetAreaLabel->SetText("Area (cm\262):"); facetVolumeLabel->SetText("Facet Volume (cm\263):"); facetMass->SetText("Sum Mass (kg):");}
		facetArea->SetText(sumArea);
		facetVolume->SetText(sumVolume);
		facetMass->SetText(sumMass);

		/*if (stickingE) {
			if (f0->userSticking.length() == 0)
				facetSticking->SetText(f0->sh.sticking);
			else facetSticking->SetText(f0->userSticking.c_str());
		}
		else facetSticking->SetText("...");*/

		if (facetDiffusionE) {
			if (f0->userDiffusion.length() == 0) {
				f0->sh.diffusionCoefficient = f0->sh.diffusionCoefficient < 0.0 ? 0.0 : f0->sh.diffusionCoefficient;
				facetDiffusion->SetText(f0->sh.diffusionCoefficient);
			}
			else facetDiffusion->SetText(f0->userDiffusion.c_str());
		}
		else { facetDiffusion->SetText("..."); }

		if (facetConcentrationE) {
			facetMass->SetEditable(true);
			if (f0->userConcentration.length() == 0) {
				f0->sh.concentration = f0->sh.concentration < 0.0 ? 0.0 : f0->sh.concentration;
				facetConcentration->SetText(f0->sh.concentration);
			}
			else facetConcentration->SetText(f0->userConcentration.c_str());
		}
		else { facetConcentration->SetText("..."); facetMass->SetEditable(false);}

		if (surfaceFactorE) {
			if (f0->userSurfacefactor.length() == 0) {
				f0->sh.effectiveSurfaceFactor= f0->sh.effectiveSurfaceFactor < 1.0 ? 1.0 : f0->sh.effectiveSurfaceFactor;
				facetSurfaceFactor->SetText(f0->sh.effectiveSurfaceFactor);
			}
			else facetSurfaceFactor->SetText(f0->userSurfacefactor.c_str());
		}
		else { facetSurfaceFactor->SetText("...");}

		if (facetDepthE) {
			facetVolume->SetEditable(true);
			if (f0->userDepth.length() == 0) {
				f0->sh.facetDepth= f0->sh.facetDepth < 0.0 ? 0.0 : f0->sh.facetDepth;
				facetDepth->SetText(f0->sh.facetDepth);
			}
			else facetDepth->SetText(f0->userDepth.c_str());
		}
		else { facetDepth->SetText("..."); facetVolume->SetEditable(false); }

		if (desrateE) {
			if (f0->userdesorption.length() == 0) {
				f0->sh.desorption = f0->sh.desorption < 0.0 ? 0.0 : f0->sh.desorption;
				facetDesRate->SetText(f0->sh.desorption);
			}
			else facetDesRate->SetText(f0->userdesorption.c_str());
		}
		else facetDesRate->SetText("...");

		if (opacityE) {
			if (f0->userOpacity.length() == 0)
				facetOpacity->SetText(f0->sh.opacity);
			else facetOpacity->SetText(f0->userOpacity.c_str());
		}
		else facetOpacity->SetText("...");

		if (coverageE) {
			if (f0->usercoverage.length() == 0)
				facetcoverage->SetText(double((f0->facetHitCache.hit.covering)/(calcNmono()/calcdNsurf())));
			else facetcoverage->SetText(f0->usercoverage.c_str());
		}
		else facetcoverage->SetText("...");

		if (coveringE) {
			if (f0->usercovering.length() == 0)
				facetcovering->SetText(llong((f0->facetHitCache.hit.covering)));
			else facetcovering->SetText(f0->usercovering.c_str());
		}
		else facetcovering->SetText("...");

		if (temperatureE) facetTemperature->SetText(f0->sh.temperature); else facetTemperature->SetText("...");
		if (is2sidedE) facetSideType->SetSelectedIndex(f0->sh.is2sided); else facetSideType->SetSelectedValue("...");
		if (desorbTypeNE) facetDesTypeN->SetText(f0->sh.desorbTypeN); else facetDesTypeN->SetText("...");
		if (recordE) facetRecType->SetSelectedIndex(f0->sh.profileType); else facetRecType->SetSelectedValue("...");

		/*
		if (selectedFacets.size() == 1) {
			facetPumping->SetFocusable(true);
			//calcFlow();
			facetPumping->SetText(0);


		}
		else {
			facetPumping->SetFocusable(false);
			facetPumping->SetText("...");
		}*/

		if (desorbTypeE) {
			facetDesType->SetSelectedIndex(f0->sh.desorbType);
			if (f0->sh.desorbType > DES_NONE) { //There is some desorption

				//facetFlow->SetEnabled(true);
				facetFIAreaLabel->SetEnabled(true);
				facetFlowArea->SetEditable(true);
				facetFILabel->SetEnabled(true);
				facetFlow->SetEditable(true);
				if (flowE) {
					if (f0->userOutgassing.length() == 0)
						facetFlow->SetText(f0->sh.outgassing*10.00); //10: Pa*m3/sec -> mbar*l/s
					else facetFlow->SetText(f0->userOutgassing);
				}
				else facetFlow->SetText("...");
				if (flowAreaE) facetFlowArea->SetText(f0->sh.outgassing / f0Area*10.00); else facetFlowArea->SetText("...");
				if (f0->sh.desorbType == 3) {
					facetDesTypeN->SetEditable(true);
					if (desorbTypeNE) facetDesTypeN->SetText(f0->sh.desorbTypeN); else facetDesTypeN->SetText("...");
				}
				else {
					facetDesTypeN->SetText("");
					facetDesTypeN->SetEditable(false);
				};

			}
			else { //No desorption
				facetFILabel->SetEnabled(false);
				facetFlow->SetEditable(false);
				facetFIAreaLabel->SetEnabled(false);
				facetFlowArea->SetEditable(false);
				facetDesTypeN->SetText("");
				facetDesTypeN->SetEditable(false);
				facetFlow->SetText("");
				facetFlowArea->SetText("");
			}
		}
		else { //Mixed state
			facetDesType->SetSelectedValue("...");
			facetFILabel->SetEnabled(false);
			facetFlow->SetEditable(false);
			facetFIAreaLabel->SetEnabled(false);
			facetFlowArea->SetEditable(false);
			facetDesTypeN->SetText("");
			facetDesTypeN->SetEditable(false);
			facetFlow->SetText("");
			facetFlowArea->SetText("");
		}

		if (facetAdvParams) facetAdvParams->Refresh(selectedFacets); //Refresh advanced facet parameters panel
		if (updateSelection) {

			if (nbSel > 1000 || geom->GetNbFacet() > 50000) { //If it would take too much time to look up every selected facet in the list
				facetList->ReOrder();
				facetList->SetSelectedRows(selectedFacets, false);
			}
			else {
				facetList->SetSelectedRows(selectedFacets, true);
			}
			facetList->lastRowSel = -1;
		}

		//Enabled->Editable
		facetDesRate->SetFocusable(true);
		facetcoverage->SetEditable(true);
		facetcovering->SetEditable(true);
		//facetSticking->SetFocusable(true);
		facetOpacity->SetEditable(true);
		facetTemperature->SetEditable(true);
		facetSideType->SetEditable(true);
		facetDesType->SetEditable(true);
		facetRecType->SetEditable(true);
		facetApplyBtn->SetEnabled(false);
		facetSurfaceFactor->SetEditable(true);
		facetDepth->SetEditable(true);
		facetDiffusion->SetEditable(true);
		facetConcentration->SetEditable(true);
	}
	else {
		ClearFacetParams();
		if (facetAdvParams) facetAdvParams->Refresh(std::vector<size_t>()); //Clear
		if (updateSelection) facetList->ClearSelection();
	}

	if (facetDetails) facetDetails->Update();
	if (facetCoordinates) facetCoordinates->UpdateFromSelection();
	if (texturePlotter) texturePlotter->Update(m_fTime, true); //Facet change
	if (outgassingMap) outgassingMap->Update(m_fTime, true);
	if (histogramSettings) histogramSettings->Refresh(selectedFacets);
}

/*
void MolFlow::LogProfile() {

Geometry *geom = worker.GetGeometry();
BYTE *buffer = worker.GetHits();
if(!buffer) return;

char filename[256];
sprintf(filename,"C:\\Temp\\dataS%d.txt",nbSt);
FILE *file = fopen(filename,"a");

SHGHITS *gHits = (SHGHITS *)buffer;
double nbAbs = (double)gHits->globalHits.hit.nbAbsEquiv;
double nbDes = (double)gHits->globalHits.hit.nbDesorbed;
double nbMCHit = (double)gHits->globalHits.hit.nbMCHit;

fprintf(file,"Time:%s Sticking=%g Des=%g\n",FormatTime(worker.simuTime),(double)nbSt/10.0,nbDes);

// Volatile profile
int nb = geom->GetNbFacet();
for(int j=0;j<nb;j++) {
Facet *f = geom->GetFacet(j);
if( f->wp.isVolatile ) {
FacetHitBuffer *fCount = (FacetHitBuffer *)(buffer + f->wp.hitOffset);
double z = geom->GetVertex(f->indices[0])->z;
fprintf(file,"%g %.10g\n",z,(double)(fCount->hit.nbAbsEquiv)/nbDes);
}
}

// Last
Facet *f = geom->GetFacet(28);
FacetHitBuffer *fCount = (FacetHitBuffer *)(buffer + f->wp.hitOffset);
double fnbAbs = (double)fCount->hit.nbAbsEquiv;
fprintf(file,"1000 %.10g\n",fnbAbs/nbDes);

fclose(file);

worker.ReleaseHits();

}
*/

// Name: FrameMove()
// Desc: Called once per frame, the call is the entry point for animating
//       the scene.

int MolFlow::FrameMove()
{	
	if (worker.isRunning && ((m_fTime - lastUpdate) >= 1.0f)) {
		if (textureScaling) textureScaling->Update();
		//if (formulaEditor && formulaEditor->IsVisible()) formulaEditor->Refresh(); //Interface::Framemove does it already
	}
	Interface::FrameMove(); //might reset lastupdate
	char tmp[256];
	if (globalSettings) globalSettings->SMPUpdate();

	if ((m_fTime - worker.startTime <= 2.0f) && worker.isRunning) {
		hitNumber->SetText("Starting...");
		desNumber->SetText("Starting...");
	}
	else {
		if (worker.wp.sMode == AC_MODE) {
			hitNumber->SetText("");
		}
		else {
			sprintf(tmp, "%s (%s)", FormatInt(worker.globalHitCache.globalHits.hit.nbMCHit, "hit"), FormatPS(hps, "hit"));
			hitNumber->SetText(tmp);
		}
		sprintf(tmp, "%s (%s)", FormatInt(worker.globalHitCache.globalHits.hit.nbDesorbed, "des"), FormatPS(dps, "des"));
		desNumber->SetText(tmp);
	}

	if (worker.calcAC) {
		sprintf(tmp, "Calc AC: %s (%zd %%)", FormatTime(worker.simuTime + (m_fTime - worker.startTime)),
			worker.calcACprg);
	}
	else {
		sprintf(tmp, "Running: %s", FormatTime(worker.simuTime + (m_fTime - worker.startTime)));
	}
	return GL_OK;
}

// Name: RestoreDeviceObjects()
// Desc: Initialize scene objects.

int MolFlow::RestoreDeviceObjects()
{
	RestoreDeviceObjects_shared();

	//Different Molflow implementations:
	RVALIDATE_DLG(facetAdvParams);
	RVALIDATE_DLG(facetDetails);
	RVALIDATE_DLG(smartSelection);
	RVALIDATE_DLG(viewer3DSettings);
	RVALIDATE_DLG(textureScaling);
	RVALIDATE_DLG(globalSettings);
	RVALIDATE_DLG(profilePlotter);
	RVALIDATE_DLG(texturePlotter);

	//Molflow only:
	RVALIDATE_DLG(importDesorption);
	RVALIDATE_DLG(timeSettings);
	RVALIDATE_DLG(movement);
	RVALIDATE_DLG(outgassingMap);
	RVALIDATE_DLG(parameterEditor);
	RVALIDATE_DLG(pressureEvolution);
	RVALIDATE_DLG(timewisePlotter);

	return GL_OK;
}

// Name: InvalidateDeviceObjects()
// Desc: Free all alocated resource

int MolFlow::InvalidateDeviceObjects()
{
	InvalidateDeviceObjects_shared();

	//Different Molflow implementations:
	IVALIDATE_DLG(facetAdvParams);
	IVALIDATE_DLG(facetDetails);
	IVALIDATE_DLG(smartSelection);
	IVALIDATE_DLG(viewer3DSettings);
	IVALIDATE_DLG(textureScaling);
	IVALIDATE_DLG(globalSettings);
	IVALIDATE_DLG(profilePlotter);
	IVALIDATE_DLG(texturePlotter);

	//Molflow only:
	IVALIDATE_DLG(importDesorption);
	IVALIDATE_DLG(timeSettings);
	IVALIDATE_DLG(movement);
	IVALIDATE_DLG(outgassingMap);
	IVALIDATE_DLG(parameterEditor);
	IVALIDATE_DLG(pressureEvolution);
	IVALIDATE_DLG(timewisePlotter);

	return GL_OK;
}

void MolFlow::ExportProfiles() {

	Geometry *geom = worker.GetGeometry();
	if (geom->GetNbSelectedFacets() == 0) {
		GLMessageBox::Display("Empty selection", "Error", GLDLG_OK, GLDLG_ICONERROR);
		return;
	}
	if (!worker.IsDpInitialized()) {
		GLMessageBox::Display("Worker Dataport not initialized yet", "Error", GLDLG_OK, GLDLG_ICONERROR);
		return;
	}

	FILENAME *fn = GLFileBox::SaveFile(currentDir, NULL, "Save File", fileProfFilters, 0);

	if (fn) {
		try {
			worker.ExportProfiles(fn->fullName);
		}
		catch (Error &e) {
			char errMsg[512];
			sprintf(errMsg, "%s\nFile:%s", e.GetMsg(), fn->fullName);
			GLMessageBox::Display(errMsg, "Error", GLDLG_OK, GLDLG_ICONERROR);
		}
	}
}

void MolFlow::ExportAngleMaps() {

	Geometry *geom = worker.GetGeometry();
	std::vector<size_t> angleMapFacetIndices;
	for (size_t i = 0; i < geom->GetNbFacet(); i++) {
		Facet* f = geom->GetFacet(i);
		if (f->selected && f->sh.anglemapParams.hasRecorded)
			angleMapFacetIndices.push_back(i);
	}
	if (angleMapFacetIndices.size() == 0) {
		GLMessageBox::Display("Select at least one facet with recorded angle map", "Error", GLDLG_OK, GLDLG_ICONERROR);
		return;
	}

	if (!worker.IsDpInitialized()) {
		GLMessageBox::Display("Worker Dataport not initialized yet", "Error", GLDLG_OK, GLDLG_ICONERROR);
		return;
	}

	FILENAME *fn = GLFileBox::SaveFile(currentDir, NULL, "Save File", fileProfFilters, 0);

	if (fn) {
		try {
			std::string fileName = fn->fullName;
			worker.ExportAngleMaps(angleMapFacetIndices, fileName);
		}
		catch (Error &e) {
			char errMsg[512];
			sprintf(errMsg, "%s\nFile:%s", e.GetMsg(), fn->fullName);
			GLMessageBox::Display(errMsg, "Error", GLDLG_OK, GLDLG_ICONERROR);
			return;
		}
	}
}

void MolFlow::ImportAngleMaps()
{
	std::vector<size_t> selFacets = worker.GetGeometry()->GetSelectedFacets();
	if (selFacets.size() == 0) {
		GLMessageBox::Display("Select at least one facet to import angle map to", "Error", GLDLG_OK, GLDLG_ICONERROR);
		return;
	}
	std::vector<FILENAME> files = GLFileBox::OpenMultipleFiles("CSV files\0*.csv\0All files\0*.*\0", "Import angle map(s)");
	if (files.size() == 0) return;
	if (selFacets.size() != files.size()) {
		GLMessageBox::Display("Select the same number of facets and files to import", "Error", GLDLG_OK, GLDLG_ICONERROR);
		return;
	}
	for (size_t i = 0; i < files.size(); i++) {
		try {
			FileReader *f = new FileReader(files[i].fullName);
			std::vector<std::vector<std::string>> table = worker.ImportCSV_string(f);
			SAFE_DELETE(f);
			AskToReset(&worker);
			worker.GetGeometry()->GetFacet(selFacets[i])->ImportAngleMap(table);
			worker.Reload();
		}
		catch (Error &e) {
				char errMsg[512];
				sprintf(errMsg, "%s\nFile:%s", e.GetMsg(), files[i].fullName);
				GLMessageBox::Display(errMsg, "Error", GLDLG_OK, GLDLG_ICONERROR);
				return;
		}
	}
}

void MolFlow::CopyAngleMapToClipboard()
{
	Geometry *geom = worker.GetGeometry();
	size_t angleMapFacetIndex;
	bool found = false;
	for (size_t i = 0; i < geom->GetNbFacet(); i++) {
		Facet* f = geom->GetFacet(i);
		if (f->selected && f->sh.anglemapParams.hasRecorded) {
			if (found) {
				GLMessageBox::Display("More than one facet with recorded angle map selected", "Error", GLDLG_OK, GLDLG_ICONERROR);
				return;
			}
			angleMapFacetIndex = i;
			found = true;
		}
	}
	if (!found) {
		GLMessageBox::Display("No facet with recorded angle map selected", "Error", GLDLG_OK, GLDLG_ICONERROR);
		return;
	}

	/*
	if (!worker.IsDpInitialized()) {
		GLMessageBox::Display("Worker Dataport not initialized yet", "Error", GLDLG_OK, GLDLG_ICONERROR);
		return;
	}*/

		try {
			std::string map = geom->GetFacet(angleMapFacetIndex)->GetAngleMap(2);
			if (map.length() > (10 * 1024 * 1024)) {
				if (GLMessageBox::Display("Angle map text over 10MB. Copy to clipboard?", "Warning", GLDLG_OK | GLDLG_CANCEL, GLDLG_ICONWARNING) != GLDLG_OK)
					return;
			}

			GLToolkit::CopyTextToClipboard(map);

		}
		catch (Error &e) {
			GLMessageBox::Display(e.GetMsg(), "Error", GLDLG_OK, GLDLG_ICONERROR);
		}
}

void MolFlow::ClearAngleMapsOnSelection() {
	//if (AskToReset()) {
		Geometry *geom = worker.GetGeometry();
		std::vector<size_t> angleMapFacetIndices;
		for (size_t i = 0; i < geom->GetNbFacet(); i++) {
			Facet* f = geom->GetFacet(i);
			if (f->selected && f->sh.anglemapParams.hasRecorded) {
				SAFE_FREE(f->angleMapCache);
				f->sh.anglemapParams.hasRecorded = false;
			}
		}
	//}
}

/*
void MolFlow::ImportDesorption_DES() {

	FILENAME *fn = GLFileBox::OpenFile(currentDir, NULL, "Import desorption File", fileDesFilters, 0);

	if (fn) {

		try {
			worker.ImportDesorption_DES(fn->fullName);
		}
		catch (Error &e) {
			char errMsg[512];
			sprintf(errMsg, "%s\nFile:%s", e.GetMsg(), fn->fullName);
			GLMessageBox::Display(errMsg, "Error", GLDLG_OK, GLDLG_ICONERROR);
		}
		worker.CalcTotalOutgassing();
		UpdateFacetParams();
	}

}
*/

void MolFlow::SaveFile() {
	if (strlen(worker.fullFileName) > 0) {

		GLProgress *progressDlg2 = new GLProgress("Saving...", "Please wait");
		progressDlg2->SetProgress(0.5);
		progressDlg2->SetVisible(true);
		//GLWindowManager::Repaint();

		try {
			worker.SaveGeometry(worker.fullFileName, progressDlg2, false);
			ResetAutoSaveTimer();
		}
		catch (Error &e) {
			char errMsg[512];
			sprintf(errMsg, "%s\nFile:%s", e.GetMsg(), worker.GetCurrentFileName());
			GLMessageBox::Display(errMsg, "Error", GLDLG_OK, GLDLG_ICONERROR);
		}
		progressDlg2->SetVisible(false);
		SAFE_DELETE(progressDlg2);
		changedSinceSave = false;

	}
	else SaveFileAs();
}

void MolFlow::LoadFile(char *fName) {

	char fullName[512];
	char shortName[512];
	strcpy(fullName, "");

	if (fName == NULL) {
		FILENAME *fn = GLFileBox::OpenFile(currentDir, NULL, "Open file", fileLFilters, 0);
		if (fn)
			strcpy(fullName, fn->fullName);
	}
	else {
		strcpy(fullName, fName);
	}

	GLProgress *progressDlg2 = new GLProgress("Preparing to load file...", "Please wait");
	progressDlg2->SetVisible(true);
	progressDlg2->SetProgress(0.0);
	//GLWindowManager::Repaint();

	if (strlen(fullName) == 0) {
		progressDlg2->SetVisible(false);
		SAFE_DELETE(progressDlg2);
		return;
	}

	char *lPart = strrchr(fullName, '\\');
	if (lPart) strcpy(shortName, lPart + 1);
	else strcpy(shortName, fullName);

	try {
		ClearFormulas();
		ClearParameters();
		ClearAllSelections();
		ClearAllViews();
		ResetSimulation(false);
		worker.LoadGeometry(fullName);
		

		Geometry *geom = worker.GetGeometry();

		
		// Default initialisation
		for (int i = 0; i < MAX_VIEWER; i++)
			viewer[i]->SetWorker(&worker);

		//UpdateModelParams();
		startSimu->SetEnabled(true);
		compACBtn->SetEnabled(modeCombo->GetSelectedIndex() == 1);
		singleACBtn->SetEnabled(modeCombo->GetSelectedIndex() == 1);
		//resetSimu->SetEnabled(true);
		ClearFacetParams();
		nbDesStart = worker.globalHitCache.globalHits.hit.nbDesorbed;
		nbHitStart = worker.globalHitCache.globalHits.hit.nbMCHit;
		AddRecent(fullName);
		geom->viewStruct = -1;

		UpdateStructMenu();
		UpdateCurrentDir(fullName);

		// Check non simple polygon
		progressDlg2->SetMessage("Checking for non simple polygons...");

		geom->CheckCollinear();
		geom->CheckNonSimple();
		geom->CheckIsolatedVertex();
		// Set up view
		// Default
		viewer[0]->SetProjection(ORTHOGRAPHIC_PROJ);
		viewer[0]->ToFrontView();
		viewer[1]->SetProjection(ORTHOGRAPHIC_PROJ);
		viewer[1]->ToTopView();
		viewer[2]->SetProjection(ORTHOGRAPHIC_PROJ);
		viewer[2]->ToSideView();
		viewer[3]->SetProjection(PERSPECTIVE_PROJ);
		viewer[3]->ToFrontView();
		SelectViewer(0);

		ResetAutoSaveTimer();
		//UpdatePlotters();

		if (timeSettings) timeSettings->RefreshMoments();
		if (momentsEditor) momentsEditor->Refresh();
		if (pressureEvolution) pressureEvolution->Reset();
		if (timewisePlotter) timewisePlotter->Refresh();
		if (histogramPlotter) histogramPlotter->Reset();
		if (histogramSettings) histogramSettings->Refresh({});
		//if (profilePlotter) profilePlotter->Refresh(); //Might have loaded views
		if (texturePlotter) texturePlotter->Update(0.0,true);
		//if (parameterEditor) parameterEditor->UpdateCombo(); //Done by ClearParameters()
		if (textureScaling) textureScaling->Update();
		if (outgassingMap) outgassingMap->Update(m_fTime, true);
		if (facetDetails) facetDetails->Update();
		if (facetCoordinates) facetCoordinates->UpdateFromSelection();
		if (vertexCoordinates) vertexCoordinates->Update();
		if (movement) movement->Update();
		if (globalSettings && globalSettings->IsVisible()) globalSettings->Update();
		if (formulaEditor) formulaEditor->Refresh();
	}
	catch (Error &e) {

		char errMsg[512];
		sprintf(errMsg, "%s\nFile:%s", e.GetMsg(), shortName);
		GLMessageBox::Display(errMsg, "Error", GLDLG_OK, GLDLG_ICONERROR);
		RemoveRecent(fName);

	}
	progressDlg2->SetVisible(false);
	SAFE_DELETE(progressDlg2);
	changedSinceSave = false;
}

void MolFlow::InsertGeometry(bool newStr, char *fName) {
	if (!AskToReset()) return;
	ResetSimulation(false);

	char fullName[1024];
	char shortName[512];
	strcpy(fullName, "");

	if (fName == NULL) {
		FILENAME *fn = GLFileBox::OpenFile(currentDir, NULL, "Open File", fileInsFilters, 0);
		if (fn)
			strcpy(fullName, fn->fullName);
	}
	else {
		strcpy(fullName, fName);
	}

	GLProgress *progressDlg2 = new GLProgress("Loading file...", "Please wait");
	progressDlg2->SetVisible(true);
	progressDlg2->SetProgress(0.0);
	//GLWindowManager::Repaint();

	if (strlen(fullName) == 0) {
		progressDlg2->SetVisible(false);
		SAFE_DELETE(progressDlg2);
		return;
	}

	char *lPart = strrchr(fullName, '\\');
	if (lPart) strcpy(shortName, lPart + 1);
	else strcpy(shortName, fullName);

	try {

		//worker.InsertGeometry(newStr, fullName);
		worker.LoadGeometry(fullName, true, newStr);

		Geometry *geom = worker.GetGeometry();
		worker.CalcTotalOutgassing();

		//Increase BB
		for (int i = 0; i < MAX_VIEWER; i++)
			viewer[i]->SetWorker(&worker);
		
		//UpdateModelParams();
		startSimu->SetEnabled(true);

		compACBtn->SetEnabled(modeCombo->GetSelectedIndex() == 1);
		singleACBtn->SetEnabled(modeCombo->GetSelectedIndex() == 1);
		//resetSimu->SetEnabled(true);
		//ClearFacetParams();
		//nbDesStart = worker.globalHitCache.globalHits.hit.nbDesorbed;
		//nbHitStart = worker.globalHitCache.globalHits.hit.nbMC;
		AddRecent(fullName);
		geom->viewStruct = -1;

		//worker.LoadTexturesGEO(fullName);
		UpdateStructMenu();
		if (profilePlotter) profilePlotter->Reset();
		if (pressureEvolution) pressureEvolution->Reset();
		if (timewisePlotter) timewisePlotter->Reset();
		//UpdateCurrentDir(fullName);

		geom->CheckCollinear();
		geom->CheckNonSimple();
		geom->CheckIsolatedVertex();

		/*
		// Set up view
		// Default
		viewer[0]->SetProjection(ORTHOGRAPHIC_PROJ);
		viewer[0]->ToFrontView();
		viewer[1]->SetProjection(ORTHOGRAPHIC_PROJ);
		viewer[1]->ToTopView();
		viewer[2]->SetProjection(ORTHOGRAPHIC_PROJ);
		viewer[2]->ToSideView();
		viewer[3]->SetProjection(PERSPECTIVE_PROJ);
		viewer[3]->ToFrontView();
		SelectViewer(0);
		*/
		
		UpdatePlotters();
		if (outgassingMap) outgassingMap->Update(m_fTime, true);
		if (facetDetails) facetDetails->Update();
		if (facetCoordinates) facetCoordinates->UpdateFromSelection();
		if (vertexCoordinates) vertexCoordinates->Update();
		if (formulaEditor) formulaEditor->Refresh();

	}
	catch (Error &e) {

		char errMsg[512];
		sprintf(errMsg, "%s\nFile:%s", e.GetMsg(), shortName);
		GLMessageBox::Display(errMsg, "Error", GLDLG_OK, GLDLG_ICONERROR);
		RemoveRecent(fName);

	}
	progressDlg2->SetVisible(false);
	SAFE_DELETE(progressDlg2);
	changedSinceSave = true;
}

void MolFlow::InsertGeometryBuffer(char *shortName) {
	if (!AskToReset()) return;
	ResetSimulation(false);

	try {

		Geometry *geom = worker.GetGeometry();
		worker.CalcTotalOutgassing();

		//Increase BB
		for (int i = 0; i < MAX_VIEWER; i++)
			viewer[i]->SetWorker(&worker);

		//UpdateModelParams();
		startSimu->SetEnabled(true);

		compACBtn->SetEnabled(modeCombo->GetSelectedIndex() == 1);
		singleACBtn->SetEnabled(modeCombo->GetSelectedIndex() == 1);
		//resetSimu->SetEnabled(true);
		//ClearFacetParams();
		//nbDesStart = worker.globalHitCache.globalHits.hit.nbDesorbed;
		//nbHitStart = worker.globalHitCache.globalHits.hit.nbMC;
		geom->viewStruct = -1;

		//worker.LoadTexturesGEO(fullName);
		UpdateStructMenu();
		if (profilePlotter) profilePlotter->Reset();
		if (pressureEvolution) pressureEvolution->Reset();
		if (timewisePlotter) timewisePlotter->Reset();
		//UpdateCurrentDir(fullName);

		geom->CheckCollinear();
		geom->CheckNonSimple();
		geom->CheckIsolatedVertex();

		/*
		// Set up view
		// Default
		viewer[0]->SetProjection(ORTHOGRAPHIC_PROJ);
		viewer[0]->ToFrontView();
		viewer[1]->SetProjection(ORTHOGRAPHIC_PROJ);
		viewer[1]->ToTopView();
		viewer[2]->SetProjection(ORTHOGRAPHIC_PROJ);
		viewer[2]->ToSideView();
		viewer[3]->SetProjection(PERSPECTIVE_PROJ);
		viewer[3]->ToFrontView();
		SelectViewer(0);
		*/

		UpdatePlotters();
		if (outgassingMap) outgassingMap->Update(m_fTime, true);
		if (facetDetails) facetDetails->Update();
		if (facetCoordinates) facetCoordinates->UpdateFromSelection();
		if (vertexCoordinates) vertexCoordinates->Update();
		if (formulaEditor) formulaEditor->Refresh();

	}
	catch (Error &e) {

		char errMsg[512];
		sprintf(errMsg, "%s\nFile:%s", e.GetMsg(), shortName);
		GLMessageBox::Display(errMsg, "Error", GLDLG_OK, GLDLG_ICONERROR);

	}
	changedSinceSave = true;
}

void MolFlow::ClearParameters() {
	worker.parameters = std::vector<Parameter>();
	if (parameterEditor) parameterEditor->UpdateCombo();
}

void MolFlow::StartStopSimulation() {
	
	if (!(worker.globalHitCache.globalHits.hit.nbMCHit > 0) && !worker.wp.calcConstantFlow && worker.moments.size() == 0) {
		bool ok = GLMessageBox::Display("Warning: in the Moments Editor, the option \"Calculate constant flow\" is disabled.\n"
			"This is useful for time-dependent simulations.\n"
			"However, you didn't define any moments, suggesting you're using steady-state mode.\n"
			"\nDo you want to continue?\n", "Strange time settings", GLDLG_OK | GLDLG_CANCEL, GLDLG_ICONWARNING) == GLDLG_OK;
		if (!ok) return;
	}
	
	worker.StartStop(m_fTime, modeCombo->GetSelectedIndex());
	if (!worker.isRunning) { //Force update on simulation stop
		UpdatePlotters();
		//if (autoUpdateFormulas) UpdateFormula();
		if (autoUpdateFormulas && formulaEditor && formulaEditor->IsVisible()) formulaEditor->ReEvaluate();
		if (particleLogger && particleLogger->IsVisible()) particleLogger->UpdateStatus();
	}

	// Frame rate measurement
	lastMeasTime = m_fTime;
	dps = 0.0;
	hps = 0.0;
	lastHps = hps;
	lastDps = dps;
	lastNbHit = worker.globalHitCache.globalHits.hit.nbMCHit;
	lastNbDes = worker.globalHitCache.globalHits.hit.nbDesorbed;
	lastUpdate = 0.0;

}

// Name: EventProc()
// Desc: Message proc function to handle key and mouse input

void MolFlow::ExportHitBufferToFile() {
	Geometry *geom = worker.GetGeometry();
	if (geom->GetNbFacet() == 0) {
		GLMessageBox::Display("Empty selection", "Error", GLDLG_OK, GLDLG_ICONERROR);
		return;
	}
	if (!worker.IsDpInitialized()) {
		GLMessageBox::Display("Worker Dataport not initialized yet", "Error", GLDLG_OK, GLDLG_ICONERROR);
		return;
	}
	FILENAME *fn = GLFileBox::SaveFile(currentDir, NULL, "Save File", fileBufferFilters, 0);
	if (fn) {
		try {
			worker.ExportHitBuffer(fn->fullName);
		}
		catch (Error &e) {
			char errMsg[512];
			sprintf(errMsg, "%s\nFile:%s", e.GetMsg(), fn->fullName);
			GLMessageBox::Display(errMsg, "Error", GLDLG_OK, GLDLG_ICONERROR);
		}
	}
	
}

void MolFlow::ExportLoadBufferToFile() {
	Geometry *geom = worker.GetGeometry();
	if (geom->GetNbFacet() == 0) {
		//GLMessageBox::Display("Empty selection", "Error", GLDLG_OK, GLDLG_ICONERROR);
		//return;
	}
	FILENAME *fn = GLFileBox::SaveFile(currentDir, NULL, "Save File", fileBufferFilters, 0);
	if (fn) {
		try {
			worker.ExportLoadBuffer(fn->fullName);
		}
		catch (Error &e) {
			char errMsg[512];
			sprintf(errMsg, "%s\nFile:%s", e.GetMsg(), fn->fullName);
			GLMessageBox::Display(errMsg, "Error", GLDLG_OK, GLDLG_ICONERROR);
		}
	}

}

void MolFlow::ImportHitBuffer(char *fName) {
	Geometry *geom = worker.GetGeometry();
	if (geom->GetNbFacet() == 0) {
		GLMessageBox::Display("No Geometry loaded. Please load corresponding file.", "Error", GLDLG_OK, GLDLG_ICONERROR);
		return;
	}
	if (!worker.IsDpInitialized()) {
		GLMessageBox::Display("Worker Dataport not initialized yet", "Error", GLDLG_OK, GLDLG_ICONERROR);
		return;
	}
	
	char fullName[512];
	char shortName[512];
	strcpy(fullName, "");

	if (fName == NULL) {
		FILENAME *fn = GLFileBox::OpenFile(currentDir, NULL, "Open file", "All files\0*.*\0", 0);
		if (fn)
			strcpy(fullName, fn->fullName);
	}
	else {
		strcpy(fullName, fName);
	}
	if (strlen(fullName) == 0) {
		return;
	}


	char *lPart = strrchr(fullName, '\\');
	if (lPart) strcpy(shortName, lPart + 1);
	else strcpy(shortName, fullName);
	try {
		worker.ImportHitBuffer(fullName);
	}
	catch (Error &e) {
		char errMsg[512];
		sprintf(errMsg, "%s\nFile:%s", e.GetMsg(), shortName);
		GLMessageBox::Display(errMsg, "Error", GLDLG_OK, GLDLG_ICONERROR);
		RemoveRecent(fName);
	}
}

void MolFlow::ImportLoadBuffer(char *fName) {

	char fullName[512];
	char shortName[512];
	strcpy(fullName, "");

	if (fName == NULL) {
		FILENAME *fn = GLFileBox::OpenFile(currentDir, NULL, "Open file", "All files\0*.*\0", 0);
		if (fn)
			strcpy(fullName, fn->fullName);
	}
	else {
		strcpy(fullName, fName);
	}

	if (strlen(fullName) == 0) {
		return;
	}

	char *lPart = strrchr(fullName, '\\');
	if (lPart) strcpy(shortName, lPart + 1);
	else strcpy(shortName, fullName);
	//try {
	//	worker.ImportLoadBuffer(fullName);
	//	int test = 0;
	//}
	//catch (Error &e) {
	//	char errMsg[512];
	//	sprintf(errMsg, "%s\nFile:%s", e.GetMsg(), shortName);
	//	GLMessageBox::Display(errMsg, "Error", GLDLG_OK, GLDLG_ICONERROR);
	//	RemoveRecent(fName);
	//}



	try {
		ClearFormulas();
		ClearParameters();
		ClearAllSelections();
		ClearAllViews();
		ResetSimulation(false);
		worker.ImportLoadBuffer(fullName);

		Geometry *geom = worker.GetGeometry();

		// Default initialisation
		for (int i = 0; i < MAX_VIEWER; i++)
			viewer[i]->SetWorker(&worker);

		//UpdateModelParams();
		startSimu->SetEnabled(true);
		compACBtn->SetEnabled(modeCombo->GetSelectedIndex() == 1);
		singleACBtn->SetEnabled(modeCombo->GetSelectedIndex() == 1);
		//resetSimu->SetEnabled(true);
		ClearFacetParams();
		nbDesStart = worker.globalHitCache.globalHits.hit.nbDesorbed;
		nbHitStart = worker.globalHitCache.globalHits.hit.nbMCHit;
		AddRecent(fullName);
		geom->viewStruct = -1;

		UpdateStructMenu();
		UpdateCurrentDir(fullName);

		// Check non simple polygon
		//geom->CheckCollinear();
		geom->CheckNonSimple();
		geom->CheckIsolatedVertex();
		
		// Set up view
		// Default
		viewer[0]->SetProjection(ORTHOGRAPHIC_PROJ);
		viewer[0]->ToFrontView();
		viewer[1]->SetProjection(ORTHOGRAPHIC_PROJ);
		viewer[1]->ToTopView();
		viewer[2]->SetProjection(ORTHOGRAPHIC_PROJ);
		viewer[2]->ToSideView();
		viewer[3]->SetProjection(PERSPECTIVE_PROJ);
		viewer[3]->ToFrontView();
		SelectViewer(0);

		ResetAutoSaveTimer();
		//UpdatePlotters();

		if (timeSettings) timeSettings->RefreshMoments();
		if (momentsEditor) momentsEditor->Refresh();
		if (pressureEvolution) pressureEvolution->Reset();
		if (timewisePlotter) timewisePlotter->Refresh();
		if (histogramPlotter) histogramPlotter->Reset();
		if (histogramSettings) histogramSettings->Refresh({});
		//if (profilePlotter) profilePlotter->Refresh(); //Might have loaded views
		if (texturePlotter) texturePlotter->Update(0.0, true);
		//if (parameterEditor) parameterEditor->UpdateCombo(); //Done by ClearParameters()
		if (textureScaling) textureScaling->Update();
		if (outgassingMap) outgassingMap->Update(m_fTime, true);
		if (facetDetails) facetDetails->Update();
		if (facetCoordinates) facetCoordinates->UpdateFromSelection();
		if (vertexCoordinates) vertexCoordinates->Update();
		if (movement) movement->Update();
		if (globalSettings && globalSettings->IsVisible()) globalSettings->Update();
		if (formulaEditor) formulaEditor->Refresh();
	}
	catch (Error &e) {

		char errMsg[512];
		sprintf(errMsg, "%s\nFile:%s", e.GetMsg(), shortName);
		GLMessageBox::Display(errMsg, "Error", GLDLG_OK, GLDLG_ICONERROR);
		RemoveRecent(fName);

	}
	changedSinceSave = false;
}

void MolFlow::ProcessMessage(GLComponent *src, int message)
{

	if (ProcessMessage_shared(src, message))
	{
		int test = 0;
		return;
	} //Already processed by common interface
	Geometry *geom = worker.GetGeometry();
	switch (message) {

		//MENU --------------------------------------------------------------------
	case MSG_MENU:
		switch (src->GetId()) {

		case MENU_FILE_IMPORTDES_SYN:

			if (geom->IsLoaded()) {
				if (!importDesorption) importDesorption = new ImportDesorption();
				importDesorption->SetGeometry(geom, &worker);
				importDesorption->SetVisible(true);
			}
			else GLMessageBox::Display("No geometry loaded.", "No geometry", GLDLG_OK, GLDLG_ICONERROR);
			break;
		/* //Deprecated
		case MENU_FILE_IMPORTDES_DES:
			ImportDesorption_DES();
			break;
		*/

		case MENU_FILE_EXPORTTEXTURE_AREA:
			ExportTextures(0, 0); break;
		case MENU_FILE_EXPORTTEXTURE_MCHITS:
			ExportTextures(0, 1); break;
		case MENU_FILE_EXPORTTEXTURE_IMPINGEMENT:
			ExportTextures(0, 2); break;
		case MENU_FILE_EXPORTTEXTURE_PART_DENSITY:
			ExportTextures(0, 3); break;
		case MENU_FILE_EXPORTTEXTURE_GAS_DENSITY:
			ExportTextures(0, 4); break;
		case MENU_FILE_EXPORTTEXTURE_PRESSURE:
			ExportTextures(0, 5); break;
		case MENU_FILE_EXPORTTEXTURE_AVG_V:
			ExportTextures(0, 6); break;
		case MENU_FILE_EXPORTTEXTURE_V_VECTOR:
			ExportTextures(1, 7); break;
		case MENU_FILE_EXPORTTEXTURE_N_VECTORS:
			ExportTextures(0, 8); break;
		case MENU_FILE_EXPORTTEXTURE_AREA_COORD:
			ExportTextures(1, 0); break;
		case MENU_FILE_EXPORTTEXTURE_MCHITS_COORD:
			ExportTextures(1, 1); break;
		case MENU_FILE_EXPORTTEXTURE_IMPINGEMENT_COORD:
			ExportTextures(1, 2); break;
		case MENU_FILE_EXPORTTEXTURE_PART_DENSITY_COORD:
			ExportTextures(1, 3); break;
		case MENU_FILE_EXPORTTEXTURE_GAS_DENSITY_COORD:
			ExportTextures(1, 4); break;
		case MENU_FILE_EXPORTTEXTURE_PRESSURE_COORD:
			ExportTextures(1, 5); break;
		case MENU_FILE_EXPORTTEXTURE_AVG_V_COORD:
			ExportTextures(1, 6); break;
		case MENU_FILE_EXPORTTEXTURE_V_VECTOR_COORD:
			ExportTextures(1, 7); break;
		case MENU_FILE_EXPORTTEXTURE_N_VECTORS_COORD:
			ExportTextures(1, 8); break;

		case MENU_FILE_EXPORTPROFILES:
			ExportProfiles();
			break;
		case MENU_FILE_EXPORTBUFFER_LOAD:
			ExportLoadBufferToFile();
			break;
		case MENU_FILE_EXPORTBUFFER_HIT:
			ExportHitBufferToFile();
			break;
		case MENU_FILE_IMPORTBUFFER_HIT:
			//ImportHitBuffer(char *fName);
			break;
		case MENU_TOOLS_MOVINGPARTS:
			if (!movement) movement = new Movement(geom, &worker);
			movement->Update();
			movement->SetVisible(true);
			break;

		case MENU_EDIT_TSCALING:
			if (!textureScaling || !textureScaling->IsVisible()) {
				SAFE_DELETE(textureScaling);
				textureScaling = new TextureScaling();
				textureScaling->Display(&worker, viewer);
			}
			break;

		case MENU_EDIT_GLOBALSETTINGS:
			if (!globalSettings) globalSettings = new GlobalSettings(&worker);
			globalSettings->Update();
			globalSettings->SetVisible(true);
			break;

		case MENU_TOOLS_PROFPLOTTER:
			if (!profilePlotter) profilePlotter = new ProfilePlotter();
			profilePlotter->Display(&worker);
			break;
		case MENU_TOOLS_TEXPLOTTER:
			if (!texturePlotter) texturePlotter = new TexturePlotter();
			texturePlotter->Display(&worker);
			break;
		case MENU_FACET_DETAILS:
			if (facetDetails == NULL) facetDetails = new FacetDetails();
			facetDetails->Display(&worker);
			break;

		case MENU_TIME_PRESSUREEVOLUTION:
			if (!pressureEvolution) pressureEvolution = new PressureEvolution(&worker);
			pressureEvolution->SetVisible(true);
			break;

		case MENU_FACET_OUTGASSINGMAP:
			if (!outgassingMap) outgassingMap = new OutgassingMap();
			outgassingMap->Display(&worker);
			break;
		case MENU_FACET_REMOVESEL:
		{
			auto selectedFacets = geom->GetSelectedFacets();
			if (selectedFacets.size() == 0) return; //Nothing selected
			if (GLMessageBox::Display("Remove selected facets?", "Question", GLDLG_OK | GLDLG_CANCEL, GLDLG_ICONINFO) == GLDLG_OK) {
				if (AskToReset()) {
					if (worker.isRunning) worker.Stop_Public();
					geom->RemoveFacets(selectedFacets);
					worker.CalcTotalOutgassing();
					//geom->CheckIsolatedVertex();
					UpdateModelParams();
					UpdatePlotters();
					if (vertexCoordinates) vertexCoordinates->Update();
					if (facetCoordinates) facetCoordinates->UpdateFromSelection();
					// Send to sub process
					worker.Reload();
				}
			}
			break;
		}
		case MENU_FACET_SELECTSTICK:
			geom->UnselectAll();
			for (int i = 0; i < geom->GetNbFacet(); i++)
				if (geom->GetFacet(i)->sh.sticking != 0.0 && !geom->GetFacet(i)->IsTXTLinkFacet())
					geom->SelectFacet(i);
			geom->UpdateSelection();
			UpdateFacetParams(true);
			break;

		case MENU_FACET_SELECTREFL:
			geom->UnselectAll();
			for (int i = 0; i < geom->GetNbFacet(); i++) {
				Facet *f = geom->GetFacet(i);
				if (f->sh.desorbType == DES_NONE && f->sh.sticking == 0.0 && f->sh.opacity > 0.0)
					geom->SelectFacet(i);
			}
			geom->UpdateSelection();
			UpdateFacetParams(true);
			break;

		case MENU_FACET_SELECTVOL:
			geom->UnselectAll();
			for (int i = 0; i < geom->GetNbFacet(); i++)
				if (geom->GetFacet(i)->sh.isVolatile)
					geom->SelectFacet(i);
			geom->UpdateSelection();
			UpdateFacetParams(true);
			break;

		case MENU_FACET_SELECTDES:
			geom->UnselectAll();
			for (int i = 0; i < geom->GetNbFacet(); i++)
				if (geom->GetFacet(i)->sh.desorbType != DES_NONE)
					geom->SelectFacet(i);
			geom->UpdateSelection();
			UpdateFacetParams(true);
			break;
		case MENU_SELECT_HASDESFILE:
			geom->UnselectAll();
			for (int i = 0; i < geom->GetNbFacet(); i++)
				if (geom->GetFacet(i)->hasOutgassingFile)
					geom->SelectFacet(i);
			geom->UpdateSelection();
			UpdateFacetParams(true);
			break;

		case MENU_TIME_SETTINGS:
			if (!timeSettings) timeSettings = new TimeSettings(&worker);
			timeSettings->SetVisible(true);
			break;
		case MENU_TIME_MOMENTS_EDITOR:
			if (momentsEditor == NULL || !momentsEditor->IsVisible()) {
				SAFE_DELETE(momentsEditor);
				momentsEditor = new MomentsEditor(&worker);
				momentsEditor->Refresh();
				momentsEditor->SetVisible(true);
			}
			break;
		case MENU_TIME_PARAMETER_EDITOR:
			if (parameterEditor == NULL) parameterEditor = new ParameterEditor(&worker);
			parameterEditor->UpdateCombo();
			parameterEditor->SetVisible(true);
			break;
		case MENU_TIMEWISE_PLOTTER:
			if (!timewisePlotter) timewisePlotter = new TimewisePlotter();
			timewisePlotter->Display(&worker);
			break;
		case MENU_VERTEX_REMOVE:
			if (geom->IsLoaded()) {
				if (GLMessageBox::Display("Remove Selected vertices?\nNote: It will also affect facets that contain them!", "Question", GLDLG_OK | GLDLG_CANCEL, GLDLG_ICONINFO) == GLDLG_OK) {
					if (AskToReset()) {
						if (worker.isRunning) worker.Stop_Public();
						geom->RemoveSelectedVertex();
						worker.CalcTotalOutgassing();
						geom->Rebuild(); //Will recalculate facet parameters
						UpdateModelParams();
						if (vertexCoordinates) vertexCoordinates->Update();
						if (facetCoordinates) facetCoordinates->UpdateFromSelection();
						if (profilePlotter) profilePlotter->Refresh();
						if (pressureEvolution) pressureEvolution->Refresh();
						if (timewisePlotter) timewisePlotter->Refresh();
						if (facetCoordinates) facetCoordinates->UpdateFromSelection();
						if (vertexCoordinates) vertexCoordinates->Update();
						// Send to sub process
						try { worker.Reload(); }
						catch (Error &e) {
							GLMessageBox::Display((char *)e.GetMsg(), "Error reloading worker", GLDLG_OK, GLDLG_ICONERROR);
						}
					}
				}

			}
			else GLMessageBox::Display("No geometry loaded.", "No geometry", GLDLG_OK, GLDLG_ICONERROR);
			break;
		}

		//TEXT --------------------------------------------------------------------
	case MSG_TEXT_UPD:
		/*if (src == facetSticking) {
			//calcFlow();
			facetApplyBtn->SetEnabled(true);
		}
		else if*/
		if (src == facetOpacity) {
			facetApplyBtn->SetEnabled(true);
		}
		else if (src == facetDesTypeN) {
			facetApplyBtn->SetEnabled(true);
		}
		else if (src == facetTemperature) {
			//calcFlow();
			if (facetcoverage->GetText() != "...") {
				//calcStickingnew();
				calcDesorptionRate();
			}
			facetApplyBtn->SetEnabled(true);
		}
		else if (src == facetFlow) {
			double outgassing;
			double area;
			facetFlow->GetNumber(&outgassing);
			facetArea->GetNumber(&area);
			if (area == 0) facetFlowArea->SetText("#DIV0");
			else facetFlowArea->SetText(outgassing / area);
			facetApplyBtn->SetEnabled(true);
			facetFILabel->SetState(true);
			facetFIAreaLabel->SetState(false);
			// //else if ( src == facetMass ) {
			//  facetApplyBtn->SetEnabled(true);
		}
		else if (src == facetFlowArea) {
			double flowPerArea;
			double area;
			facetFlowArea->GetNumber(&flowPerArea);
			facetArea->GetNumber(&area);
			facetFlow->SetText(flowPerArea*area);
			facetApplyBtn->SetEnabled(true);
			facetFIAreaLabel->SetState(true);
			facetFILabel->SetState(false);
			// //else if ( src == facetMass ) {
			//  facetApplyBtn->SetEnabled(true);
		}
		/*else if (src == facetPumping) {// do nothing, we don't need pumping here
			//calcSticking();
			//facetApplyBtn->SetEnabled(true);
			
			// //else if ( src == facetMass ) {
			//  facetApplyBtn->SetEnabled(true);
		}*/
		else if (src == facetcoverage) {
			calcCovering();
			//calcStickingnew();
			calcDesorptionRate();
			facetApplyBtn->SetEnabled(true);
		}
		else if (src == facetcovering) {
			calcCoverage();
			//calcStickingnew();
			calcDesorptionRate();
			facetApplyBtn->SetEnabled(true);
		}
		else if (src == facetSurfaceFactor) {
			facetApplyBtn->SetEnabled(true);
		}
		else if (src == facetDepth) {
			calcVolume();
			calcMass();
			facetApplyBtn->SetEnabled(true);
		}
		else if (src == facetVolume) {
			calcDepth();
			calcMass();
			facetApplyBtn->SetEnabled(true);
		}
		else if (src == facetDiffusion) {
			facetApplyBtn->SetEnabled(true);
		}

		else if (src == facetConcentration) {
			calcMass();
			facetApplyBtn->SetEnabled(true);
		}
		else if (src == facetMass) {
			calcConcentration();
			facetApplyBtn->SetEnabled(true);
		}

		break;

	case MSG_TEXT:
		/*if (src == facetSticking) {
			ApplyFacetParams();
		}
		else if*/
		if (src == facetDesTypeN) {
			ApplyFacetParams();
		}
		else if (src == facetOpacity) {
			ApplyFacetParams();
		}
		else if (src == facetTemperature) {
			ApplyFacetParams();
		}
		/*else if (src == facetPumping) {
			ApplyFacetParams();
		}*/
		else if (src == facetFlow) {
			ApplyFacetParams();
		}
		else if (src == facetFlowArea) {
			ApplyFacetParams();
		}
		else if (src == facetcoverage) {
			ApplyFacetParams();
		}
		else if (src == facetcovering) {
			ApplyFacetParams();
		}
		else if (src == facetSurfaceFactor) {
			ApplyFacetParams();
		}
		else if (src == facetDepth) {
			ApplyFacetParams();
		}
		break;

		//COMBO -------------------------------------------------------------------
	case MSG_COMBO:

		if (src == facetDesType) {
			facetApplyBtn->SetEnabled(true);
			bool hasDesorption = !facetDesType->GetSelectedIndex() == 0;

			facetFlow->SetEditable(hasDesorption);
			facetFlowArea->SetEditable(hasDesorption);

			int color = (hasDesorption) ? 0 : 110;
			facetFILabel->SetTextColor(color, color, color);
			facetFIAreaLabel->SetTextColor(color, color, color);
			facetFILabel->SetEnabled(hasDesorption);
			facetFIAreaLabel->SetEnabled(hasDesorption);
			facetDesTypeN->SetEditable(facetDesType->GetSelectedIndex() == 3);
		}
		else if (src == facetRecType) {
			facetApplyBtn->SetEnabled(true);
		}
		else if (src == facetSideType) {
			facetApplyBtn->SetEnabled(true);
		}
		else if (src == modeCombo) {

			compACBtn->SetEnabled(modeCombo->GetSelectedIndex() == 1);

			singleACBtn->SetEnabled(modeCombo->GetSelectedIndex() == 1);
			UpdateFacetHits();
		}
		break;

		//TOGGLE ------------------------------------------------------------------
	case MSG_TOGGLE:
		if (src == facetFILabel) {
			facetFILabel->SetState(true);
			facetFIAreaLabel->SetState(false);
		}
		else if (src == facetFIAreaLabel) {
			facetFILabel->SetState(false);
			facetFIAreaLabel->SetState(true);
		}
		else if (src == autoFrameMoveToggle) {
			autoFrameMove = autoFrameMoveToggle->GetState();
			forceFrameMoveButton->SetEnabled(!autoFrameMove);
		}
		else {
			// Update viewer flags
			UpdateViewerFlags();
		}
		break;

		//BUTTON ------------------------------------------------------------------
	case MSG_BUTTON:
		if (src == startSimu) {
			if(!history) history = new HistoryWin(&worker);
			changedSinceSave = true;
			StartStopSimulation();
			resetSimu->SetEnabled(!worker.isRunning);
			if (!worker.isRunning) {
					history->UpdateList();
			}
		}

		else if (src == facetApplyBtn) {
			//changedSinceSave=true;
			ApplyFacetParams();
		}
		else if (src == facetDetailsBtn) {
			if (facetDetails == NULL) facetDetails = new FacetDetails();
			facetDetails->Display(&worker);
		}
		else if (src == viewerMoreButton) {
			if (!viewer3DSettings)	viewer3DSettings = new Viewer3DSettings();
			viewer3DSettings->SetVisible(!viewer3DSettings->IsVisible());
			viewer3DSettings->Reposition();
			viewer3DSettings->Refresh(geom, viewer[curViewer]);

		}
		else if (src == textureScalingBtn) {
			if (!textureScaling || !textureScaling->IsVisible()) {
				SAFE_DELETE(textureScaling);
				textureScaling = new TextureScaling();
				textureScaling->Display(&worker, viewer);
			}
			else {
				textureScaling->SetVisible(false);
			}
		}
		else if (src == profilePlotterBtn) {
			if (!profilePlotter) profilePlotter = new ProfilePlotter();
			if (!profilePlotter->IsVisible()) profilePlotter->Display(&worker);
			else profilePlotter->SetVisible(false);
		}
		else if (src == texturePlotterBtn) {
			if (!texturePlotter) texturePlotter = new TexturePlotter();
			if (!texturePlotter->IsVisible()) texturePlotter->Display(&worker);
			else {
				texturePlotter->SetVisible(false);
				SAFE_DELETE(texturePlotter);
			}
		}
		else if (src == globalSettingsBtn) {
			if (!globalSettings) globalSettings = new GlobalSettings(&worker);
			if (!globalSettings->IsVisible()) {
				globalSettings->Update();
				globalSettings->SetVisible(true);
			}
			else globalSettings->SetVisible(false);
		}

		else if (src == historyBtn) {
			if (!history) history = new HistoryWin(&worker);
			if (!history->IsVisible()) {
				history->SetVisible(true);
			}
			else history->SetVisible(false);
		}

		else if (src == facetAdvParamsBtn) {
			if (!facetAdvParams) {
				facetAdvParams = new FacetAdvParams(&worker);
				facetAdvParams->Refresh(geom->GetSelectedFacets());
			}
			facetAdvParams->SetVisible(!facetAdvParams->IsVisible());
			facetAdvParams->Reposition();
		}

		else if (src == compACBtn) {
			try { lastUpdate = 0.0; worker.ComputeAC(m_fTime); }
			catch (Error &e) {
				GLMessageBox::Display((char *)e.GetMsg(), "Error", GLDLG_OK, GLDLG_ICONERROR);
				return;
			}
			break;
		}
		else if (src == singleACBtn) {
			try { lastUpdate = 0.0; worker.StepAC(m_fTime); }
			catch (Error &e) {
				GLMessageBox::Display((char *)e.GetMsg(), "Error", GLDLG_OK, GLDLG_ICONERROR);
				return;
			}
			break;
		}
		/*else {
			ProcessFormulaButtons(src);
		}*/
		break;
	}
}

void MolFlow::BuildPipe(double ratio, int steps) {

	char tmp[256];
	MolflowGeometry *geom = worker.GetMolflowGeometry();

	double R = 1.0;
	double L = ratio * R;
	int    step;
	if (steps) step = 5; //Quick Pipe
	else {
		sprintf(tmp, "100");
		//sprintf(title,"Pipe L/R = %g",L/R);
		char *nbF = GLInputBox::GetInput(tmp, "Number of facet", "Build Pipe");
		if (!nbF) return;
		if ((sscanf(nbF, "%d", &step) <= 0) || (step < 3)) {
			GLMessageBox::Display("Invalid number", "Error", GLDLG_OK, GLDLG_ICONERROR);
			return;
		}
	}
	std::ostringstream temp;
	temp << "PIPE" << L / R;
	geom->UpdateName(temp.str().c_str());
	ResetSimulation(false);

	try {
		geom->BuildPipe(L, R, 0, step);
		worker.CalcTotalOutgassing();
		//default values
		worker.wp.enableDecay = false;
		worker.wp.halfLife = 1;
		worker.wp.gasMass = 28;
		worker.ResetMoments();
	}
	catch (Error &e) {
		GLMessageBox::Display((char *)e.GetMsg(), "Error building pipe", GLDLG_OK, GLDLG_ICONERROR);
		geom->Clear();
		return;
	}
	//worker.globalHitCache.globalHits.hit.nbDesorbed = 0; //Already done by ResetWorkerStats
	//sprintf(tmp,"L|R %g",L/R);
	worker.SetCurrentFileName("");
	nbDesStart = 0;
	nbHitStart = 0;
	
	for (int i = 0; i < MAX_VIEWER; i++)
		viewer[i]->SetWorker(&worker);
	
	//UpdateModelParams();
	startSimu->SetEnabled(true);
	compACBtn->SetEnabled(modeCombo->GetSelectedIndex() == 1);
	//resetSimu->SetEnabled(true);
	ClearFacetParams();
	ClearFormulas();
	ClearParameters();
	ClearAllSelections();
	ClearAllViews();

	/*GLParser *f = new GLParser();
	f->SetExpression("A2/SUMDES");
	f->SetName("Trans. Prob.");
	f->Parse();*/
	AddFormula("Trans.prob.", "A2/SUMDES");

	UpdateStructMenu();
	// Send to sub process
	worker.Reload();

	//UpdatePlotters();

	if (timeSettings) timeSettings->RefreshMoments();
	if (momentsEditor) momentsEditor->Refresh();
	if (pressureEvolution) pressureEvolution->Reset();
	if (timewisePlotter) timewisePlotter->Refresh();
	if (profilePlotter) profilePlotter->Refresh();
	if (histogramSettings) histogramSettings->Refresh({});
	if (histogramPlotter) histogramPlotter->Reset();
	if (texturePlotter) texturePlotter->Update(0.0, true);
	//if (parameterEditor) parameterEditor->UpdateCombo(); //Done by ClearParameters()
	if (textureScaling) textureScaling->Update();
	if (outgassingMap) outgassingMap->Update(m_fTime, true);
	if (facetDetails) facetDetails->Update();
	if (facetCoordinates) facetCoordinates->UpdateFromSelection();
	if (vertexCoordinates) vertexCoordinates->Update();
	if (movement) movement->Update();
	if (globalSettings && globalSettings->IsVisible()) globalSettings->Update();
	if (formulaEditor) formulaEditor->Refresh();
	UpdateTitle();
	changedSinceSave = false;
	ResetAutoSaveTimer();
}

void MolFlow::EmptyGeometry() {

	Geometry *geom = worker.GetGeometry();
	ResetSimulation(false);

	try {
		_ASSERTE(_CrtCheckMemory());
		geom->EmptyGeometry();
		_ASSERTE(_CrtCheckMemory());
		worker.CalcTotalOutgassing();
		//default values
		worker.wp.enableDecay = false;
		worker.wp.halfLife = 1;
		worker.wp.gasMass = 28;
		worker.ResetMoments();
	}
	catch (Error &e) {
		GLMessageBox::Display((char *)e.GetMsg(), "Error resetting geometry", GLDLG_OK, GLDLG_ICONERROR);
		geom->Clear();
		return;
	}
	worker.SetCurrentFileName("");
	nbDesStart = 0;
	nbHitStart = 0;

	for (int i = 0; i < MAX_VIEWER; i++)
		viewer[i]->SetWorker(&worker);

	//UpdateModelParams();
	startSimu->SetEnabled(true);
	compACBtn->SetEnabled(modeCombo->GetSelectedIndex() == 1);
	//resetSimu->SetEnabled(true);
	ClearFacetParams();
	ClearFormulas();
	ClearParameters();
	ClearAllSelections();
	ClearAllViews();

	/*GLParser *f = new GLParser();
	f->SetExpression("A2/SUMDES");
	f->SetName("Trans. Prob.");
	f->Parse();*/
	//AddFormula("Trans.prob.", "A2/SUMDES");

	UpdateStructMenu();
	// Send to sub process
	worker.Reload();

	//UpdatePlotters();

	if (timeSettings) timeSettings->RefreshMoments();
	if (momentsEditor) momentsEditor->Refresh();
	if (pressureEvolution) pressureEvolution->Refresh();
	if (timewisePlotter) timewisePlotter->Refresh();
	if (profilePlotter) profilePlotter->Refresh();
	if (histogramSettings) histogramSettings->Refresh({});
	if (histogramPlotter) histogramPlotter->Reset();
	if (texturePlotter) texturePlotter->Update(0.0, true);
	//if (parameterEditor) parameterEditor->UpdateCombo(); //Done by ClearParameters()
	if (outgassingMap) outgassingMap->Update(m_fTime, true);
	if (movement) movement->Update();
	if (globalSettings && globalSettings->IsVisible()) globalSettings->Update();
	if (formulaEditor) formulaEditor->Refresh();
	
	if (textureScaling) textureScaling->Update();
	if (facetDetails) facetDetails->Update();
	if (facetCoordinates) facetCoordinates->UpdateFromSelection();
	if (vertexCoordinates) vertexCoordinates->Update();
	
	UpdateTitle();
	changedSinceSave = false;
	ResetAutoSaveTimer();
	UpdatePlotters();
}

void MolFlow::LoadConfig() {

	FileReader *f = NULL;
	char *w;
	nbRecent = 0;

	try {

		f = new FileReader("molflow.cfg");
		MolflowGeometry *geom = worker.GetMolflowGeometry();

		f->ReadKeyword("showRules"); f->ReadKeyword(":");
		for (int i = 0; i < MAX_VIEWER; i++)
			viewer[i]->showRule = f->ReadInt();
		f->ReadKeyword("showNormals"); f->ReadKeyword(":");
		for (int i = 0; i < MAX_VIEWER; i++)
			viewer[i]->showNormal = f->ReadInt();
		f->ReadKeyword("showUV"); f->ReadKeyword(":");
		for (int i = 0; i < MAX_VIEWER; i++)
			viewer[i]->showUV = f->ReadInt();
		f->ReadKeyword("showLines"); f->ReadKeyword(":");
		for (int i = 0; i < MAX_VIEWER; i++)
			viewer[i]->showLine = f->ReadInt();
		f->ReadKeyword("showLeaks"); f->ReadKeyword(":");
		for (int i = 0; i < MAX_VIEWER; i++)
			viewer[i]->showLeak = f->ReadInt();
		f->ReadKeyword("showHits"); f->ReadKeyword(":");
		for (int i = 0; i < MAX_VIEWER; i++)
			viewer[i]->showHit = f->ReadInt();
		f->ReadKeyword("showVolume"); f->ReadKeyword(":");
		for (int i = 0; i < MAX_VIEWER; i++)
			viewer[i]->showVolume = f->ReadInt();
		f->ReadKeyword("showTexture"); f->ReadKeyword(":");
		for (int i = 0; i < MAX_VIEWER; i++)
			viewer[i]->showTexture = f->ReadInt();
		f->ReadKeyword("showFilter"); f->ReadKeyword(":");
		for (int i = 0; i < MAX_VIEWER; i++)
			viewer[i]->showFilter = f->ReadInt();
		f->ReadKeyword("showIndices"); f->ReadKeyword(":");
		for (int i = 0; i < MAX_VIEWER; i++)
			viewer[i]->showIndex = f->ReadInt();
		f->ReadKeyword("showVertices"); f->ReadKeyword(":");
		for (int i = 0; i < MAX_VIEWER; i++)
			viewer[i]->showVertex = f->ReadInt();
		f->ReadKeyword("showMode"); f->ReadKeyword(":");
		for (int i = 0; i < MAX_VIEWER; i++)
			viewer[i]->showBack = f->ReadInt();
		f->ReadKeyword("showMesh"); f->ReadKeyword(":");
		for (int i = 0; i < MAX_VIEWER; i++)
			viewer[i]->showMesh = f->ReadInt();
		f->ReadKeyword("showHidden"); f->ReadKeyword(":");
		for (int i = 0; i < MAX_VIEWER; i++)
			viewer[i]->showHidden = f->ReadInt();
		f->ReadKeyword("showHiddenVertex"); f->ReadKeyword(":");
		for (int i = 0; i < MAX_VIEWER; i++)
			viewer[i]->showHiddenVertex = f->ReadInt();
		f->ReadKeyword("showTimeOverlay"); f->ReadKeyword(":");
		for (int i = 0; i < MAX_VIEWER; i++)
			viewer[i]->showTime = f->ReadInt();
		f->ReadKeyword("texColormap"); f->ReadKeyword(":");
		for (int i = 0; i < MAX_VIEWER; i++)
			//viewer[i]->showColormap = 
			f->ReadInt();
		f->ReadKeyword("translation"); f->ReadKeyword(":");
		for (int i = 0; i < MAX_VIEWER; i++)
			viewer[i]->transStep = f->ReadDouble();
		f->ReadKeyword("dispNumLines"); f->ReadKeyword(":");
		for (int i = 0; i < MAX_VIEWER; i++)
			viewer[i]->dispNumHits = f->ReadLLong();
		f->ReadKeyword("dispNumLeaks"); f->ReadKeyword(":");
		for (int i = 0; i < MAX_VIEWER; i++)
			viewer[i]->dispNumLeaks = f->ReadLLong();
		f->ReadKeyword("dirShow"); f->ReadKeyword(":");
		for (int i = 0; i < MAX_VIEWER; i++)
			viewer[i]->showDir = f->ReadInt();
		f->ReadKeyword("dirNorme"); f->ReadKeyword(":");
		geom->SetNormeRatio((float)f->ReadDouble());
		f->ReadKeyword("dirAutoNormalize"); f->ReadKeyword(":");
		geom->SetAutoNorme(f->ReadInt());
		f->ReadKeyword("dirCenter"); f->ReadKeyword(":");
		geom->SetCenterNorme(f->ReadInt());
		f->ReadKeyword("angle"); f->ReadKeyword(":");
		for (int i = 0; i < MAX_VIEWER; i++)
			viewer[i]->angleStep = f->ReadDouble();
		f->ReadKeyword("autoScale"); f->ReadKeyword(":");
		geom->texAutoScale = f->ReadInt();
		f->ReadKeyword("autoScale_include_constant_flow"); f->ReadKeyword(":");
		geom->texAutoScaleIncludeConstantFlow = f->ReadInt();

		f->ReadKeyword("textures_min_pressure_all"); f->ReadKeyword(":");
		geom->texture_limits[0].autoscale.min.all = f->ReadDouble();
		f->ReadKeyword("textures_min_pressure_moments_only"); f->ReadKeyword(":");
		geom->texture_limits[0].autoscale.min.moments_only = f->ReadDouble();
		f->ReadKeyword("textures_max_pressure_all"); f->ReadKeyword(":");
		geom->texture_limits[0].autoscale.max.all = f->ReadDouble();
		f->ReadKeyword("textures_max_pressure_moments_only"); f->ReadKeyword(":");
		geom->texture_limits[0].autoscale.max.moments_only = f->ReadDouble();

		f->ReadKeyword("textures_min_impingement_all"); f->ReadKeyword(":");
		geom->texture_limits[1].autoscale.min.all = f->ReadDouble();
		f->ReadKeyword("textures_min_impingement_moments_only"); f->ReadKeyword(":");
		geom->texture_limits[1].autoscale.min.moments_only = f->ReadDouble();
		f->ReadKeyword("textures_max_impingement_all"); f->ReadKeyword(":");
		geom->texture_limits[1].autoscale.max.all = f->ReadDouble();
		f->ReadKeyword("textures_max_impingement_moments_only"); f->ReadKeyword(":");
		geom->texture_limits[1].autoscale.max.moments_only = f->ReadDouble();

		f->ReadKeyword("textures_min_density_all"); f->ReadKeyword(":");
		geom->texture_limits[2].autoscale.min.all = f->ReadDouble();
		f->ReadKeyword("textures_min_density_moments_only"); f->ReadKeyword(":");
		geom->texture_limits[2].autoscale.min.moments_only = f->ReadDouble();
		f->ReadKeyword("textures_max_density_all"); f->ReadKeyword(":");
		geom->texture_limits[2].autoscale.max.all = f->ReadDouble();
		f->ReadKeyword("textures_max_density_moments_only"); f->ReadKeyword(":");
		geom->texture_limits[2].autoscale.max.moments_only = f->ReadDouble();

		f->ReadKeyword("processNum"); f->ReadKeyword(":");
		nbProc = f->ReadLLong();
#ifdef _DEBUG
		nbProc = 1;
#endif
		if (nbProc <= 0) nbProc = 1;
		f->ReadKeyword("recents"); f->ReadKeyword(":"); f->ReadKeyword("{");
		w = f->ReadString();
		while (strcmp(w, "}") != 0 && nbRecent < MAX_RECENT) {
			recents[nbRecent] = _strdup(w);
			nbRecent++;
			w = f->ReadString();
		}

		f->ReadKeyword("cdir"); f->ReadKeyword(":");
		strcpy(currentDir, f->ReadString());
		f->ReadKeyword("cseldir"); f->ReadKeyword(":");
		strcpy(currentSelDir, f->ReadString());
		f->ReadKeyword("autonorme"); f->ReadKeyword(":");
		geom->SetAutoNorme(f->ReadInt());
		f->ReadKeyword("centernorme"); f->ReadKeyword(":");
		geom->SetCenterNorme(f->ReadInt());
		f->ReadKeyword("normeratio"); f->ReadKeyword(":");
		geom->SetNormeRatio((float)(f->ReadDouble()));
		f->ReadKeyword("autoSaveFrequency"); f->ReadKeyword(":");
		autoSaveFrequency = f->ReadDouble();
		f->ReadKeyword("autoSaveSimuOnly"); f->ReadKeyword(":");
		autoSaveSimuOnly = f->ReadInt();
		f->ReadKeyword("checkForUpdates"); f->ReadKeyword(":");
		/*checkForUpdates =*/ f->ReadInt(); //Old checkforupdates
		f->ReadKeyword("autoUpdateFormulas"); f->ReadKeyword(":");
		autoUpdateFormulas = f->ReadInt();
		f->ReadKeyword("compressSavedFiles"); f->ReadKeyword(":");
		compressSavedFiles = f->ReadInt();
		f->ReadKeyword("sh.gasMass"); f->ReadKeyword(":");
		worker.wp.gasMass = f->ReadDouble();
		f->ReadKeyword("expandShortcutPanel"); f->ReadKeyword(":");
		bool isOpen = f->ReadInt();
		if (isOpen) shortcutPanel->Open();
		else shortcutPanel->Close();
		f->ReadKeyword("hideLot"); f->ReadKeyword(":");
		for (int i = 0; i < MAX_VIEWER; i++)
			viewer[i]->hideLot = f->ReadInt();
		f->ReadKeyword("lowFluxMode"); f->ReadKeyword(":");
		worker.ontheflyParams.lowFluxMode = f->ReadInt();
		f->ReadKeyword("lowFluxCutoff"); f->ReadKeyword(":");
		worker.ontheflyParams.lowFluxCutoff = f->ReadDouble();
		f->ReadKeyword("leftHandedView"); f->ReadKeyword(":");
		leftHandedView = f->ReadInt();
	}
	catch (...) {
		/*std::ostringstream tmp;
		tmp << err.GetMsg() << "\n\nThis is normal on the first launch and if you upgrade from an earlier version\n";
		tmp << "MolFlow will use default program settings.\nWhen you quit, a correct config file will be written\n";
		GLMessageBox::Display(tmp.str().c_str(), "Error loading config file", GLDLG_OK, GLDLG_ICONINFO);*/
		SAFE_DELETE(f);
		return;
	}
	SAFE_DELETE(f);
}


#define WRITEI(name,var) {             \
	f->Write(name);                      \
	f->Write(":");                       \
	for(size_t i=0;i<MAX_VIEWER;i++)        \
	f->Write(viewer[i]->var," ");   \
	f->Write("\n");                      \
}


#define WRITED(name,var) {             \
	f->Write(name);                      \
	f->Write(":");                       \
	for(size_t i=0;i<MAX_VIEWER;i++)        \
	f->Write(viewer[i]->var," ");\
	f->Write("\n");                      \
}


void MolFlow::SaveConfig() {

	FileWriter *f = NULL;

	try {

		f = new FileWriter("molflow.cfg");
		MolflowGeometry *geom = worker.GetMolflowGeometry();

		// Save flags
		WRITEI("showRules", showRule);
		WRITEI("showNormals", showNormal);
		WRITEI("showUV", showUV);
		WRITEI("showLines", showLine);
		WRITEI("showLeaks", showLeak);
		WRITEI("showHits", showHit);
		WRITEI("showVolume", showVolume);
		WRITEI("showTexture", showTexture);
		WRITEI("showFilter", showFilter);
		WRITEI("showIndices", showIndex);
		WRITEI("showVertices", showVertex);
		WRITEI("showMode", showBack);
		WRITEI("showMesh", showMesh);
		WRITEI("showHidden", showHidden);
		WRITEI("showHiddenVertex", showHiddenVertex);
		WRITEI("showTimeOverlay", showTime);
		//WRITEI("texColormap", showColormap);
		f->Write("texColormap:1 1 1 1\n");
		WRITED("translation", transStep);
		WRITEI("dispNumLines", dispNumHits);
		WRITEI("dispNumLeaks", dispNumLeaks);
		WRITEI("dirShow", showDir);
		f->Write("dirNorme:"); f->Write(geom->GetNormeRatio(), "\n");
		f->Write("dirAutoNormalize:"); f->Write(geom->GetAutoNorme(), "\n");
		f->Write("dirCenter:"); f->Write(geom->GetCenterNorme(), "\n");

		WRITED("angle", angleStep);
		f->Write("autoScale:"); f->Write(geom->texAutoScale, "\n");
		f->Write("autoScale_include_constant_flow:"); f->Write(geom->texAutoScaleIncludeConstantFlow, "\n");

		f->Write("textures_min_pressure_all:");
		f->Write(geom->texture_limits[0].autoscale.min.all, "\n");
		f->Write("textures_min_pressure_moments_only:");
		f->Write(geom->texture_limits[0].autoscale.min.moments_only, "\n");
		f->Write("textures_max_pressure_all:");
		f->Write(geom->texture_limits[0].autoscale.max.all, "\n");
		f->Write("textures_max_pressure_moments_only:");
		f->Write(geom->texture_limits[0].autoscale.max.moments_only, "\n");

		f->Write("textures_min_impingement_all:");
		f->Write(geom->texture_limits[1].autoscale.min.all, "\n");

		f->Write("textures_min_impingement_moments_only:");
		f->Write(geom->texture_limits[1].autoscale.min.moments_only, "\n");
		f->Write("textures_max_impingement_all:");
		f->Write(geom->texture_limits[1].autoscale.max.all, "\n");
		f->Write("textures_max_impingement_moments_only:");
		f->Write(geom->texture_limits[1].autoscale.max.moments_only, "\n");

		f->Write("textures_min_density_all:");
		f->Write(geom->texture_limits[2].autoscale.min.all, "\n");
		f->Write("textures_min_density_moments_only:");
		f->Write(geom->texture_limits[2].autoscale.min.moments_only, "\n");
		f->Write("textures_max_density_all:");
		f->Write(geom->texture_limits[2].autoscale.max.all, "\n");
		f->Write("textures_max_density_moments_only:");
		f->Write(geom->texture_limits[2].autoscale.max.moments_only, "\n");

#ifdef _DEBUG
		f->Write("processNum:");f->Write(numCPU, "\n");
#else
		f->Write("processNum:"); f->Write(worker.GetProcNumber(), "\n");
#endif
		f->Write("recents:{\n");
		for (int i = 0; i < nbRecent; i++) {
			f->Write("\"");
			f->Write(recents[i]);
			f->Write("\"\n");
		}
		f->Write("}\n");

		f->Write("cdir:\""); f->Write(currentDir); f->Write("\"\n");
		f->Write("cseldir:\""); f->Write(currentSelDir); f->Write("\"\n");
		f->Write("autonorme:"); f->Write(geom->GetAutoNorme(), "\n");
		f->Write("centernorme:"); f->Write(geom->GetCenterNorme(), "\n");
		f->Write("normeratio:"); f->Write((double)(geom->GetNormeRatio()), "\n");

		f->Write("autoSaveFrequency:"); f->Write(autoSaveFrequency, "\n");
		f->Write("autoSaveSimuOnly:"); f->Write(autoSaveSimuOnly, "\n");
		f->Write("checkForUpdates:"); f->Write(/*checkForUpdates*/ 0, "\n"); //Deprecated
		f->Write("autoUpdateFormulas:"); f->Write(autoUpdateFormulas, "\n");
		f->Write("compressSavedFiles:"); f->Write(compressSavedFiles, "\n");
		f->Write("sh.gasMass:"); f->Write(worker.wp.gasMass, "\n");
		f->Write("expandShortcutPanel:"); f->Write(!shortcutPanel->IsClosed(), "\n");

		WRITEI("hideLot", hideLot);
		f->Write("lowFluxMode:"); f->Write(worker.ontheflyParams.lowFluxMode, "\n");
		f->Write("lowFluxCutoff:"); f->Write(worker.ontheflyParams.lowFluxCutoff, "\n");
		f->Write("leftHandedView:"); f->Write(leftHandedView, "\n");
	}
	catch (Error &err) {
		GLMessageBox::Display(err.GetMsg(), "Error saving config file", GLDLG_OK, GLDLG_ICONWARNING);
	}

	SAFE_DELETE(f);

}
/*
void MolFlow::calcFlow() {
	double sticking=1;
	double area;
	double outgassing;
	double temperature;
	//double mass;

	//facetSticking->GetNumber(&sticking);
	facetArea->GetNumber(&area);
	facetTemperature->GetNumber(&temperature);
	//facetMass->GetNumber(&mass);

	outgassing = 1 * sticking*area / 10.0 / 4.0*sqrt(8.0*8.31*temperature / PI / (worker.wp.gasMass*0.001));
	facetPumping->SetText(outgassing);
	return;
}*/

/*void MolFlow::calcSticking() {
	double sticking;
	double area;
	double outgassing;
	double temperature;
	//double mass;

	facetPumping->GetNumber(&outgassing);
	facetArea->GetNumber(&area); //area is in units of [cm^2]
	facetTemperature->GetNumber(&temperature);
	//facetMass->GetNumber(&mass);

	sticking = abs(outgassing / (area / 10.0)*4.0*sqrt(1.0 / 8.0 / 8.31 / (temperature)*PI*(worker.wp.gasMass*0.001)));
	//outgassing is here the variable, where the pumping speed [l/s] is stored.
	//liters are converted in m^3 and cm^3, therefore factor 10 remains!
	

	//if (sticking<=1.0) {
	facetSticking->SetText(sticking);
	//}
	//else { //Sticking: max. 1.0
	//	SetParam(facetSticking,1.0);
	//	calcFlow();
	//}
	return;
}*/
/*
void MolFlow::calcStickingnew() {
	llong covering;
	double sticking;
	facetcovering->GetNumberSizeT(&covering);
	if (covering) {
		if (covering < 100) {
			sticking = 0;
		}
		else
		{
			sticking = 1;
		}
	}


	
	double s1 = 0.1;
	double s2 = 0.2;
	double E_ad = pow(10, -21);
	//double E_de = 1.5*pow(10, -21);
	double kb = 1.38 * pow(10, -23);

	double sticking=0.0;
	double temperature;
	double coverage;

	facetcoverage->GetNumber(&coverage);
	facetTemperature->GetNumber(&temperature);
	if (coverage){
		if (coverage < 1) { 
			sticking = (s1*(1 - coverage) + s2 * coverage)*(1 - exp(-E_ad / (kb*temperature)));
		}
		else
		{
			sticking = s2 * (1 - exp(-E_ad / (kb*temperature)));
		}
	}
	facetSticking->SetText(sticking);
}*/

void MolFlow::calcCoverage() {
	llong covering;
	facetcovering->GetNumberSizeT(&covering);
	facetcoverage->SetText(double((covering) / (calcNmono() / calcdNsurf())));

}

void MolFlow::calcCovering() {
	double coverage;
	facetcoverage->GetNumber(&coverage);
	facetcovering->SetText(llong(abs(coverage * calcNmono() / calcdNsurf())));

}

double MolFlow::calcNmono() {//Calculates the Number of (carbon equivalent) particles of one monolayer.
	double area;
	facetArea->GetNumber(&area); //area is in units of [cm^2] => has to be converted to [m^2]
	double Nmono = (area * 1E-4 /(pow(carbondiameter, 2)));
	return Nmono;
}

double MolFlow::calcdNsurf() {// Calculates the (carbon equivalent realtive) mass factor 
	return (worker.wp.gasMass / 12.011);
}

//Funktionen neu machen:
//double calcDesorption("Facet"){}
//double calcDesorptionRate("Facet"){}

double MolFlow::calcDesorption() {
	double d = 1;
	double E_de = 1.5E-21;
	double coverage;
	double temperature;

	facetcoverage->GetNumber(&coverage);
	facetTemperature->GetNumber(&temperature);

	double desorption = 1.0 / tau * pow(coverage, d) *exp(-E_de / (kb*temperature));
	return desorption;
}

void MolFlow::calcDesorptionRate() {
	double temperature;
	facetTemperature->GetNumber(&temperature);

	double desorption = calcDesorption();
	double desorptionRate = desorption * (calcNmono() / calcdNsurf()) * 1.38E-23* temperature;
	desorptionRate = desorptionRate<0.0? 0.0 : desorptionRate;
	facetDesRate->SetText(desorptionRate);
	
}

void MolFlow::calcVolume() {

	Facet *f;
	double area = 0.0; double facetdepth;


	Geometry *geom = worker.GetGeometry();
	std::vector<size_t> selectedFacets = geom->GetSelectedFacets();
	for (size_t sel = 0; sel < selectedFacets.size();sel++) {
		f = geom->GetFacet(selectedFacets[sel]);
		area += f->GetArea() * (f->sh.is2sided ? 0.5 : 1);
	}

	facetDepth->GetNumber(&facetdepth);
	facetVolume->SetText(facetdepth*area);

}

void MolFlow::calcDepth() {
	Facet *f;
	double area = 0.0; double volume;

	Geometry *geom = worker.GetGeometry();
	std::vector<size_t> selectedFacets = geom->GetSelectedFacets();
	for (size_t sel = 0; sel < selectedFacets.size();sel++) {
		f = geom->GetFacet(selectedFacets[sel]);
		area += f->GetArea() * (f->sh.is2sided ? 0.5 : 1);
	}

	facetVolume->GetNumber(&volume);
	facetDepth->SetText(volume / area);
}

//TODO: correct transformation?
void MolFlow::calcConcentration() {
	double mass, volume;
	facetVolume->GetNumber(&volume);
	facetMass->GetNumber(&mass);
	if (volume == 0)
		facetConcentration->SetText(0);
	else
		facetConcentration->SetText(mass / volume);

}

void MolFlow::calcMass() {
	double concentration, volume;
	facetVolume->GetNumber(&volume);
	facetConcentration->GetNumber(&concentration);
	facetMass->SetText(concentration*volume);
}



void MolFlow::CrashHandler(Error *e) {
	char tmp[1024];
	sprintf(tmp, "Well, that's emberassing. Molflow crashed and will exit now.\nBefore that, an autosave will be attempted.\nHere is the error info:\n\n%s", (char *)e->GetMsg());
	GLMessageBox::Display(tmp, "Main crash handler", GLDLG_OK, GLDGL_ICONDEAD);
	try {
		if (AutoSave(true))
			GLMessageBox::Display("Good news, autosave worked!", "Main crash handler", GLDLG_OK, GLDGL_ICONDEAD);
		else
			GLMessageBox::Display("Sorry, I couldn't even autosave.", "Main crash handler", GLDLG_OK, GLDGL_ICONDEAD);
	}
	catch (Error &e) {
		e.GetMsg();
		GLMessageBox::Display("Sorry, I couldn't even autosave.", "Main crash handler", GLDLG_OK, GLDGL_ICONDEAD);
	}
}

bool MolFlow::EvaluateVariable(VLIST *v) {
	bool ok = true;
	Geometry* geom = worker.GetGeometry();
	size_t nbFacet = geom->GetNbFacet();
	int idx;

	if ((idx = GetVariable(v->name, "A")) > 0) {
		ok = (idx > 0 && idx <= nbFacet);
		if (ok) v->value = geom->GetFacet(idx - 1)->facetHitCache.hit.nbAbsEquiv;
	}
	else if ((idx = GetVariable(v->name, "D")) > 0) {
		ok = (idx > 0 && idx <= nbFacet);
		if (ok) v->value = (double)geom->GetFacet(idx - 1)->facetHitCache.hit.nbDesorbed;
	}
	else if ((idx = GetVariable(v->name, "MCH")) > 0) {
		ok = (idx > 0 && idx <= nbFacet);
		if (ok) v->value = (double)geom->GetFacet(idx - 1)->facetHitCache.hit.nbMCHit;
	}
	else if ((idx = GetVariable(v->name, "H")) > 0) {
		ok = (idx > 0 && idx <= nbFacet);
		if (ok) v->value = (double)geom->GetFacet(idx - 1)->facetHitCache.hit.nbHitEquiv;
	}
	else if ((idx = GetVariable(v->name, "P")) > 0) {
		ok = (idx > 0 && idx <= nbFacet);
		if (ok) v->value = geom->GetFacet(idx - 1)->facetHitCache.hit.sum_v_ort *
			worker.GetMoleculesPerTP(worker.displayedMoment)*1E4 / geom->GetFacet(idx - 1)->GetArea() * (worker.wp.gasMass / 1000 / 6E23)*0.0100;
	}
	else if ((idx = GetVariable(v->name, "DEN")) > 0) {
		ok = (idx > 0 && idx <= nbFacet);
		if (ok) {
			Facet *f = geom->GetFacet(idx - 1);
			v->value = f->DensityCorrection() * f->facetHitCache.hit.sum_1_per_ort_velocity /
				f->GetArea() *
				worker.GetMoleculesPerTP(worker.displayedMoment)*1E4;
		}
	}
	else if ((idx = GetVariable(v->name, "Z")) > 0) {
		ok = (idx > 0 && idx <= nbFacet);
		if (ok) v->value = geom->GetFacet(idx - 1)->facetHitCache.hit.nbHitEquiv /
			geom->GetFacet(idx - 1)->GetArea() *
			worker.GetMoleculesPerTP(worker.displayedMoment)*1E4;
	}
	else if ((idx = GetVariable(v->name, "V")) > 0) {
		ok = (idx > 0 && idx <= nbFacet);
		if (ok) /*v->value = 4.0*(double)(geom->GetFacet(idx - 1)->facetHitCache.hit.nbMCHit + geom->GetFacet(idx - 1)->facetHitCache.hit.nbDesorbed) /
			geom->GetFacet(idx - 1)->facetHitCache.hit.sum_1_per_ort_velocity;*/
			v->value = (geom->GetFacet(idx - 1)->facetHitCache.hit.nbHitEquiv + static_cast<double>(geom->GetFacet(idx - 1)->facetHitCache.hit.nbDesorbed)) / geom->GetFacet(idx - 1)->facetHitCache.hit.sum_1_per_velocity;
	}
	else if ((idx = GetVariable(v->name, "T")) > 0) {
		ok = (idx > 0 && idx <= nbFacet);
		if (ok) v->value = geom->GetFacet(idx - 1)->sh.temperature;
	}
	else if ((idx = GetVariable(v->name, "AR")) > 0) {
		ok = (idx > 0 && idx <= nbFacet);
		if (ok) v->value = geom->GetFacet(idx - 1)->sh.area;
	}
	else if (_stricmp(v->name, "SUMDES") == 0) {
		v->value = (double)worker.globalHitCache.globalHits.hit.nbDesorbed;
	}
	else if (_stricmp(v->name, "SUMABS") == 0) {
		v->value = worker.globalHitCache.globalHits.hit.nbAbsEquiv;
	}
	else if (_stricmp(v->name, "SUMMCHIT") == 0) {
		v->value = (double)worker.globalHitCache.globalHits.hit.nbMCHit;
	}
	else if (_stricmp(v->name, "SUMHIT") == 0) {
		v->value = worker.globalHitCache.globalHits.hit.nbHitEquiv;
	}
	else if (_stricmp(v->name, "MPP") == 0) {
		v->value = worker.globalHitCache.distTraveled_total / (double)worker.globalHitCache.globalHits.hit.nbDesorbed;
	}
	else if (_stricmp(v->name, "MFP") == 0) {
		v->value = worker.globalHitCache.distTraveledTotal_fullHitsOnly / worker.globalHitCache.globalHits.hit.nbHitEquiv;
	}
	else if (_stricmp(v->name, "DESAR") == 0) {
		double sumArea = 0.0;
		for (int i2 = 0; i2 < geom->GetNbFacet(); i2++) {
			Facet *f_tmp = geom->GetFacet(i2);
			if (f_tmp->sh.desorbType) sumArea += f_tmp->GetArea();
		}
		v->value = sumArea;
	}
	else if (_stricmp(v->name, "ABSAR") == 0) {
		double sumArea = 0.0;

		for (int i2 = 0; i2 < geom->GetNbFacet(); i2++) {
			Facet *f_tmp = geom->GetFacet(i2);
			if (f_tmp->sh.sticking > 0.0) sumArea += f_tmp->GetArea()*f_tmp->sh.opacity;
		}
		v->value = sumArea;
	}
	else if (_stricmp(v->name, "QCONST") == 0) {
		v->value = worker.wp.finalOutgassingRate_Pa_m3_sec*10.00; //10: Pa*m3/sec -> mbar*l/s
	}
	else if (_stricmp(v->name, "QCONST_N") == 0) {
		v->value = worker.wp.finalOutgassingRate;
	}
	else if (_stricmp(v->name, "NTOT") == 0) {
		v->value = worker.wp.totalDesorbedMolecules;
	}
	else if (_stricmp(v->name, "sh.gasMass") == 0) {
		v->value = worker.wp.gasMass;
	}
	else if (_stricmp(v->name, "KB") == 0) {
		v->value = 1.3806504e-23;
	}
	else if (_stricmp(v->name, "R") == 0) {
		v->value = 8.314472;
	}
	else if (_stricmp(v->name, "Na") == 0) {
		v->value = 6.02214179e23;
	}
	else if ((beginsWith(v->name, "SUM(") || beginsWith(v->name, "sum(") || beginsWith(v->name, "AVG(") || beginsWith(v->name, "avg(")) && endsWith(v->name, ")")) {
		bool avgMode = (beginsWith(v->name, "AVG(") || beginsWith(v->name, "avg(")); //else SUM mode
		std::string inside = v->name; inside.erase(0, 4); inside.erase(inside.size() - 1, 1);
		std::vector<std::string> tokens = SplitString(inside,',');
		if (!Contains({ 2,3 }, tokens.size()))
			return false;
		if (avgMode) {
			if (!Contains({ "P","DEN","Z","p","den","z" }, tokens[0]))
				return false;
		}
		else {
			if (!Contains({ "MCH","H","D","A","mch","h","d","a" }, tokens[0]))
				return false;
		}
		std::vector<size_t> facetsToSum;
		if (tokens.size() == 3) { // Like SUM(H,3,6) = H3 + H4 + H5 + H6
			size_t startId, endId, pos;
			try {
				startId = std::stol(tokens[1], &pos); if (pos != tokens[1].size() || startId > geom->GetNbFacet() || startId == 0) return false;
				endId = std::stol(tokens[2], &pos); if (pos != tokens[2].size() || endId > geom->GetNbFacet() || endId == 0) return false;
			}
			catch (...) {
				return false;
			}
			if (!(startId < endId)) return false;
			facetsToSum = std::vector<size_t>(endId-startId+1);
			std::iota(facetsToSum.begin(), facetsToSum.end(), startId-1);
		}
		else { //Selection group
			if (!(beginsWith(tokens[1], "S") || beginsWith(tokens[1], "s"))) return false;
			std::string selIdString = tokens[1]; selIdString.erase(0, 1);
			if (Contains({ "EL","el" }, selIdString)) { //Current selections
				facetsToSum = geom->GetSelectedFacets();
			}
			else {
				size_t selGroupId, pos;
				try {
					selGroupId = std::stol(selIdString, &pos); if (pos != selIdString.size() || selGroupId > selections.size() || selGroupId == 0) return false;
				}
				catch (...) {
					return false;
				}
				facetsToSum = selections[selGroupId - 1].selection;
			}
		}
		llong sumLL=0;
		double sumD=0.0;
		double sumArea = 0.0; //We average by area
		for (auto& sel : facetsToSum) {
			if (Contains({"MCH", "mch"},tokens[0])) {
				sumLL+=geom->GetFacet(sel)->facetHitCache.hit.nbMCHit;
			} 
			else if (Contains({ "H", "h" }, tokens[0])) {
				sumD += geom->GetFacet(sel)->facetHitCache.hit.nbHitEquiv;
			}
			else if (Contains({ "D", "d" }, tokens[0])) {
				sumLL+=geom->GetFacet(sel)->facetHitCache.hit.nbDesorbed;
			} else if (Contains({ "A", "a" }, tokens[0])) {
				sumD += geom->GetFacet(sel)->facetHitCache.hit.nbAbsEquiv;
			} else if (Contains({ "AR", "ar" }, tokens[0])) {
				sumArea += geom->GetFacet(sel)->GetArea();
			}
			else if (Contains({ "P", "p" }, tokens[0])) {
				sumD+= geom->GetFacet(sel)->facetHitCache.hit.sum_v_ort *
					(worker.wp.gasMass / 1000 / 6E23)*0.0100;
				sumArea += geom->GetFacet(sel)->GetArea();
			} else if (Contains({ "DEN", "den" }, tokens[0])) {
				Facet *f = geom->GetFacet(sel);
				sumD += f->DensityCorrection() * f->facetHitCache.hit.sum_1_per_ort_velocity;
				sumArea += geom->GetFacet(sel)->GetArea();
			} else if (Contains({ "Z", "z" }, tokens[0])) {
				sumD += geom->GetFacet(sel)->facetHitCache.hit.nbHitEquiv;
				sumArea += geom->GetFacet(sel)->GetArea();
			} else return false;
		}
		if (avgMode) v->value=sumD * worker.GetMoleculesPerTP(worker.displayedMoment)*1E4 / sumArea;
		else if (Contains({ "AR", "ar" }, tokens[0])) v->value = sumArea;
		else if (Contains({ "H", "h", "A", "a" }, tokens[0])) v->value = sumD;
		else v->value = static_cast<double>(sumLL); //Only one conversion at the end (instead of at each summing operation)
	}
	else ok = false;
	return ok;
}

void MolFlow::UpdatePlotters() {
	if (pressureEvolution) pressureEvolution->Update(m_fTime, true);
	if (timewisePlotter) timewisePlotter->Update(m_fTime, true);
	if (profilePlotter) profilePlotter->Update(m_fTime, true);
	if (texturePlotter) texturePlotter->Update(m_fTime, true);
	if (histogramPlotter) histogramPlotter->Update(m_fTime,true);
}

void MolFlow::UpdateFacetHits(bool allRows) {
	char tmp[256];
	Geometry *geom = worker.GetGeometry();

	try {
		// Facet list
		if (geom->IsLoaded()) {

			int sR, eR;
			if (allRows)
			{
				sR = 0;
				eR = (int)facetList->GetNbRow() - 1;
			}
			else
			{
				facetList->GetVisibleRows(&sR, &eR);
			}

			if (worker.displayedMoment == 0) {
				int colors[] = { COLOR_BLACK, COLOR_BLACK, COLOR_BLACK, COLOR_BLACK };
				facetList->SetColumnColors(colors);
			}
			else
			{
				int colors[] = { COLOR_BLACK, COLOR_BLUE, COLOR_BLUE, COLOR_BLUE };
				facetList->SetColumnColors(colors);
			}

			for (int i = sR; i <= eR; i++) {
				int facetId = facetList->GetValueInt(i, 0) - 1;
				if (facetId == -2) facetId = (int)i;
				if (i >= geom->GetNbFacet()) {
					char errMsg[512];
					sprintf(errMsg, "Molflow::UpdateFacetHits()\nError while updating facet hits. Was looking for facet #%d in list.\nMolflow will now autosave and crash.", i + 1);
					GLMessageBox::Display(errMsg, "Error", GLDLG_OK, GLDLG_ICONERROR);
					AutoSave();
				}
				Facet *f = geom->GetFacet(facetId);
				sprintf(tmp, "%d", facetId + 1);
				facetList->SetValueAt(0, i, tmp);
				switch (modeCombo->GetSelectedIndex()) {
				case MC_MODE:
					facetList->SetColumnLabel(1, "Hits");
					sprintf(tmp, "%I64d", f->facetHitCache.hit.nbMCHit);
					facetList->SetValueAt(1, i, tmp);
					sprintf(tmp, "%I64d", f->facetHitCache.hit.nbDesorbed);
					facetList->SetValueAt(2, i, tmp);
					sprintf(tmp, "%g", f->facetHitCache.hit.nbAbsEquiv);
					facetList->SetValueAt(3, i, tmp);
					break;
				case AC_MODE:
					facetList->SetColumnLabel(1, "Density");
					sprintf(tmp, "%g", f->facetHitCache.density.value);
					facetList->SetValueAt(1, i, tmp);

					sprintf(tmp, "%g", f->facetHitCache.density.desorbed);
					facetList->SetValueAt(2, i, tmp);
					sprintf(tmp, "%g", f->facetHitCache.density.absorbed);
					facetList->SetValueAt(3, i, tmp);
					break;
				}
			}

		}
	}
	catch (Error &e) {
		char errMsg[512];
		sprintf(errMsg, "%s\nError while updating facet hits", e.GetMsg());
		GLMessageBox::Display(errMsg, "Error", GLDLG_OK, GLDLG_ICONERROR);
	}

}
