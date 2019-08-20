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
#ifndef _PARAMETERH_
#define _PARAMETERH_

#include "Distributions.h"
#include <string>
#include <vector>

class Parameter:public Distribution2D {
public:
	std::string name;
	bool fromCatalog;
	
	Parameter();

	template<class Archive>
	void serialize(Archive & archive)
	{
		archive(name,fromCatalog);
	}
};

/*
class StringClass:public Distribution2D {
public:
	std::string name;
	bool fromCatalog;
	StringClass();
};*/
#endif