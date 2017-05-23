/*
File:        GeometryViewer.cpp
Description: Geometry 3D Viewer component
Program:     MolFlow
Author:      R. KERSEVAN / J-L PONS / M ADY
Copyright:   E.S.R.F / CERN

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.
*/
#include "GLApp/GLApp.h"
#include "GLApp/GLWindowManager.h"
#include "GeometryViewer.h"
#include "GLApp/GLToolkit.h"
#include "GLApp/GLMatrix.h"
#include "GLApp/MathTools.h"
#include <math.h>

#include "GLApp/GLButton.h"
#include "GLApp/GLCombo.h"
#include "GLApp/GLLabel.h"
#include "GLApp\GLGradient.h"
//#include <malloc.h>

#ifdef MOLFLOW
#include "MolFlow.h"
extern MolFlow *mApp;
#endif

#ifdef SYNRAD
#include "SynRad.h"
extern SynRad *mApp;
#endif

void GeometryViewer::SetBounds(int x, int y, int width, int height) {
	if (this->posX != x || this->posY != y || this->width != width || this->height != height) {
		GLComponent::SetBounds(x, y, width, height);

		frontBtn->SetBounds(posX + 5, posY + height - 22, 33, 19);
		topBtn->SetBounds(posX + 40, posY + height - 22, 33, 19);
		sideBtn->SetBounds(posX + 75, posY + height - 22, 33, 19);
		projCombo->SetBounds(posX + 110, posY + height - 22, 50, 19);
		coordLab->SetBounds(posX + 162, posY + height - 20, 100, 18);
		facetSearchState->SetBounds(posX + 10, posY + 10, 90, 19);
		timeLabel->SetBounds(posX + width - 200, posY + height - 100, 100, 50);

		autoBtn->SetBounds(posX + width - 122, posY + height - 22, 19, 19);
		selBtn->SetBounds(posX + width - 102, posY + height - 22, 19, 19);

		selVxBtn->SetBounds(posX + width - 82, posY + height - 22, 19, 19);
		zoomBtn->SetBounds(posX + width - 62, posY + height - 22, 19, 19);
		handBtn->SetBounds(posX + width - 42, posY + height - 22, 19, 19);
		sysBtn->SetBounds(posX + width - 22, posY + height - 22, 19, 19);

		hideLotlabel->SetBounds(posX + 10, posY + height - 67, 0, 19);
		capsLockLabel->SetBounds(posX + 10, posY + height - 47, 0, 19);

		toolBack->SetBounds(posX + 1, posY + height - DOWN_MARGIN, width - 2, DOWN_MARGIN - 1);
		SetCurrentView(view);
	}
}

void GeometryViewer::DrawLinesAndHits() {

	// Lines
	if (showLine && mApp->worker.hitCacheSize) {

		glDisable(GL_TEXTURE_2D);
		glDisable(GL_LIGHTING);
		glDisable(GL_CULL_FACE);



		int count = 0;
		while (count < MIN(dispNumHits, mApp->worker.hitCacheSize) && mApp->worker.hitCache[count].type != 0) {

			if (mApp->whiteBg) { //whitebg
				glColor3f(0.2f, 0.7f, 0.2f);
			}
			else {
				glColor3f(0.5f, 1.0f, 0.5f);
			}

			if (mApp->antiAliasing) {
				glEnable(GL_BLEND);
				glEnable(GL_LINE_SMOOTH);
			}

			glBegin(GL_LINE_STRIP);
			while (count < MIN(dispNumHits, mApp->worker.hitCacheSize) && mApp->worker.hitCache[count].type != HIT_ABS) {

				//teleport routine
				if (mApp->worker.hitCache[count].type == HIT_TELEPORT) {
					glVertex3d(mApp->worker.hitCache[count].pos.x, mApp->worker.hitCache[count].pos.y, mApp->worker.hitCache[count].pos.z);
					glEnd();
					if (showTP) {
						if (!mApp->whiteBg) {
							glColor3f(1.0f, 0.7f, 0.2f);
						}
						else {
							glColor3f(1.0f, 0.0f, 1.0f);
						}
						glPushAttrib(GL_ENABLE_BIT);

						glLineStipple(1, 0x0101);
						glEnable(GL_LINE_STIPPLE);
						glBegin(GL_LINE_STRIP);
						glVertex3d(mApp->worker.hitCache[count].pos.x, mApp->worker.hitCache[count].pos.y, mApp->worker.hitCache[count].pos.z); //source point
						count++;
						glVertex3d(mApp->worker.hitCache[count].pos.x, mApp->worker.hitCache[count].pos.y, mApp->worker.hitCache[count].pos.z);  //teleport dest.
						glEnd();
						glPopAttrib();

						if (mApp->whiteBg) { //whitebg
							glColor3f(0.2f, 0.7f, 0.2f);
						}
						else {
							glColor3f(0.5f, 1.0f, 0.5f);
						}
						glBegin(GL_LINE_STRIP);
					}
					else {

						//glVertex3d(hitCache[count].pos.x , hitCache[count].pos.y , hitCache[count].pos.z); //source point
						count++;
						glBegin(GL_LINE_STRIP);
						glVertex3d(mApp->worker.hitCache[count].pos.x, mApp->worker.hitCache[count].pos.y, mApp->worker.hitCache[count].pos.z);  //teleport dest.

					}
				} //end treating teleport 

				if (mApp->worker.hitCache[count].type == HIT_LAST) { //pen up at cache refresh border
									glEnd();
									count++;
									glBegin(GL_LINE_STRIP);
				}
				else {
					glVertex3d(mApp->worker.hitCache[count].pos.x, mApp->worker.hitCache[count].pos.y, mApp->worker.hitCache[count].pos.z);
					count++;
				}
			}
			if (count < MIN(dispNumHits, mApp->worker.hitCacheSize) && mApp->worker.hitCache[count].type != 0) {
				glVertex3d(mApp->worker.hitCache[count].pos.x, mApp->worker.hitCache[count].pos.y, mApp->worker.hitCache[count].pos.z);
				count++;
			}
			glEnd();
			if (mApp->antiAliasing) {
				glDisable(GL_LINE_SMOOTH);
				glDisable(GL_BLEND);
			}
		}

	}

	// Hit
	if (showHit) {

		glDisable(GL_TEXTURE_2D);
		glDisable(GL_LIGHTING);
		glDisable(GL_BLEND);
		glDisable(GL_CULL_FACE);

		// Refl
		float pointSize = (bigDots) ? 2.0f : 1.0f;
		glPointSize(pointSize);
		if (mApp->whiteBg) { //whitebg
			glColor3f(0.2f, 0.2f, 0.2f);
		}
		else {
			glColor3f(0.0f, 1.0f, 0.0f);
		}
		glBegin(GL_POINTS);
		for (int i = 0; i < MIN(dispNumHits, mApp->worker.hitCacheSize); i++)
			if (mApp->worker.hitCache[i].type == HIT_REF)
				glVertex3d(mApp->worker.hitCache[i].pos.x, mApp->worker.hitCache[i].pos.y, mApp->worker.hitCache[i].pos.z);
		glEnd();

		// Moving Refl

		glPointSize(pointSize);
		/*if (mApp->whiteBg) { //whitebg
			glColor3f(0.2f, 0.2f, 0.2f);
		}
		else {*/
		glColor3f(1.0f, 0.0f, 1.0f);
		//}
		glBegin(GL_POINTS);
		for (int i = 0; i < MIN(dispNumHits, mApp->worker.hitCacheSize); i++)
			if (mApp->worker.hitCache[i].type == HIT_MOVING)
				glVertex3d(mApp->worker.hitCache[i].pos.x, mApp->worker.hitCache[i].pos.y, mApp->worker.hitCache[i].pos.z);
		glEnd();

		// Trans
		pointSize = (bigDots) ? 3.0f : 2.0f;
		glPointSize(pointSize);
		glColor3f(0.5f, 1.0f, 1.0f);
		glBegin(GL_POINTS);
		for (int i = 0; i < MIN(dispNumHits, mApp->worker.hitCacheSize); i++)
			if (mApp->worker.hitCache[i].type == HIT_TRANS)
				glVertex3d(mApp->worker.hitCache[i].pos.x, mApp->worker.hitCache[i].pos.y, mApp->worker.hitCache[i].pos.z);
		glEnd();

		// Teleport
		if (showTP) {
			//pointSize=(bigDots)?3.0f:2.0f;
			glPointSize(pointSize);
			if (!mApp->whiteBg) {
				glColor3f(1.0f, 0.7f, 0.2f);
			}
			else {
				glColor3f(1.0f, 0.0f, 1.0f);
			}
			glBegin(GL_POINTS);
			for (int i = 0; i < MIN(dispNumHits, mApp->worker.hitCacheSize); i++)
				if (mApp->worker.hitCache[i].type == HIT_TELEPORT)
					glVertex3d(mApp->worker.hitCache[i].pos.x, mApp->worker.hitCache[i].pos.y, mApp->worker.hitCache[i].pos.z);
			glEnd();
		}

		// Abs
		glPointSize(pointSize);
		glColor3f(1.0f, 0.0f, 0.0f);
		glBegin(GL_POINTS);
		for (int i = 0; i < MIN(dispNumHits, mApp->worker.hitCacheSize); i++)
			if (mApp->worker.hitCache[i].type == HIT_ABS)
				glVertex3d(mApp->worker.hitCache[i].pos.x, mApp->worker.hitCache[i].pos.y, mApp->worker.hitCache[i].pos.z);
		glEnd();

		// Des
		glColor3f(0.3f, 0.3f, 1.0f);
		glBegin(GL_POINTS);
		for (int i = 0; i < MIN(dispNumHits, mApp->worker.hitCacheSize); i++)
			if (mApp->worker.hitCache[i].type == HIT_DES)
				glVertex3d(mApp->worker.hitCache[i].pos.x, mApp->worker.hitCache[i].pos.y, mApp->worker.hitCache[i].pos.z);
		glEnd();

	}

}