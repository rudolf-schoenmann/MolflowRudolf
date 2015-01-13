/*
  File:        Movement.h
  Description: Define moving parts in system
*/
#ifndef _MOVEMENTH_
#define _MOVEMENTH_

#include "GLApp/GLWindow.h"
#include "GLApp/GLButton.h"
#include "GLApp/GLTextField.h"
#include "GLApp/GLLabel.h"
#include "GLApp/GLToggle.h"
#include "GLApp/GLTitledPanel.h"

#include "Geometry.h"
#include "Worker.h"

#define MODE_NOMOVE 0
#define MODE_FIXED 1
#define MODE_ROTATING 2

class Movement : public GLWindow {

public:
  // Construction
	Movement(Geometry *geom, Worker *work);
  void ProcessMessage(GLComponent *src,int message);
  void Update();

  // Implementation
private:

  void UpdateToggle(GLComponent *src);
  
  GLLabel	*label1;
  GLTitledPanel	*groupBox1;
  GLLabel	*label16;
  GLLabel	*label14;
  GLLabel	*label15;
  GLTextField	*hzText;
  GLTextField	*degText;
  GLTextField	*rpmText;
  GLLabel	*label17;
  GLButton	*button2;
  GLLabel	*label10;
  GLLabel	*label11;
  GLTextField	*ryText;
  GLTextField	*rzText;
  GLLabel	*label12;
  GLTextField	*rxText;
  GLLabel	*label13;
  GLButton	*button1;
  GLLabel	*label6;
  GLLabel	*label7;
  GLTextField	*ayText;
  GLTextField	*azText;
  GLLabel	*label8;
  GLTextField	*axText;
  GLLabel	*label9;
  GLToggle	*checkBox3;
  GLLabel	*label5;
  GLLabel	*label4;
  GLTextField	*vyText;
  GLTextField	*vzText;
  GLLabel	*label3;
  GLTextField	*vxText;
  GLLabel	*label2;
  GLToggle	*checkBox2;
  GLToggle	*checkBox1;
  GLButton	*button3;
  GLButton	*button4;

  std::vector<GLTextField*> group1;
  std::vector<GLTextField*> group2;

  int mode;

  Geometry     *geom;
  Worker	   *work;

};

#endif /* _MovementH_ */
