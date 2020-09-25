#ifndef INTERACTIONS_H
#define INTERACTIONS_H
#define W 920 
#define H 920
#define DELTA_P 0.1f
#define TITLE_STRING "Stability"
int sys = 0;
float param = 0.1f;
float z0 = 0.0;
int nlm = 320;
int nlmIndex = 0;
float zdelta = 0.05;
bool auto_move = false;
void keyboard(unsigned char key, int x, int y) {

  int nlmValues[7] = { 310, 320, 321, 322, 211, 210, 100 };

  if (key == 27) exit(0);
  if (key == '0') sys = 0;
  if (key == '1') sys = 1;
  if (key == '2') sys = 2;
  if (key == 'a') z0 += 0.05;
  if (key == 'd') z0 -= 0.05;
  if (key == 'n') {
    nlm = nlmValues[nlmIndex];
    nlmIndex = (nlmIndex + 1) % 7;
    z0 = 0.0;
  }
  if (key == 'w') auto_move = auto_move ? false : true; 
  glutPostRedisplay();
}

void handleSpecialKeypress(int key, int x, int y) {
  if (key == GLUT_KEY_DOWN) param -= DELTA_P;
  if (key == GLUT_KEY_UP) param += DELTA_P;
  glutPostRedisplay();
}

void idle(void) {
  if (auto_move == true) {
    if (z0 > 5.0 || z0 < -5.0) zdelta *= -1;
    z0 += zdelta;
    glutPostRedisplay();
  }
}


// no mouse interactions implemented for this app
void mouseMove(int x, int y) { return; }
void mouseDrag(int x, int y) { return; }

void printInstructions() {
  printf("Atom 3d\n");
  printf(" a and d keys to move the slize accros z axis:\n");
  printf(" w automate z-axis move\n");
  printf(" n to switch between different nlm\n");
  printf("Equation source from https://chem.libretexts.org");
}

#endif
