/*
 * Fibracion de Hopf, proyeccion de circulos de Hopf del hiperespacio CxC a CxR
 * esto es para mostrarselo a Ximena :) 
 * compila 
 * gcc -lm -lGL -lglut -lGLU hopf.c -o hopf 
 *
 * beck@math.co.ro
 */

#include <GL/glut.h>
#include <complex.h>
#include <math.h>
#include <stdio.h>
float rot = 0.0f;

/* Vector en R4 */
typedef struct r4v
{
  double x1, x2, x3, x4;
} r4v_t;

/* Vector en R3 */
typedef struct r3v
{
  double x1, x2, x3;
} r3v_t;

/* Vector en CxR */
typedef struct crv
{
  complex z;
  double x;
} crv_t;

/* Vector en CxC */
typedef struct c2v
{
  complex z0, z1;
} c2v_t;

/* Identificamos RxRxRxR = R4 con el hiperespacio complejo CxC=C2 */

c2v_t
r4v_a_c2v (r4v_t u)
{
  c2v_t v;
  v.z0 = u.x1 + I * u.x2;
  v.z1 = u.x3 + I * u.x4;
  return v;
}

/* Identificamos RxRxR=R3 con el espacio de dimension 3 CxR */
crv_t
r3v_a_crv (r3v_t u)
{
  crv_t v;
  v.z = u.x1 + I * u.x2;
  v.x = u.x3;
  return v;
}

r3v_t
crv_a_r3v (crv_t u)
{
  r3v_t v;
  v.x1 = creal (u.z);
  v.x2 = cimag (u.z);
  v.x3 = u.x;
  return v;
}


/* Fibracion de Hopf rhoR:R4 |-> R3 */

r3v_t
rhoR (r4v_t u)
{
  r3v_t r;
  r.x1 = 2 * ((u.x1 * u.x2) + (u.x3 * u.x4));
  r.x2 = 2 * ((u.x1 * u.x4) - (u.x2 * u.x3));
  r.x3 = ((u.x1 * u.x1) + (u.x3 * u.x3)) - ((u.x2 * u.x2) + (u.x4 * u.x4));
  return r;
}

/* Fibracion de Hopf rho:CxC |-> CxR */

crv_t
rho (c2v_t u)
{
  crv_t v;
  v.z = 2 * u.z0 * conj (u.z1);
  v.x = (cabs (u.z0) * cabs (u.z0)) - (cabs (u.z1) * cabs (u.z1));
  return v;
}

r3v_t
proy_estereografica_s3_r3 (r4v_t u)
{
  r3v_t r;
  r.x1 = u.x1 / (u.x4 - 1);
  r.x2 = u.x2 / (u.x4 - 1);
  r.x3 = u.x3 / (u.x4 - 1);
  return r;
}


void
fibra (r3v_t p)
{
  double a, b, c, x, y, z, w, phi, theta, alpha, beta, r, k,s;
  r4v_t h;
  r3v_t rp;
  a = p.x1;
  b = p.x2;
  c = p.x3;
  s=2;
  alpha = sqrt ((1 + c) / 2);
  beta = sqrt ((1 - c) / 2);
  for (phi = 0; phi <= 2 * M_PI; phi += 0.005)
    {
      theta = atan (-a / b) - phi;
      w = alpha * cos (theta);
      r = acos (w) / M_PI;
      k = r / sqrt (1 - (w * w));
      h.x1 = w;
      h.x2 = k * alpha * sin (theta);
      h.x3 = k * beta * cos (phi);
      h.x4 = k * beta * sin (phi);
     glVertex3f (s*h.x2, s*h.x3, s*h.x4);
      h.x2 = -1*k * alpha * sin (theta);
      h.x3 = -1*k * beta * cos (phi);
      h.x4 = -1*k * beta * sin (phi);
     glVertex3f (s*h.x2, s*h.x3, s*h.x4);
    }
  return;
}


/* parametrizacion de la hiperesfera S3 psi y theta en [0,pi], phi en [0,2pi] */
r4v_t
punto_s3 (double psi, double theta, double phi)
{
  r4v_t z;
  z.x1 = cos (psi);
  z.x2 = sin (psi) * cos (theta);
  z.x3 = sin (psi) * sin (theta) * cos (phi);
  z.x4 = sin (psi) * sin (theta) * sin (phi);
  return z;
}

/* Imprime vector de R4 */
void
print_r4v (r4v_t v)
{
  printf ("(%.2f,%.2f,%.2f,%.2f)\n", v.x1, v.x2, v.x3, v.x4);
  return;
}

/* Imprime vector de R3 */
void
print_r3v (r3v_t v)
{
  printf ("(%.2f,%.2f,%.2f)\n", v.x1, v.x2, v.x3);
  return;
}

/* Imprime vector de C2 */
void
print_c2v (c2v_t v)
{
  printf ("(%.2f+%.2fi,%.2f+%.2fi)\n", creal (v.z0), cimag (v.z0),
	  creal (v.z1), cimag (v.z1));
  return;
}

/* Imprime vector de CxR */
void
print_crv (crv_t v)
{
  printf ("(%.2f+%.2fi,%.2f)\n", creal (v.z), cimag (v.z), v.x);
  return;
}

r3v_t
prod_escalar (double k, r3v_t u)
{
  r3v_t r;
  r.x1 = u.x1 * k;
  r.x2 = u.x2 * k;
  r.x3 = u.x3 * k;
  return r;
}

void
display (void)
{
  double psi, theta, phi,k1,a,c;
  r4v_t p, q;
  c2v_t z;
  crv_t w;
  r3v_t h, f;
  double k = 2.5, d, v;
  psi = theta = phi = 0;
  glLoadIdentity ();
  glTranslatef (0.0f, -0.0f, -8.0f);
  glRotatef (rot, 0.1f, 0.1f, 0.1f);
  glClear (GL_COLOR_BUFFER_BIT);

  glBegin (GL_POINTS);
  
phi=M_PI/2;
a=2;

/* descomenta si quieres dibujar la imagen inversa de varias circunferencias */
//for (phi = 0.0; phi < M_PI; phi += 0.5) {
        srand(phi*(phi+1));
        k1 = rand();
    for (theta = 0; theta < 2*M_PI; theta += 0.1)
	  {
            h.x1 = cos(theta)*sin(phi);
            h.x2 = sin(theta)*sin(phi);
            h.x3 = cos(phi);

            /* Loxodromica 
            c=atan(a*theta);
            h.x1 = cos(theta)/sqrt(1+(a*a*theta*theta));
            h.x2 = sin(theta)/sqrt(1+(a*a*theta*theta));
            h.x3 = (-1*a*theta)/sqrt(1+(a*a*theta*theta));
            */
        /*solo para el color */
            k1=theta;
	    glColor3f (cos (k1*k1*k1), sin (k1 / 2), 1);
             fibra(h);
	    glVertex3f (-3+h.x1, h.x2, h.x3);
	  }
/* Descomenta tambien si descomentaste arriba */
//}
	glutSwapBuffers ();



  glEnd ();
  rot += 1;			// velocidad
  glutSwapBuffers ();
}

void
reshape (int w, int h)
{
  glViewport (0, 0, (GLsizei) w, (GLsizei) h);
  glMatrixMode (GL_PROJECTION);
  glLoadIdentity ();
  gluPerspective (40, (GLfloat) w / (GLfloat) (h), 1.5, 15.0);
  glMatrixMode (GL_MODELVIEW);
  glLoadIdentity ();
  gluLookAt (0, 0, 0, 0, 0, 0, 0.0, 0, 0.0);
}

int
main (int argc, char **argv)
{
  glutInit (&argc, argv);
 // glPointSize(1);
  glutInitDisplayMode (GLUT_SINGLE | GLUT_RGB);
  glutInitWindowSize (1300, 800);
  glutInitWindowPosition (0, 0);
  glutCreateWindow ("test");
  glClearColor (0, 0, 0, 0);
  glShadeModel (GL_FLAT);
  glMatrixMode (GL_PROJECTION);
  glutReshapeFunc (reshape);
  glutDisplayFunc (display);
  glutIdleFunc (display);
  glutMainLoop ();
  return 0;
}
