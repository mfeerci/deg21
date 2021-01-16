/* 'sim2d' (c) F.Masset - 1995 - 2020. Codigo 2D de simulaciones
   N-cuerpos con metodo PIC y evaluacion del potencial por FFTs,
   enfocado a dinamica galactica. Codigo adaptado para fines de
   ensenyanza de la materia Dinamica y Estructura de Galaxias del
   Posgrado en Astrofisica de la UNAM. Para evitar problemas de
   formato, todos los comentarios en espanyol en las fuentes de este
   codigo son voluntariamente sin acentos */

#include "sim2d.h"

static int      nn[5];
extern boolean  Symmetrized;

void FourierDirecte(Grid *array)
{
  real           *tab;
  nn[1] = array->Nx;
  nn[2] = array->Ny;
  tab = array->Field;
  Fourn(tab, nn, 2, 1);
}

void FourierInverse(Grid *array)
{
  real           *tab;
  nn[1] = array->Nx;
  nn[2] = array->Ny;
  tab = array->Field;
  Fourn(tab, nn, 2, -1);
}

void ErrorMessage (char *string)
{
  fprintf(stderr, "%s\n", string);
  exit(1);
}

void Message (char *msg) 
{
	fprintf (stdout, "%s", msg);
	fflush (stdout);
}

real LHalo (real x) 
{
	return (x-atan(x));
}

void Multiply(Grid *array1, Grid *array2, real normalization)
{
  int             i, j, nx, ny, l;
  real            a1, a2, b1, b2, *tab1, *tab2, *r1, *r2, *i1, *i2;
  nx = array1->Nx;
  ny = array2->Ny;
  tab1 = array1->Field;
  tab2 = array2->Field;
  for (i = 1; i <= nx; i++) {
    for (j = 1; j <= ny; j++) {
      l = 2 * (nx * (i - 1) + j) - 1;
      r1 = tab1 + l;
      i1 = r1 + 1;
      r2 = tab2 + l;
      i2 = r2 + 1;
      a1 = *r1;
      b1 = *i1;
      a2 = *r2;
      b2 = *i2;
      *r1 = a1 * a2 - b1 * b2;
      *r1 *= normalization;
      *i1 = a1 * b2 + a2 * b1;
      *i1 *= normalization;
    }
  }
}

Grid *CreateGrid(int Nx, real Ax, int Ny, real Ay, char *name)
{
  Grid           *array;
  real           *field;
  char           *string;
  int             i, j, k;

  array = (Grid *) malloc(sizeof(Grid));
  if (array == NULL)
    ErrorMessage ("Insufficient memory for Grid creation");
  field = (real *) malloc(sizeof(real) * (Nx + 1) * (Ny + 1) * 2);
  if (field == NULL)
    ErrorMessage ("Insufficient memory for Grid creation");
  string = (char *) malloc(sizeof(char) * 80);
  if (string == NULL)
    ErrorMessage ("Insufficient memory for Grid creation");
  strcpy(string, name);
  array->Field = field;
  array->Name = string;
  array->Nx = Nx;
  array->Ny = Ny;
  array->Ax = Ax;
  array->Ay = Ay;
  for (i = 1; i <= Nx; i++) {
    for (j = 1; j <= Ny; j++) {
      k = 2 * (Nx * (i - 1) + j) - 1;
      field[k] = 0.;
    }
  }
  return array;
}

ParticleSet *CreateParticleSet(int Npart, real mass, char *name)
{
  ParticleSet    *set;
  real           *x, *y, *vx, *vy;
  char           *string;
  int             i;
  x = (real *) malloc(sizeof(real) * Npart);
  if (x == NULL)
    ErrorMessage ("Insufficient memory for Particle Storage");
  y = (real *) malloc(sizeof(real) * Npart);
  if (y == NULL)
    ErrorMessage ("Insufficient memory for Particle Storage");
  vx = (real *) malloc(sizeof(real) * Npart);
  if (vx == NULL)
    ErrorMessage ("Insufficient memory for Particle Storage");
  vy = (real *) malloc(sizeof(real) * Npart);
  if (vy == NULL)
    ErrorMessage ("Insufficient memory for Particle Storage");
  string = (char *) malloc(sizeof(char) * 80);
  if (string == NULL)
    ErrorMessage ("Insufficient memory for Particle Storage");
  strcpy(string, name);
  set = (ParticleSet *) malloc(sizeof(ParticleSet));
  if (set == NULL)
    ErrorMessage ("Insufficient memory for Particle Storage");
  for (i = 0; i < Npart; i++) {
    x[i] = y[i] = vx[i] = vy[i] = 0.;
  }
  set->NumberParticles = Npart;
  set->Mass = mass;
  set->Xparticles = x;
  set->Yparticles = y;
  set->VXparticles = vx;
  set->VYparticles = vy;
  set->Name = string;
  return set;
}

void ForgetGrid (Grid *array)
{
  if (array != NULL) {
    if (array->Name  != NULL) free(array->Name);
    if (array->Field != NULL) free(array->Field);
    free (array);
  }
}

void CopyGrid (Grid *arraysrc, Grid *arraydest)
{
  int i, Nx, Ny;
  real *fieldsrc, *fielddest;
  /* We suppose that arrays have the same size */
  Nx = arraysrc->Nx;
  Ny = arraysrc->Ny;
  fieldsrc = arraysrc->Field;
  fielddest= arraydest->Field;
  for (i = 0; i < 2*(Nx+1)*(Ny+1); i++) {
    fielddest[i] = fieldsrc[i];
  }
}
