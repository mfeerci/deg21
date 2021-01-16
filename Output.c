/* 'sim2d' (c) F.Masset - 1995 - 2020. Codigo 2D de simulaciones
   N-cuerpos con metodo PIC y evaluacion del potencial por FFTs,
   enfocado a dinamica galactica. Codigo adaptado para fines de
   ensenyanza de la materia Dinamica y Estructura de Galaxias del
   Posgrado en Astrofisica de la UNAM. Para evitar problemas de
   formato, todos los comentarios en espanyol en las fuentes de este
   codigo son voluntariamente sin acentos */

#include "sim2d.h"

extern boolean  Symmetrized;

void WriteDisk(Grid *array, int number)
{
  int             i, j, k, Nx, Ny,l=0;
  FILE           *dump;
  char 		name[80];
  real 		*ptr;
  real		temp[16384];
  ptr = array->Field;
  sprintf (name, "%s%s%d.dat", OUTPUTDIR, array->Name, number);
  Nx = array->Nx;
  Ny = array->Ny;
  dump = fopen(name, "w");
  if (dump == NULL) {
    fprintf(stderr, "Unable to open '%s'.\n", name);
    exit(1);
  }
  fprintf(stdout, "Writing '%s/%s%d.dat'...", OUTPUTDIR, array->Name, number);
  fflush (stdout);
  for (j = Ny/4+1; j <= 3*Ny/4; j++) {
    for (i = Nx/4+1; i <= 3*Nx/4; i++) {
      k = 2 * (Nx * (i - 1) + j) - 1;
      temp[l++] = (real)ptr[k];
    }
  }
  fwrite (temp, sizeof(real), Nx*Ny/4,dump);
  fclose(dump);
  Message ("done\n");
}


/*
 * The function below dumps to the disk the position and velocities of a
 * particles set as stored in memory. In the case of a symmetrized
 * distribution, we see one particle over two. This saves disk space and
 * improves lisibility of plots
 */

void DumpParticles(ParticleSet *set, int numero)
{
  int             i;
  FILE           *dump;
  char            name[80];
  real           *x, *y, *vx, *vy;
  x = set->Xparticles;
  y = set->Yparticles;
  vx = set->VXparticles;
  vy = set->VYparticles;
  sprintf(name, "%s%s%d.dat", OUTPUTDIR, set->Name, numero);
  dump = fopen(name, "w");
  if (dump == NULL) {
    fprintf(stderr, "Unable to open %s.\n", name);
    exit(1);
  }
  fprintf(stdout, "Writing '%s'...", name);
  fflush (stdout);
  for (i = 0; i < set->NumberParticles; i++)
    fprintf(dump, "%f\t%f\t%f\t%f\n", x[i], y[i], vx[i], vy[i]);
  fclose(dump);
  fprintf(stdout, "done\n");
}

/* The function below gives correct results only if XRESOL == YRESOL */

void DumpSurfaceDensity(ParticleSet *set, int numero)
{
  int             i, bin;
  FILE           *dump;
  char            name[80];
  real            r, *x, *y, surfdens[1000];
  x = set->Xparticles;
  y = set->Yparticles;
  sprintf(name, "%ssd%s%d.dat", OUTPUTDIR, set->Name, numero);
  dump = fopen(name, "w");
  if (dump == NULL) {
    fprintf(stderr, "Unable to open %s.\n", name);
    exit(1);
  }
  fprintf(stdout, "Writing '%s'...", name);
  fflush (stdout);
  for (i = 0; i < XSIZE; i++)
    surfdens[i] = 0.;
  for (i = 0; i < set->NumberParticles; i++) {
    r = sqrt(x[i] * x[i] + y[i] * y[i]);
    bin = (int) (r / XRESOL);
    if (bin < XSIZE)
      surfdens[bin] += 1;

  }
  for (i = 0; i < XSIZE; i++) {
    surfdens[i] *= set->Mass / M_PI / XRESOL / XRESOL / (2. * (real) i + 1.);
    if (Symmetrized == YES)
      surfdens[i] *= 2.;
    fprintf(dump, "%f\t%f\n", XRESOL * ((real) i + .5), surfdens[i]);
  }
  fclose(dump);
  fprintf(stdout, "done.\n");
}

void WriteDim () {
  char filename[200];
  FILE 	*dim;
  sprintf (filename, "%sdims.dat", OUTPUTDIR);
  if ((dim = fopen (filename, "w")) == NULL) {
    fprintf (stderr, "Unable to open %s. Program stopped\n", filename);
    exit (1);
  }
  fprintf (dim, "%d\t%d\n", XSIZE/2, YSIZE/2);
  fclose (dim);
}
