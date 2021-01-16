/* 'sim2d' (c) F.Masset - 1995 - 2020. Codigo 2D de simulaciones
   N-cuerpos con metodo PIC y evaluacion del potencial por FFTs,
   enfocado a dinamica galactica. Codigo adaptado para fines de
   ensenyanza de la materia Dinamica y Estructura de Galaxias del
   Posgrado en Astrofisica de la UNAM. Para evitar problemas de
   formato, todos los comentarios en espanyol en las fuentes de este
   codigo son voluntariamente sin acentos */


#include "sim2d.h"

void InitKernelPotential(Grid *array)
{
  int             i, j, k, ii, jj, Nx, Ny;
  real            radius, radiusshort, Ax, Ay, *kernel;
  kernel = array->Field;
  Ax = array->Ax;
  Ay = array->Ay;
  Nx = array->Nx;
  Ny = array->Ny;
  Message("Potential kernel initialization...");
  kernel[0] = 0.;
  for (i = 1; i <= Nx; i++) {
    for (j = 1; j <= Ny; j++) {
      radiusshort = Ax * (float)(i - Nx / 2) - Ax/2.;
      radius = radiusshort * radiusshort;
      radiusshort = Ay * (float)(j - Ny / 2) - Ay/2.;
      radius += radiusshort * radiusshort;
      radius += ZSOFT * ZSOFT;
      /*
       * We use here the trick known as softened gravity in
       * order to smooth the potential on short ranges.
       * Zsoft has to be taken about the thickness of a
       * stellar disk
       */
      radius = sqrt(radius);
      ii = i - Nx / 2;
      jj = j - Ny / 2;
      /*
       * The translation above is made so as to put the
       * central point of the potential kernel in the (1,1)
       * corner the array, in order to prevent that the
       * convolution of density by this kernel leads also
       * to a translation
       */
      if (ii < 1)
	ii += Nx;
      if (jj < 1)
	jj += Ny;
      if (ii > Nx)
	ii -= Nx;
      if (jj > Ny)
	jj -= Ny;
      k = 2 * ((ii - 1) * Nx + jj - 1) + 1;
      kernel[k + 1] = 0.;
      kernel[k] = -G / radius;
    }
  }
  FourierDirecte (array);
  Message("done.\n");
}
