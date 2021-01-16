/* 'sim2d' (c) F.Masset - 1995 - 2020. Codigo 2D de simulaciones
   N-cuerpos con metodo PIC y evaluacion del potencial por FFTs,
   enfocado a dinamica galactica. Codigo adaptado para fines de
   ensenyanza de la materia Dinamica y Estructura de Galaxias del
   Posgrado en Astrofisica de la UNAM. Para evitar problemas de
   formato, todos los comentarios en espanyol en las fuentes de este
   codigo son voluntariamente sin acentos */

#include "sim2d.h"

void AdvanceParticles(ParticleSet *set, real dt)
{
  real *x, *y, *vx, *vy;
  int  Nb, i;
  x = set->Xparticles;
  y = set->Yparticles;
  vx = set->VXparticles;
  vy = set->VYparticles;
  Nb = set->NumberParticles;
  for (i = 0; i < Nb; i++) {
    /* POSGRADO: aqui escribir las dos lineas faltantes para
       actualizar las coordenadas de la particula i en funcion del
       paso de tiempo 'dt' y de las componentes 'vx[i]' y 'vy[i]' de
       su velocidad.  */
  }
}
