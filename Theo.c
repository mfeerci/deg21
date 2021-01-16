/* 'sim2d' (c) F.Masset - 1995 - 2020. Codigo 2D de simulaciones
   N-cuerpos con metodo PIC y evaluacion del potencial por FFTs,
   enfocado a dinamica galactica. Codigo adaptado para fines de
   ensenyanza de la materia Dinamica y Estructura de Galaxias del
   Posgrado en Astrofisica de la UNAM. Para evitar problemas de
   formato, todos los comentarios en espanyol en las fuentes de este
   codigo son voluntariamente sin acentos */

#include "sim2d.h"

static real     *BinneyTremaineS;
static real     *Circular;
static real     *Epicyclic;
static real     kmax;

static real     StaticEpicyclic[5000];
static real     StaticCircular[5000];

/* See Binney and Tremaine, Galactic Dynamics, p. 76 */

real Sigma(real r)
{
  return SIGMA0 * exp(-r / Rd);
}

void *AllocWorkArrays () {
  BinneyTremaineS = (real *)malloc ((NINTERV * 2 + 1) * sizeof(real));
  Circular        = (real *)malloc ((NINTERV * 2 + 1) * sizeof(real));
  Epicyclic       = (real *)malloc ((NINTERV * 2 + 1) * sizeof(real));
  if ((BinneyTremaineS == NULL) || (Circular == NULL) || (Epicyclic == NULL))
    return NULL;
  else
    return (void *)BinneyTremaineS;
}

void FreeWorkArrays () {
  free(BinneyTremaineS);
  free(Circular);
  free(Epicyclic);
}

real TheoMass ()
{
  real r, rmax, dr, sum = 0.;
  rmax = XRESOL * XSIZE * .25;
  dr = rmax/1000.;
  for (r = dr; r < rmax; r += dr)
    sum += Sigma(r) * r;
  sum += .5 * Sigma(r) * r;
  sum *= 2. * M_PI * dr;
  return sum;
}
	
real ComputeS(real k)
{
  int             i;
  real            sum = 0., h, hh, Rmin, Rmax, r;
  Rmax = XRESOL * XSIZE / 4.;
  Rmin = Rmax / 10. / NINTERV;
  h = Rmax / NINTERV;
  r = h;
  hh = h / 2;
  for (i = 0; i < NINTERV; i++) {
    sum += Bessj0(k * r) * r * Sigma(r);
    sum += 2. * Bessj0(k * (r - hh)) * (r - hh) * Sigma(r - hh);
    r += h;
  }
  sum += (Sigma(Rmin) * Rmin * Bessj0(k * Rmin) - Sigma(Rmax) * Rmax * Bessj0(k * Rmax)) * .5;
  sum *= h / 3.;
  return -2. * M_PI * G * sum;
  /*
   * This function corresponds to the evaluation of S(k) with a Simpson
   * method
   */
}

void FillArrayS()
{
  int             i;
  Message ("Hankel Transforming initial potential...");
  kmax = 2000. / XRESOL / XSIZE;
  for (i = 0; i <= 2 * NINTERV; i++)
    BinneyTremaineS[i] = ComputeS(kmax * (real) i / 2. / (real) NINTERV);
  Message ("done.\n");
}

void FillCircularVelocity()
{
  int             i, j;
  real            sum1, sum2, k, dk, radius, gammacentre, halo, k2halo,
    r;
  real            a, a2, b, bulge, k2bulge, vcbulge2;
  b = BULGERADIUS;
  dk = kmax / (real) NINTERV;
  Message ("Tabulating circular velocity and epicyclic frequency...");
  for (j = 0; j <= 2. * NINTERV; j++) {
    radius = XRESOL * XSIZE * (real) j / 8. / (real) NINTERV;
    r = sqrt(radius * radius + ZSOFT * ZSOFT);
    a2 = b * b + radius * radius;
    a = sqrt(a2);
    sum1 = sum2 = 0.;
    for (i = 1; i <= 2 * NINTERV; i += 2) {
      k = kmax * (real) i / 2. / (real) NINTERV;
      sum1 += BinneyTremaineS[i] * 2. * k * Bessj1(k * radius);
      sum2 += BinneyTremaineS[i] * 2. * k * k * Bessj0(k * radius);
      k += dk / 2;
      sum1 += BinneyTremaineS[i + 1] * k * Bessj1(k * radius);
      sum2 += BinneyTremaineS[i + 1] * k * k * Bessj0(k * radius);
    }
    sum1 -= BinneyTremaineS[2 * NINTERV] * .5 * kmax * Bessj1(kmax * radius);
    sum2 -= BinneyTremaineS[2 * NINTERV] * .5 * kmax * kmax * Bessj0(kmax * radius);
    sum1 *= dk / 3.;
    sum2 *= dk / 3.;
    sum2 = -2. * sum1 / r - sum2;
    sum1 *= -radius;
    gammacentre = G * CENTRALMASS / r / r;
    halo = HALOSPEEDLIM * HALOSPEEDLIM * HALOCORE / r / r * LHalo(r / HALOCORE);
    k2halo = halo / r;
    k2halo += HALOSPEEDLIM * HALOSPEEDLIM / (HALOCORE * HALOCORE + radius * radius);
    bulge = G * BULGEMASS / (b + a) / (b + a) / a;
    vcbulge2 = bulge * radius * radius;
    k2bulge = bulge * (4. - radius * radius / a2 * (b + 3. * a) / (b + a));
    Circular[j] = sqrt(sum1 + gammacentre * radius + halo * radius + vcbulge2);
    Epicyclic[j] = sqrt(sum2 + gammacentre / r + k2halo + k2bulge);
    StaticCircular[j] = sqrt(gammacentre * radius + halo * radius + vcbulge2);
    StaticEpicyclic[j] = sqrt(gammacentre / r + k2halo + k2bulge);
    /*
     * In the formulae above we should have only 'radius' in the
     * denominator but we smooth it in 'r' so as to avoid
     * problems near the center
     */
  }
  Message("done.\n");
}

real CircularVelocity(real radius)
{
  int             i;
  real            ii;
  ii = 8. * NINTERV * radius / XRESOL / XSIZE;
  i = (int) ii;
  if (i >= 2 * NINTERV) {
    i = 2 * NINTERV - 1;
    ii = 2. * (real) NINTERV - 1.;
  }
  return Circular[i] * (1. - ii + (real) i) + Circular[i + 1] * (ii - (real) i);
}

real EpicyclicFrequency(real radius)
{
  int             i;
  real            ii;
  ii = 8. * NINTERV * radius / XRESOL / XSIZE;
  i = (int) ii;
  if (i >= 2 * NINTERV) {
    i = 2 * NINTERV - 1;
    ii = 2. * (real) NINTERV - 1.;
  }
  return Epicyclic[i] * (1. - ii + (real) i) + Epicyclic[i + 1] * (ii - (real) i);
}

real StaticCircularVelocity(real radius)
{
  int             i;
  real            ii;
  ii = 8. * NINTERV * radius / XRESOL / XSIZE;
  i = (int) ii;
  if (i >= 2 * NINTERV) {
    i = 2 * NINTERV - 1;
    ii = 2. * (real) NINTERV - 1.;
  }
  return StaticCircular[i] * (1. - ii + (real) i) + StaticCircular[i + 1] * (ii - (real) i);
}


real StaticEpicyclicFrequency(real radius)
{
  int             i;
  real            ii;
  ii = 8. * NINTERV * radius / XRESOL / XSIZE;
  i = (int) ii;
  if (i >= 2 * NINTERV) {
    i = 2 * NINTERV - 1;
    ii = 2. * (real) NINTERV - 1.;
  }
  return StaticEpicyclic[i] * (1. - ii + (real) i) + StaticEpicyclic[i + 1] * (ii - (real) i);
}
