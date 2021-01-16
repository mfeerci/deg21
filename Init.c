/* 'sim2d' (c) F.Masset - 1995 - 2020. Codigo 2D de simulaciones
   N-cuerpos con metodo PIC y evaluacion del potencial por FFTs,
   enfocado a dinamica galactica. Codigo adaptado para fines de
   ensenyanza de la materia Dinamica y Estructura de Galaxias del
   Posgrado en Astrofisica de la UNAM. Para evitar problemas de
   formato, todos los comentarios en espanyol en las fuentes de este
   codigo son voluntariamente sin acentos */

#include "sim2d.h"

extern boolean  Symmetrized;

real DispersionProfile (real radius)
{
  return 3.36*TOOMRECENTRAL*G*Sigma(radius)/EpicyclicFrequency(radius);
}

real Q_Toomre (real radius)
{
  (void) (radius);		/* Avoid unused parameter warning */
  return TOOMRECENTRAL;
}

real TotalMass(ParticleSet    *particles)
{
  real            mass;
  int             number;

  number = particles->NumberParticles;
  mass = particles->Mass;
  if (Symmetrized == YES)
    mass *= 2.;
  return mass * number;
}

void InitStars (ParticleSet *stars)
{
  real            radius, phi_orbit, velocity, dispersion, MeanDispersion, varalea;
  real            phi_offset, Omega, kappa, phi_epicyclic, tempwork, DerivDisp;
  int             i, TotalStars=0, NbStars, Nx;
  real		*x_stars, *y_stars, *px_stars, *py_stars, Ax, velocityt;

  x_stars = stars->Xparticles;
  y_stars = stars->Yparticles;
  px_stars= stars->VXparticles;
  py_stars= stars->VYparticles;

  Nx = XSIZE;
  Ax = XRESOL;

  NbStars = stars->NumberParticles;

  Message("Initializing stars positions and velocities...");

  for (i = 0; i < NbStars; i++) {

    do {
      do {
	radius = drand48();
	radius *= drand48();
      } while ((radius < 1.0e-20) && ((1.-radius) < 1.0e-20));
      radius = -log(radius) * Rd;
      TotalStars++;
    } while (radius > 0.25 * Ax * Nx);

    phi_orbit = 2 * M_PI * drand48();
    phi_epicyclic = 2* M_PI * drand48();
    x_stars[i] = radius * cos(phi_orbit);
    y_stars[i] = radius * sin(phi_orbit);
    velocity = CircularVelocity(radius);
    kappa = EpicyclicFrequency (radius);
    Omega = velocity / radius;
    MeanDispersion = DispersionProfile (radius);
    phi_offset = velocity / radius * DT / 2.;
    velocityt = velocity * velocity - radius / Rd * MeanDispersion * MeanDispersion;
    tempwork = kappa * kappa / 4. / Omega / Omega - 1.;
    velocityt -= MeanDispersion * MeanDispersion * tempwork;
    /*	DerivDisp = DispersionProfile (radius+1.)-DispersionProfile(radius); */
    DerivDisp = - DispersionProfile (radius) / Rd;
    velocityt += radius * 2. * DispersionProfile(radius) * DerivDisp;
    if (velocityt < 0.) {
      MeanDispersion  = DispersionProfile (radius);
      velocityt = 0.;
    } 
    velocity = sqrt(velocityt);
    do {
      do
	varalea = drand48();
      while (varalea <1e-20);
      dispersion = 1.08 * MeanDispersion * sqrt(-log(varalea)*2.);
    } while (dispersion > 2.5 * MeanDispersion);
    px_stars[i] = -velocity * REDUCTION * sin(phi_orbit - phi_offset);
    py_stars[i] = velocity * REDUCTION * cos(phi_orbit - phi_offset);
    px_stars[i] += dispersion * sin(phi_epicyclic)*cos(phi_orbit);
    px_stars[i] += dispersion * cos(phi_epicyclic) * sin(phi_orbit) * kappa/2./Omega;
    py_stars[i] += dispersion * sin(phi_epicyclic)*sin(phi_orbit);
    py_stars[i] -= dispersion * cos(phi_epicyclic)*cos(phi_orbit) * kappa/2./Omega;
  }
  stars->Mass = TheoMass () / TotalStars;
  /* The numerator represents the total mass of the galaxy */
  if (Symmetrized == YES) {
    stars->Mass /= 2.;
  }
  Message("done\n");
}

void DumpAngularVelocities()
{
  FILE           *dump;
  real            velocity, radius, Oort_A, kappa;
  real            dispersion, Omega, ILR, OLR;
  char 		filename[100];
  sprintf (filename, "%sangular.dat", OUTPUTDIR);
  Message("Writing 'angular.dat'...");
  dump = fopen(filename, "w");
  if (dump == NULL) {
    fprintf(stderr, "Unable to open 'angular.dat'\n");
    exit(1);
  }
  for (radius = XRESOL / 10.; radius < XRESOL * XSIZE * .25; radius += XRESOL / 5.) {
    velocity = CircularVelocity(radius);
    Oort_A = CircularVelocity(radius + 0.001) - velocity;
    Oort_A *= radius;
    Oort_A -= velocity / radius;
    Oort_A *= .5;
    kappa = Oort_A + velocity / radius;
    kappa *= 4. * velocity / radius;
    kappa = EpicyclicFrequency (radius);
    Omega = velocity / radius;
    ILR = Omega - kappa / 2.;
    OLR = Omega + kappa / 2.;
    dispersion = Q_Toomre(radius) * 3.36 * G * Sigma(radius);
    dispersion /=  kappa;
    fprintf(dump, "%f\t%f\t%f", radius, velocity, Omega);
    fprintf(dump, "\t%f\t%f\t%f\t", ILR, OLR, dispersion);
    fprintf(dump, "%f\t", dispersion * sqrt(1. + kappa * kappa / 4. / Omega / Omega));
    fprintf(dump, "%f\n", Q_Toomre(radius));
  }
  fclose(dump);
  Message("done\n");
}
