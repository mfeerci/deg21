/* 'sim2d' (c) F.Masset - 1995 - 2020. Codigo 2D de simulaciones
   N-cuerpos con metodo PIC y evaluacion del potencial por FFTs,
   enfocado a dinamica galactica. Codigo adaptado para fines de
   ensenyanza de la materia Dinamica y Estructura de Galaxias del
   Posgrado en Astrofisica de la UNAM. Para evitar problemas de
   formato, todos los comentarios en espanyol en las fuentes de este
   codigo son voluntariamente sin acentos */

#include "sim2d.h"

#define  MAXSTEPS 20000 	/* Resolution of External Force Array */

extern boolean  Symmetrized;

static real	ExtForceStep;
static real	InvExtForceStep;
static real	ExtForce[MAXSTEPS];

void PrepareExtForce () {
  int i;
  real radius, sqradius, central, CENTRALtemp1, ZSOFT2;
  real HALOSPEEDtemp1, b, a, bulge, halo;
  ExtForceStep = XRESOL * XSIZE / 2. / MAXSTEPS;
  InvExtForceStep = 1./ExtForceStep;
  CENTRALtemp1 = G * CENTRALMASS;
  HALOSPEEDtemp1 = HALOSPEEDLIM * HALOSPEEDLIM * HALOCORE;
  ZSOFT2 = ZSOFT*ZSOFT;
  fprintf (stdout, "Initializing external forces...");

  for (i = 1; i < MAXSTEPS; i++) {
    radius = (real)i * ExtForceStep;
    sqradius = radius * radius;
    central = CENTRALtemp1 / (sqradius + ZSOFT2);
    halo = HALOSPEEDtemp1/sqradius;
    halo *= LHalo(radius/HALOCORE);
    b = BULGERADIUS;
    a = sqrt(b * b + sqradius);
    bulge = G * BULGEMASS * radius / ((b + a) * (b + a) * a);
    central += halo + bulge;
    ExtForce[i] = central / (radius);
  }
  ExtForce[0] = 0.;
  fprintf (stdout, "done.\n");
}

void ParticleInCell(Grid *array, ParticleSet *particles)
{
  int   k, l, i, j, N2x, N2y, Nx, Ny, Number;
  real  resx, resy, x, y, sx, sy, xr, yr, quantum, sqradius;
  real  a1, a2, a3, a4, ParticleMass, sqradiusmax;
  real *density, *x_particles, *y_particles, Ax, Ay;
  real  invAx, invAy;
  
  ParticleMass = particles->Mass;
  Nx = array->Nx;
  Ny = array->Ny;
  density = array->Field;
  Ax = array->Ax;
  Ay = array->Ay;
  invAx = 1./Ax;
  invAy = 1./Ay;
  Number = particles->NumberParticles;
  x_particles = particles->Xparticles;
  y_particles = particles->Yparticles;

  quantum = ParticleMass / Ax / Ay;
  for (k = 0; k < (Nx + 1) * (Ny + 1) * 2; k++)
    density[k] = 0.;
  N2x = Nx / 2;
  N2y = Ny / 2;
  sqradiusmax = .0625 * Ax * Ax * Nx * Nx;
  for (k = 0; k < Number; k++) {
    x = x_particles[k];
    y = y_particles[k];
    xr = x * invAx;
    yr = y * invAy;
    sx = floor(xr);
    sy = floor(yr);
    i = N2x + (int) sx;
    j = N2y + (int) sy;
    resx = xr - sx;
    resy = yr - sy;
    sqradius = x * x + y * y;
    if (sqradius < sqradiusmax) {
      l = 2 * ((i - 1) * Nx + j) - 1;
      a1 = a2 = quantum * (1. - resy);
      a1 *= 1. - resx;
      a2 *= resx;
      a3 = a4 = quantum * resy;
      a3 *= 1. - resx;
      a4 *= resx;
      /* POSGRADO: en las 4 lineas siguientes, se afectan 4 valores
	 (a1, a2, a3 y a4) a cuatro celdas vecinas de la malla de
	 densidad, usando un metodo PIC. Son 4 y no 2 como lo hemos
	 visto en clase, ya que aqui el problema es en 2 dimensiones
	 en lugar de una dimension.  Tal como lo hemos visto en clase,
	 a las coordenadas (i,j) corresponde el indice 'l'. A las
	 coordenadas (i+1,j), 'l+2'. A las coordenadas (i,j+1),
	 'l+2*Nx'. Que indice corresponde a (i+1,j+1) ? Escribalo en
	 lugar del comentario en la linea donde se usa 'a4' */
      density[l] += a1;
      density[l + 2 * Nx] += a2;
      density[l + 2] += a3;
      density[ l + 2 + 2*Nx  ] += a4;


      /*
       * We fill the grid according to the first-order
       * method PIC (Particle in Cell)
       */
      if (Symmetrized == YES) {
	i = Nx - i - 1;
	j = Ny - j - 1;
	/*
	 * The two equations above are valid only if
	 * x and y are not exactly integers. However,
	 * the probability that they are exactly
	 * integers is almost null, hence we can use
	 * these rapid formulae for symmetrization
	 */
	l = 2 * ((i - 1) * Nx + j) - 1;
	density[l] += a4;
	density[l + 2 * Nx] += a3;
	density[l + 2] += a2;
	density[l + 2 * Nx + 2] += a1;
      }
    }
  }
}


Pair GetForce(real x, real y, real invAx, real invAy,\
	       int Nx, int Ny, real *pot, boolean gravit)
{
  real		phi00, phi01, phi02, phi10, phi11, phi12, phi20, phi21, phi22;
  real            sx, sy, xr, yr;
  real		central=0.0, radius, sqradius;
  real		temp2x, temp2y, x1my, y1mx, crosscr, crossr;
  Pair            dphi;
  int             i, j, l;
  real		integradius, resi, resj;
  int		ii, jj;
  register real 	resx;
  register real	resy;
  register real	dxphi;
  register real	dyphi;
  register real	cresx;
  register real	cresy;

  if (gravit == YES) {
    sqradius = x * x + y * y;
    radius = sqrt(sqradius);
    integradius = radius * InvExtForceStep;
    ii = integradius;
    jj = ii + 1;
    resi = integradius - ii;
    resj = 1.-resi;
    central = ExtForce[ii]*resj+ExtForce[jj]*resi;
  }
  xr = x * invAx;
  yr = y * invAy;
  sx = floor(xr);
  sy = floor(yr);
  resx = xr - sx;
  resy = yr - sy;
  cresx = 1.-resx;
  cresy = 1.-resy;
  i = (Nx>>1) + (int)sx;
  j = (Ny>>1) + (int)sy;
  l = 2 * (Ny * (i - 1) + j) - 1;

  phi11 = pot[l];
  phi21 = pot[l + 2*Nx];
  phi22 = pot[l + 2 + 2 * Nx];
  phi12 = pot[l + 2];
  phi01 = pot[l - 2*Nx];
  phi20 = pot[l - 2 + 2 * Nx];
  phi10 = pot[l - 2];
  phi00 = pot[l - 2 - 2 * Nx];
  phi02 = pot[l + 2 - 2 * Nx];

  temp2x = 1. - 2.*resx;
  temp2y = 1. - 2.*resy;
  x1my   = resx * cresy;
  y1mx   = resy * cresx;
  crosscr= cresx * cresy;
  crossr = resx * resy; 

  dxphi = phi11*temp2x+phi10*resy*temp2x+phi12*cresy*temp2x-phi01*cresx+phi21*resx;
  dxphi+= phi22*x1my + phi20*crossr - phi02*crosscr - phi00*y1mx;

  dyphi = phi11*temp2y+phi01*resx*temp2y+phi21*cresx*temp2y-phi10*cresy+phi12*resy;
  dyphi+= phi22*y1mx + phi02*crossr - phi20*crosscr - phi00*x1my;

  dphi.x = - dxphi * 0.5 * invAx;
  dphi.y = - dyphi * 0.5 * invAy;
  if ((gravit == YES) && (central != 0.)) {
    dphi.x -= central * x;
    dphi.y -= central * y;
  }
  return dphi;
}

real Interpole (Grid *array, real x, real y)
{
  real ii, jj, resolx, resoly, p1, p2, p3, p4;
  real p, resx, resy, *val;
  int Nx, Ny, i, j, l;
  resolx = array->Ax;
  resoly = array->Ay;
  Nx = array->Nx;
  Ny = array->Ny;
  val = array->Field;
  ii = x / resolx + Nx/2;
  jj = y / resoly + Ny/2;
  i = (int)floor(ii);
  j = (int)floor(jj);
  l = 2 * (Ny * (i - 1) + j) - 1;
  resx = ii-(real)i;
  resy = jj-(real)j;
  p1 = val[l];
  p2 = val[l+2*Nx];
  p3 = val[l+2];
  p4 = val[l+2*Nx+2];
  p = p1*(1.-resy)*(1.-resx)+p2*resx*(1.-resy)+p3*resy*(1.-resx);
  p += resx*resy*p4;
  return p;
}
	
Pair KeplerianCatch(real x, real y)
{
  real            sqradius, central, radius;
  Pair            dphi;
  int 		ii, jj;
  real 		integradius, resi, resj;
  
  sqradius = x * x + y * y;
  radius = sqrt(sqradius);
  integradius = radius * InvExtForceStep;
  ii = (int)integradius;
  jj = ii + 1;
  if (jj < MAXSTEPS-2) {
    resi = integradius - ii;
    resj = 1.-resi;
    central = ExtForce[ii]*resj+ExtForce[jj]*resi;
  } else {
    central = ExtForce[MAXSTEPS-2]*(real)(MAXSTEPS-2)/InvExtForceStep/radius;
  }	
  dphi.x = -central * x;
  dphi.y = -central * y;
  return dphi;
}

real AdvanceVelocities(ParticleSet *set, Grid *potential, real dt, boolean gravit)
{
  int  i, outside, Nx, Ny, NbStars;
  Pair dphi;
  real sqradius, mass, sqradiusmax, Ax, *x, *y, *vx, *vy;
  real OutsidePercentage, TheoreticalMass, Ay, invAx, invAy;
  real *pot;
  x = set->Xparticles;
  y = set->Yparticles;
  vx = set->VXparticles;
  vy = set->VYparticles;
  mass = set->Mass;
  NbStars = set->NumberParticles;
  Nx = potential->Nx;
  Ny = potential->Ny;
  Ay = potential->Ay;
  Ax = potential->Ax;
  invAx = 1./Ax;
  invAy = 1./Ay;
  pot = potential->Field;
  outside = 0;
  sqradiusmax = .0625 * Nx * Nx * Ax * Ax;
  TheoreticalMass = TheoMass ();
  /* Valid only if square grid */
  for (i = 0; i < NbStars; i++) {
    sqradius = x[i] * x[i] + y[i] * y[i];

    if (sqradius < sqradiusmax) {
      /*
       * The star is in the greatest disk inside the active
       * part of the grid
       */
      dphi = GetForce(x[i], y[i], invAx, invAy, Nx, Ny, pot, gravit);
      vx[i] += dt * dphi.x;
      vy[i] += dt * dphi.y;
    } else {
      /*
       * The star is not in the greatest disk inside the
       * active part of the grid
       */
      outside++;
      if (gravit == YES) {
	dphi = KeplerianCatch(x[i], y[i]);
	vx[i] += dt * dphi.x;
	vy[i] += dt * dphi.y;
      }
    }

  }
  OutsidePercentage = (real) outside / (real) NbStars;
  return OutsidePercentage;
}
