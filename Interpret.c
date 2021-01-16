/* 'sim2d' (c) F.Masset - 1995 - 2020. Codigo 2D de simulaciones
   N-cuerpos con metodo PIC y evaluacion del potencial por FFTs,
   enfocado a dinamica galactica. Codigo adaptado para fines de
   ensenyanza de la materia Dinamica y Estructura de Galaxias del
   Posgrado en Astrofisica de la UNAM. Para evitar problemas de
   formato, todos los comentarios en espanyol en las fuentes de este
   codigo son voluntariamente sin acentos */

#include "sim2d.h"
#define MAXVARIABLES 200

extern boolean  Symmetrized;
static Param    VariableSet[MAXVARIABLES];
static int      VariableIndex = 0;
static int	FirstStep = YES;
static clock_t  First, Preceeding, Current, FirstUser, CurrentUser, PreceedingUser;
static long	Ticks;

int		XSIZE;		
int		YSIZE;		
real		XRESOL;	
real		YRESOL;	
int		NUMBEROFSTARS;
real		ZSOFT;  
real		DT;
real		SIGMA0;	
real		Rd;
int		NINTERM;
real		TOOMRECENTRAL;
int		NTOT;
real 		REDUCTION;
int		NINTERV;
real		CENTRALMASS;
real		HALOCORE;
real		HALOSPEEDLIM;
real		BULGERADIUS;
real		BULGEMASS;
char		OUTPUTDIR[81];
char		SYMMETRY[81];

void InitVariables()
{
  Var("XSIZE", (char *)&XSIZE, INT, YES, "128.");
  Var("YSIZE", (char *)&YSIZE, INT, YES, "128.");
  Var("XRESOL", (char *)&XRESOL, REAL, YES, "800.");
  Var("YRESOL", (char *)&YRESOL, REAL, YES, "800.");
  Var("NUMBEROFSTARS", (char *)&NUMBEROFSTARS, INT, YES, "20000.");
  Var("ZSOFT", (char *)&ZSOFT, REAL, YES, "300.");
  Var("DT", (char *)&DT, REAL, YES, "1.0");
  Var("SIGMA0", (char *)&SIGMA0, REAL, YES, "173.");
  Var("RD", (char *)&Rd, REAL, YES, "3500.");
  Var("NINTERM", (char *)&NINTERM, INT, YES, "10.");
  Var("TOOMRECENTRAL", (char *)&TOOMRECENTRAL, REAL, NO, "1.5");
  Var("NTOT", (char *)&NTOT, INT, YES, "1501.");
  Var("REDUCTION", (char *)&REDUCTION, REAL, NO, "1.0");
  Var("NINTERV", (char *)&NINTERV, INT, YES, "500.");
  Var("CENTRALMASS", (char *)&CENTRALMASS, REAL, NO, "0.");
  Var("HALOCORE", (char *)&HALOCORE, REAL, NO, "2000.");
  Var("HALOSPEEDLIM", (char *)&HALOSPEEDLIM, REAL, NO, "0.");
  Var("BULGERADIUS", (char *)&BULGERADIUS, REAL, NO, "1000.");
  Var("BULGEMASS", (char *)&BULGEMASS, REAL, NO, "0.");
  Var("OUTPUTDIR", OUTPUTDIR, STRING, NO, "out");
  Var("SYMMETRY", SYMMETRY, STRING, NO, "NO");
}

void Var(char *name, char *ptr, int type, int necessary, char *deflt)
{
  real            valuer;
  int             valuei;
  float		  temp;
  sscanf (deflt, "%f", &temp);
  valuer = (real) (temp);
  valuei = (int) valuer;
  strcpy(VariableSet[VariableIndex].name, name);
  VariableSet[VariableIndex].variable = ptr;
  VariableSet[VariableIndex].type = type;
  VariableSet[VariableIndex].necessary = necessary;
  VariableSet[VariableIndex].read = NO;
  if (necessary == NO) {
    if (type == INT) {
      *((int *) ptr) = valuei;
    } else if (type == REAL) {
      *((real *) ptr) = valuer;
    } else if (type == STRING) {
      strcpy (ptr, deflt);
    }
  }
  VariableIndex++;
}

void ReadVariables(char *filename)
{
  char            nm[80], s[200],stringval[81];
  char           *s1;
  float           temp;
  real            valuer;
  int             i, found, valuei, success, type;
  int            *ptri;
  real           *ptrr;
  FILE           *input;

  InitVariables();
  input = fopen(filename, "r");
  if (input == NULL) {
    fprintf(stderr, "Unable to read '%s'. Program stopped.\n",filename);
    exit(1);
  }
  fprintf (stderr, "Reading parameters file '%s'.\n", filename);
  while (fgets(s, 199, input) != NULL) {
    success = sscanf(s, "%s ", nm);
    if ((nm[0] != '#') && (success == 1)) {	/* # begins a comment
						 * line */
      s1 = s + strlen(nm);
      sscanf(s1 + strspn(s1, "\t :=>_"), "%f", &temp);
      sscanf(s1 + strspn(s1, "\t :=>_"), "%80s ", stringval);
      valuer = (real) temp;
      valuei = (int) temp;
      for (i = 0; i < (int)strlen(nm); i++) {
	nm[i] = (char) toupper(nm[i]);
      }
      found = NO;
      for (i = 0; i < VariableIndex; i++) {
	if (strcmp(nm, VariableSet[i].name) == 0) {
	  if (VariableSet[i].read == YES) {
	    fprintf(stderr, "Warning : %s defined more than once.\n", nm);
	  }
	  found = YES;
	  VariableSet[i].read = YES;
	  ptri = (int *) (VariableSet[i].variable);
	  ptrr = (real *) (VariableSet[i].variable);
	  if (VariableSet[i].type == INT) {
	    *ptri = valuei;
	  } else if (VariableSet[i].type == REAL) {
	    *ptrr = valuer;
	  } else if (VariableSet[i].type == STRING) {
	    strcpy (VariableSet[i].variable, stringval);
	  }
	}
      }
      if (found == NO) {
	fprintf(stderr, "Warning : variable %s defined but non-existent in code.\n", nm);
      }
    }
  }

  found = NO;
  for (i = 0; i < VariableIndex; i++) {
    if ((VariableSet[i].read == NO) && (VariableSet[i].necessary == YES)) {
      if (found == NO) {
	fprintf(stderr, "Fatal error : undefined necessary variables:\n");
	found = YES;
      }
      fprintf(stderr, "%s\n", VariableSet[i].name);
    }
    if (found == YES)
      exit(1);

  }
  found = NO;
  for (i = 0; i < VariableIndex; i++) {
    if (VariableSet[i].read == NO) {
      if (found == NO) {
	fprintf(stderr, "Secondary variables omitted :\n");
	found = YES;
      }
      if ((type = VariableSet[i].type) == REAL)
	fprintf(stderr, "%s ;\t Default Value : %.5g\n", VariableSet[i].name, *((real *) VariableSet[i].variable));
      if (type == INT)
	fprintf(stderr, "%s ;\t Default Value : %d\n", VariableSet[i].name, *((int *) VariableSet[i].variable));
      if (type == STRING)
	fprintf(stderr, "%s ;\t Default Value : %s\n", VariableSet[i].name, VariableSet[i].variable);
    }
  }
  if (*SYMMETRY == 'Y') 	Symmetrized = YES;
  else			Symmetrized = NO;
}

void PrintUsage (char *execname)
{
  fprintf (stderr, "Usage : %s [-nvt] <parameters file>\n", execname);
  fprintf (stderr, "\n-n : Disables simulation. The program just read parameters file\n");
  fprintf (stderr, "-v : Verbose mode. Tells everything about parameters file\n");
  fprintf (stderr, "-t : Gives time information at each time step\n");
  exit (EXIT_FAILURE);
}

void TellEverything () {
  real masshalo;
  printf("\n*** GRID CHARACTERISTICS:\n");
  printf("%d * %d cells, each %f * %f parsecs wide\n",XSIZE/2,YSIZE/2,XRESOL,YRESOL);
  printf("This gives %f * %f kpc for the total size\n",XSIZE*XRESOL/2000,YSIZE*YRESOL/2000);
  printf("The maximum radius is %f kpc.\n",(real)XSIZE*XRESOL/4000.);
  printf("The softened gravity scale is %f pc.\n", ZSOFT);
  printf("(Note that the actual grid is twice as big in each dimension\n");
  printf("but only the fourth is used to avoid aliasing in FFTs)\n");
  printf("\n*** OUTPUTS CHARACTERISTICS:\n");
  printf("Time step : %f (in Million Years)\n", DT);
  printf("Total number of time steps : %d\n", NTOT);
  printf("Time steps between outputs : %d\n", NINTERM);
  printf("Total time elapsed : %f (in Million Years)\n", (real)NTOT*DT);
  printf("Total number of outputs : %d\n", NTOT/NINTERM);
  printf("Outputs directory : %s\n", OUTPUTDIR);
  printf("\n*** GALAXY CHARACTERISTICS:\n");
  printf("Total mass of stars : %.2g Solar Masses\n", TheoMass ());
  printf ("Exponential Radius of disk : %f kpc\n", Rd/1000.);
  printf ("Surfacic density at GC: %f Sol. Mas. per pc^2\n", SIGMA0);
  printf ("Total number of Stars : %d\n", (Symmetrized==YES) ? 2 * NUMBEROFSTARS : NUMBEROFSTARS);
  printf ("Distribution %ssymmetrized\n", (Symmetrized == YES) ? "":"non-");
  printf ("Mass per 'star' : %g Solar Masses\n", TheoMass()/(real)NUMBEROFSTARS*((Symmetrized == YES) ? .5 : 1.0));
  printf ("Central Toomre Parameter : %f\n", TOOMRECENTRAL);
  if (REDUCTION != 0.0) {
    printf ("Reduction of circular velocity : %f\n", REDUCTION);
  }
  printf ("Number of intervals for Hankel integrals computation : %d\n", NINTERV);

  if (HALOSPEEDLIM != 0.0) {
    printf ("\n*** HALO CHARACTERISTICS:\n");
    printf ("Limit circular velocity : %f\n", HALOSPEEDLIM);
    printf ("Core radius : %f kpc\n", HALOCORE / 1000.);
    printf ("Mass (inside optical boundaries) : %g Solar Masses\n", masshalo = HALOCORE*HALOSPEEDLIM*HALOSPEEDLIM/G*LHalo(XRESOL*XSIZE/4./HALOCORE));
    printf ("Ratio Halo/Disk : %f\n", masshalo / TheoMass());
  } else printf ("\n*** NO HALO\n");

  if (BULGEMASS != 0.0) {
    printf ("\n*** BULGE CHARACTERISTICS:\n");
    printf ("Mass : %g Solar Masses\n", BULGEMASS);
    printf ("Mass ratio Bulge/Disk : %f\n", BULGEMASS/TheoMass());
    printf ("Bulge radius : %f\n", BULGERADIUS);
  } else printf ("\n*** NO BULGE\n");

  if (CENTRALMASS != 0.0) {
    printf ("\n*** CENTRAL OBJECT CHARACTERISTICS:\n");
    printf ("Mass : %g Solar Masses\n", CENTRALMASS);
    printf ("Mass ratio Center/Disk : %f\n", CENTRALMASS/TheoMass());
  } else printf ("\n*** NO CENTRAL OBJECT\n");
  if (Symmetrized == YES) {
    printf ("\n*** CENTRALLY SYMMETRIC DISTRIBUTION\n");
  } else {
    printf ("\n*** NOT CENTRALLY SYMMETRIC DISTRIBUTION\n");
  }
}

void GiveTimeInfo (int number)
{
  struct tms buffer;
  real total, last, mean, totalu;
  Current = times (&buffer);
  CurrentUser = buffer.tms_utime;
  if (FirstStep == YES) {
    First = Current;
    FirstUser = CurrentUser;
    fprintf (stderr, "Time counters initialized\n");
    FirstStep = NO;
    Ticks = sysconf (_SC_CLK_TCK);
  }
  else {
    total = (real)(Current - First)/Ticks;
    totalu= (real)(CurrentUser-FirstUser)/Ticks;
    last  = (real)(CurrentUser - PreceedingUser)/Ticks;
    mean  = totalu / number;
    fprintf (stderr, "Total Real Time elapsed    : %.3f s\n", total);
    fprintf (stderr, "Total CPU Time of process  : %.3f s (%.1f %%)\n", totalu, 100.*totalu/total);
    fprintf (stderr, "CPU Time since last time step : %.3f s\n", last);
    fprintf (stderr, "Mean CPU Time between time steps : %.3f s\n", mean);
    fprintf (stderr, "CPU Load on last time step : %.1f %% \n", (real)(CurrentUser-PreceedingUser)/(real)(Current-Preceeding)*100.);

  }	
  PreceedingUser = CurrentUser;
  Preceeding = Current;
}
