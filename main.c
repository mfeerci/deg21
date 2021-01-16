/* 'sim2d' (c) F.Masset - 1995 - 2020. Codigo 2D de simulaciones
   N-cuerpos con metodo PIC y evaluacion del potencial por FFTs,
   enfocado a dinamica galactica. Codigo adaptado para fines de
   ensenyanza de la materia Dinamica y Estructura de Galaxias del
   Posgrado en Astrofisica de la UNAM. Para evitar problemas de
   formato, todos los comentarios en espanyol en las fuentes de este
   codigo son voluntariamente sin acentos */

#include "sim2d.h"

boolean         Symmetrized = YES;

int main(int argc, char *argv[])
{
  ParticleSet    *stars;
  Grid           *star_density;
  Grid           *kernel;
  Grid 	         *potential;
  real		  outstar;
  FILE	         *ControlEnergy;	
  int             i, verbose = NO, disable = NO, TimeInfo = NO, temp;
  char		  EnergyFile[100], ParameterFile[100], command[512];
  char            msg[512];

  if (argc == 1) PrintUsage (argv[0]);
  strcpy (ParameterFile, "");
  for (i = 1; i < argc; i++) {
    if (*(argv[i]) == '-') {
      if (strchr (argv[i], 'n'))
	disable = YES;
      if (strchr (argv[i], 'v'))
	verbose = YES;
      if (strchr (argv[i], 't'))
	TimeInfo = YES;
      if (strspn (argv[i], "-nvt") != strlen (argv[i]))
	PrintUsage (argv[0]);
    }
    else strcpy (ParameterFile, argv[i]);
  }
  if (ParameterFile[0] == 0) PrintUsage (argv[0]);
	
  ReadVariables (ParameterFile);

  if (verbose == YES) 
    TellEverything ();
  if (disable == YES)
    exit (0);

  stars         = CreateParticleSet(NUMBEROFSTARS, 0., "stars");
  kernel        = CreateGrid(XSIZE, XRESOL, YSIZE, YRESOL, "kernel");
  star_density  = CreateGrid(XSIZE, XRESOL, YSIZE, YRESOL, "stardens");
  potential     = CreateGrid(XSIZE, XRESOL, YSIZE, YRESOL, "potential");

  sprintf (command, "mkdir -p %s", OUTPUTDIR);
  temp = system (command);
 
  sprintf (EnergyFile, "%senergy.dat", OUTPUTDIR);
  ControlEnergy = fopen (EnergyFile, "w");
  if (ControlEnergy == NULL) {
    fprintf (stdout, "Unable to open 'energy.dat'.\nProgram stopped.\n");
    fprintf (stdout, "Line OUTPUTDIR in parameters file may be incorrect.\n");
    exit (1);
  }
  fclose (ControlEnergy);

  if (AllocWorkArrays () == NULL) {
    fprintf (stderr, "Insufficient memory. Program stopped.\n");
    exit (1);
  }
  WriteDim ();
  FillArrayS ();
  FillCircularVelocity ();
  DumpAngularVelocities ();	
  InitStars(stars);
  InitKernelPotential(kernel);
  PrepareExtForce ();
  
  fprintf (stdout, "\nINTEGRATION STARTS\n\n");
  
  for (i = 0; i < NTOT; i++) {
    
    ParticleInCell(star_density, stars);

    /*
     * The lines below correspond to the computation of the
     * gravitational potential due to the stars
     */

    CopyGrid (star_density, potential);
    FourierDirecte(potential);

    /* POSGRADO: la siguiente linea multiplica 'potential' (que a
       estas alturas es la FFT de la densidad de superficie, vease
       lineas anteriores) por 'kernel', que es la FFT del kernel de
       convolucion. La meta es realizar una convolucion, que es una
       mera multiplicacion en el espacio de Fourier. Despues de llamar
       a esta rutina, la malla potential habr'a sido multiplicada por
       la malla kernel, y entonces s'i contendr'a la FFT del potencial
       galactico. Sin embargo, es indispensable normalizar
       correctamente este producto. Por eso, el ultimo y tercer
       argumento de la funcion 'Multiply' tiene que ser una variable
       flotante (y no una malla como los dos primeros), que se usa
       como factor adicional en la multiplicacion de las dos
       mallas. Se tiene que encontrar su valor correcto, a sabiendas
       de que depende solamente de XSIZE, XRESOL, YSIZE y YRESOL. Si
       estan atorados con la teoria, podran proceder con ensayo y
       error y tratar de reproducir la magnitud del potencial a t=0
       que aparece en el jupyter notebook 'PresentacionSIM2D'. Para
       ello tendran que haber resuelto la pregunta de 'Algo.c', ya que
       se usa la funcion 'ParticleInCell' mas arriba. */
        
    Multiply(potential, kernel,  XRESOL*YRESOL/XSIZE/YSIZE);
    FourierInverse(potential);

    if (i <= (NTOT/NINTERM)*NINTERM)
      sprintf (msg, "%d", i/NINTERM);
    else
      sprintf (msg, "N/A");

    fprintf (stdout,"\r--- TIMESTEP %d --- NEXT OUTPUT: %s ---", i, msg);
    fflush  (stdout);
    if (TimeInfo == YES) {
      printf ("\n");
      GiveTimeInfo (i/NINTERM);
    }
    if (NINTERM * (i / NINTERM) == i) {	/* Outputs if i multiple of NINTERM */ 
      printf ("\n\n");
      WriteDisk(star_density, i / NINTERM); 
      WriteDisk(potential, i / NINTERM);
      DumpParticles (stars, i / NINTERM);
      printf ("\n");
    }

    /* We advance the stars velocities */

    outstar = AdvanceVelocities (stars, potential, DT, YES);

    /* We advance the stars positions */

    AdvanceParticles (stars, DT);
  }
  FreeWorkArrays ();
  return EXIT_SUCCESS;
}
