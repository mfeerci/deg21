#include <stdio.h>
#include <math.h>
#include <string.h>
#include <sys/times.h>
#include <unistd.h>
#include <stdlib.h>
#include <ctype.h>
#include "types.h"
#include "param.h"


/* Constante de la gravedad en las unidades del codigo, que son: 
        - masa solar 
        - parsec 
        - million de anyos */

#define		G	4.4988e-3 //Msol^-1 pc^3 Myr^-2



/* Prototipos de todas las funciones del codigo */

void            AdvanceParticles(ParticleSet*, real);
real 		AdvanceVelocities (ParticleSet*, Grid*, real, boolean);
void*           AllocWorkArrays (void);
real            Bessj1(real);
real            Bessj0(real);
real            Bessi0(real);
real            Bessi1(real);
real            Bessk0(real);
real            Bessk1(real);
real            CircularVelocity(real);
real            ComputeS(real);
void            CopyGrid (Grid*, Grid*);
Grid*           CreateGrid(int,real,int,real,char*);
ParticleSet*    CreateParticleSet(int,real,char*);
real            DispersionProfile (real);
void            DumpAngularVelocities(void);
void            DumpParticles(ParticleSet *, int);
void            DumpSurfaceDensity(ParticleSet *, int);
real            EpicyclicFrequency(real);
void            ErrorMessage(char *);
void            FillArrayS(void);
void            FillCircularVelocity(void);
void            ForgetGrid (Grid *);
void            FourierDirecte(Grid *);
void            FourierInverse(Grid *);
void            Fourn(float*,int*,int,int);
void            FreeWorkArrays (void);
Pair            GetForce(real,real,real,real,int,int,real*,boolean);
void            GiveTimeInfo (int);
void            InitStars (ParticleSet*);
void            InitVariables(void);
void            InitKernelPotential(Grid*);
real            Interpole (Grid*,real,real);
Pair            KeplerianCatch(real,real);
real            LHalo (real);
int             main(int, char **);
void            Message (char *);
void            Multiply(Grid *, Grid *, real);
void            ParticleInCell(Grid*, ParticleSet*);
void            PrepareExtForce (void);
void            PrintUsage (char*);
real            Q_Toomre (real);
void            ReadVariables(char*);
real            Sigma(real);
real            StaticCircularVelocity(real);
real            StaticEpicyclicFrequency(real);
void            TellEverything (void);
real 		TheoMass (void);
real            TotalMass(ParticleSet*);
void            Var(char*,char*,int,int,char*);
void            WriteDim (void);
void            WriteDisk(Grid *, int);
////////////////////////////////////
// boolean         Symmetrized = YES;
