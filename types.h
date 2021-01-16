/* Aqui se definen los tipos usados en el codigo. Algunos solamente
   son sinonimos (real, boolean), otros son mas sofisticados (por
   ejemplo Grid, que embebe toda los datos necesarios para definir una
   malla, o ParticleSet, que contiene toda la informacion relativa a
   las estrellas: x, y, vx, vy, asi como el numero de estas y su
   masa. */

typedef int     boolean;
typedef float	real;

struct pair {
  real            x;
  real            y;
};

typedef struct pair Pair;

struct grid {		/**< Estructura describiendo una malla */
  int             Nx;        	/**< Numero de celdas en x */
  int             Ax;		/**< Resolucion (tamanyo de UNA celda) en x */
  int             Ny;		/**< Numero de celdas en y */
  int             Ay;		/**< Resolucion (tamanyo en UNA celda) en y */
  real           *Field;	/**< Puntero al campo almacenado en la
				   malla (p.e. densidad de
				   superficie) */
  char           *Name;		/**< Un nombre que se da a la malla, y
				   que se usa de manera automatica
				   cuando se guarda en un archivo
				   (usando Output.c::WriteDisk()) */
};

typedef struct grid Grid;	/**< Podremos usar el tipo "Grid" en
				   lugar de la version larga mas
				   tediosa "struct grid" */

struct particle_set {		/**< Estructura describiendo el conjunto de estrellas  */
  int             NumberParticles;   /**< Numero de estrellas */
  real            Mass;		     /**< Masa de las estrellas (la misma para todas: no es un arreglo)  */
  real           *Xparticles;	     /**< Arreglo (puntero) de la coordenada X de las estrellas */
  real           *Yparticles;        /**< Arreglo (puntero) de la coordenada Y de las estrellas */
  real           *VXparticles;	     /**< Arreglo (puntero) de la
					coordenada X de la velocidad
					de las estrellas */
  real           *VYparticles;	/**< Arreglo (puntero) de la
				   coordenada Y de la velocidad de las
				   estrellas */
  char           *Name;		/**< Nombre de la estructura, que se
				   usa de manera automatica cuando se
				   guarda en un archivo (usando
				   Output.c::DumpParticles()) */
};

typedef struct particle_set ParticleSet; /**< Podemos usar el alias breve "ParticleSet" para el tipo  */

/* Algunas definiciones */

#define		YES	1
#define		NO	0
#define		REAL	1
#define		INT	0
#define		STRING  2

/* La estructura siguiente se usa solamente para leer el archivo de
   parametros antes de arrancar la simulacion. Entender su
   funcionamiento no es imprescindible para entender la parte
   astrofisica del codigo */

struct param {			/**< Estructura para leer parametros */
  char name[80];		/**< Nombre del parametro (p.e. XRESOL) */
  int  type;			/**< Tipo del parametro (REAL, INT o STRING) */
  char *variable;		/**< Puntero a la variable
				   correspondiente. Se cambiara el
				   tipo del puntero segun el tipo al
				   momento de la lectura */
  int read;			/**< Para saber si ya se leyo el
				   parametro en el archivo. Se usa
				   para imprimir un mensaje de error
				   si se define por error dos veces el
				   mismo parametro en un archivo de
				   parametros. */
  int necessary;		/**< Para saber si es imperativo que
				   el usuario defina el parametro. Si
				   no es el caso y no esta definido,
				   se usa un valor por default
				   definido en la funcion
				   Interpret.c::InitVariables() */
};

typedef struct param Param;	/**< Se puede usar el alias "Param" para el tipo  */
