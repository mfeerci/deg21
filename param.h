/* Caracteristica de la malla: XSIZE y YSIZE son su numero de celdas
   en X y Y */
int		XSIZE;		
int		YSIZE;
/* XRESOL y YRESOL son el tamanyo de una celda en X y Y (en parsec) */
real		XRESOL;	
real		YRESOL;

/* Numero total de estrellas */
int		NUMBEROFSTARS;

/* Parametro de suavizado del potencial (la "epsilon" de la clase) */
real		ZSOFT;

/* Caracteristicas del disco de estrellas */
/* \Sigma(r) = SIGMA0 * exp (-radius/Rd) */
real		SIGMA0;	
real		Rd;

/* DT: paso de tiempo elemental */
real		DT;
/* NINTERM: numero de pasos de tiempo entre dos salidas (escrituras al disco) */
int		NINTERM;
/* NTOT: numero total de pasos de tiempo (despues la simulacion se detiene) */
int		NTOT;
/* Carpeta ("dir": directory) donde se escriben los datos de la simulacion */
char		OUTPUTDIR[81];

real		TOOMRECENTRAL;	/* Parametro de Toomre del disco */
real 		REDUCTION;	/* Reduccion de la velocidad orbital
				   de las estrellas para tomar en
				   cuenta la deriva asimetrica */
int		NINTERV;	/* Resolucion de la transformada de
				   Hankel que se usa para inicializar
				   el disco */
real		CENTRALMASS;	/* Masa del hipotetico cuerpo central
				   puntual (hoyo negro supermasivo) */
real		HALOCORE;	/* Radio del nucleo del halo */
real		HALOSPEEDLIM;	/* Velocidad orbital asintotica a gran
				   distancia del potencial del halo */
real		BULGERADIUS;	/* Radio del bulbo */
real		BULGEMASS;	/* Masa del bulbo */
char		SYMMETRY[81];	/* Si el usuario activa el parametro
				   SYMMETRY, a cada estrella en (x,y)
				   con velocidad (vx,vy) corresponde
				   tambien (en la rutina de PIC) otra
				   estrella en (-x,-y) con velocidad
				   (-vx,-vy). El numero efectivo de
				   estrellas es por ende
				   2*NUMBEROFSTARS, y por razones de
				   simetria la galaxia solamente puede
				   desarrollar un numero par de brazos
				   espirales */
