OBJ= main.o Advance.o Init.o Interpret.o Algo.o Kernels.o\
     Theo.o Output.o LowTasks.o NR_Fourier.o NR_Bessel.o 

LIBS=-lm #Libreria para compilar el codigo (m : maths)

OPT=-O3 -Wall -Wextra # Opciones para compilar el codigo (O3:		\
optimizacion agresiva y todos los warnings -advertencias- posibles)
#OPT= -g    #Esta opcion se usaria para depurar

NAME= sim2d    # Nombre el archivo ejecutable

$(NAME): $(OBJ)
	gcc $(OBJ) $(OPT) -o $(NAME) $(LIBS)

$(OBJ): param.h types.h Makefile
.c.o:
	gcc $*.c  $(OPT) -c

clean:
	/bin/rm -f $(OBJ) $(NAME)
	/bin/rm -f *~
