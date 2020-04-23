//
// Suma los cuadrados de los n primeros numeros naturales
//
// Cada proceso calcula un cuadrado, el primer proceso recibe todos y los suma
//

#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<float.h>
#include<stdbool.h>
#include<cputils.h>
#include<mpi.h>
#include<stddef.h>

typedef struct {
	float pos_row, pos_col;		// Position
	float mov_row, mov_col;		// Direction of movement
	float choose_mov[3];		// Genes: Probabilities of 0 turning-left; 1 advance; 2 turning-right
	float storage;			// Food/Energy stored
	int age;			// Number of steps that the cell has been alive
	unsigned short random_seq[3];	// Status value of its particular random sequence
	bool alive;			// Flag indicating if the cell is still alive
} Cell;

int main(int argc, char *argv[]) {
	float i = 24.6;
	Cell *cells = (Cell *)malloc( sizeof(Cell) * (size_t)10 );

	for (int s=0 ; s<10; s++){

		x=cells[i].alive;
		 printf(x ? "true" : "false");
	}
}
