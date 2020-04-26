/*
 * Simplified simulation of life evolution
 *
 * Computacion Paralela, Grado en Informatica (Universidad de Valladolid)
 * 2019/2020
 *
 * v1.2
 *
 * (c) 2020 Arturo Gonzalez Escribano
 */
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<float.h>
#include<stdbool.h>
#include<cputils.h>
#include<mpi.h>
#include<stddef.h>

/* Structure to store data of a cell */
typedef struct {
	float pos_row, pos_col;		// Position
	float mov_row, mov_col;		// Direction of movement
	float choose_mov[3];		// Genes: Probabilities of 0 turning-left; 1 advance; 2 turning-right
	float storage;			// Food/Energy stored
	int age;			// Number of steps that the cell has been alive
	unsigned short random_seq[3];	// Status value of its particular random sequence
	bool alive;			// Flag indicating if the cell is still alive
} Cell;


/* Structure for simulation statistics */
typedef struct {
	int history_total_cells;	// Accumulated number of cells created
	int history_dead_cells;		// Accumulated number of dead cells
	int history_max_alive_cells;	// Maximum number of cells alive in a step
	int history_max_new_cells;	// Maximum number of cells created in a step
	int history_max_dead_cells;	// Maximum number of cells died in a step
	int history_max_age;		// Maximum age achieved by a cell
	float history_max_food;		// Maximum food level in a position of the culture
} Statistics;


/* 
 * Macro function to simplify accessing with two coordinates to a flattened array
 * 	This macro-function can be changed and/or optimized by the students
 *
 */
#define accessMat( arr, exp1, exp2 )	arr[ (int)(exp1) * columns + (int)(exp2) ]

/*
 * Function: Choose a new direction of movement for a cell
 * 	This function can be changed and/or optimized by the students
 */
/*
 * Function: Mutation of the movement genes on a new cell
 * 	This function can be changed and/or optimized by the students
 */
void cell_mutation( Cell *cell ) {
	/* 1. Select which genes change:
	 	0 Left grows taking part of the Advance part
	 	1 Advance grows taking part of the Left part
	 	2 Advance grows taking part of the Right part
	 	3 Right grows taking part of the Advance part
	*/
	int mutation_type = (int)(4 * erand48( cell->random_seq ));
	/* 2. Select the amount of mutation (up to 50%) */
	float mutation_percentage = (float)(0.5 * erand48( cell->random_seq ));
	/* 3. Apply the mutation */
	float mutation_value;
	switch( mutation_type ) {
		case 0:
			mutation_value = cell->choose_mov[1] * mutation_percentage;
			cell->choose_mov[1] -= mutation_value;
			cell->choose_mov[0] += mutation_value;
			break;
		case 1:
			mutation_value = cell->choose_mov[0] * mutation_percentage;
			cell->choose_mov[0] -= mutation_value;
			cell->choose_mov[1] += mutation_value;
			break;
		case 2:
			mutation_value = cell->choose_mov[2] * mutation_percentage;
			cell->choose_mov[2] -= mutation_value;
			cell->choose_mov[1] += mutation_value;
			break;
		case 3:
			mutation_value = cell->choose_mov[1] * mutation_percentage;
			cell->choose_mov[1] -= mutation_value;
			cell->choose_mov[2] += mutation_value;
			break;
		default:
			fprintf(stderr,"Error: Imposible type of mutation\n");
			MPI_Abort( MPI_COMM_WORLD, EXIT_FAILURE );
	}
	/* 4. Correct potential precision problems */
	cell->choose_mov[2] = 1.0f - cell->choose_mov[1] - cell->choose_mov[0];
}

#ifdef DEBUG
/* 
 * Function: Print the current state of the simulation 
 */
void print_status( int iteration, int rows, int columns, float *culture, int num_cells, Cell *cells, int num_cells_alive, Statistics sim_stat ) {
	/* 
	 * You don't need to optimize this function, it is only for pretty printing and debugging purposes.
	 * It is not compiled in the production versions of the program.
	 * Thus, it is never used when measuring times in the leaderboard
	 */
	int i,j;

	printf("Iteration: %d\n", iteration );
	printf("+");
	for( j=0; j<columns; j++ ) printf("---");
	printf("+\n");
	for( i=0; i<rows; i++ ) {
		printf("|");
		for( j=0; j<columns; j++ ) {
			char symbol;
			if ( accessMat( culture, i, j ) >= 20 ) symbol = '+';
			else if ( accessMat( culture, i, j ) >= 10 ) symbol = '*';
			else if ( accessMat( culture, i, j ) >= 5 ) symbol = '.';
			else symbol = ' ';

			int t;
			int counter = 0;
			for( t=0; t<num_cells; t++ ) {
				int row = (int)(cells[t].pos_row);
				int col = (int)(cells[t].pos_col);
				if ( cells[t].alive && row == i && col == j ) {
					counter ++;
				}
			}
			if ( counter > 9 ) printf("(M)" );
			else if ( counter > 0 ) printf("(%1d)", counter );
			else printf(" %c ", symbol );
		}
		printf("|\n");
	}
	printf("+");
	for( j=0; j<columns; j++ ) printf("---");
	printf("+\n");
	printf("Num_cells_alive: %04d\nHistory( Cells: %04d, Dead: %04d, Max.alive: %04d, Max.new: %04d, Max.dead: %04d, Max.age: %04d, Max.food: %6f )\n\n", 
		num_cells_alive, 
		sim_stat.history_total_cells, 
		sim_stat.history_dead_cells, 
		sim_stat.history_max_alive_cells, 
		sim_stat.history_max_new_cells, 
		sim_stat.history_max_dead_cells, 
		sim_stat.history_max_age,
		sim_stat.history_max_food
	);
}

#endif

/*
 * Function: Print usage line in stderr
 */
void show_usage( char *program_name ) {
	fprintf(stderr,"Usage: %s ", program_name );
	fprintf(stderr,"<rows> <columns> <maxIter> <max_food> <food_density> <food_level> <short_rnd1> <short_rnd2> <short_rnd3> <num_cells>\n");
	fprintf(stderr,"\tOptional arguments for special food spot: [ <row> <col> <size_rows> <size_cols> <density> <level> ]\n");
	fprintf(stderr,"\n");
}


/*
 * MAIN PROGRAM
 */
int main(int argc, char *argv[]) {
	double tiempoTotalBucle1=0;
	double tiempoTotalBucle2=0;
	double tiempoTotalBucle3=0;
	double tiempoTotalBucle4=0;
	double tiempoTotalBucle5=0;
	double tiempoTotalBucle6=0;
	double tiempoTotalBucle7=0;
	double tiempoTotalBucle8=0;
	double tiempoTotalBucle9=0;
	double tiempoTotalBucle10=0;
	double tiempoBucle1;
	double tiempoBucle2;
	double tiempoBucle3;
	double tiempoBucle4;
	double tiempoBucle5;
	double tiempoBucle6;
	double tiempoBucle7;
	double tiempoBucle8;
	double tiempoBucle9;
	double tiempoBucle10;
	int i,j;

	// Simulation data
	int max_iter;			// Maximum number of simulation steps
	int rows, columns;		// Cultivation area sizes
	float *culture;			// Cultivation area values
	short *culture_cells;		// Ancillary structure to count the number of cells in a culture space

	float max_food;			// Maximum level of food on any position
	float food_density;		// Number of food sources introduced per step
	float food_level;		// Maximum number of food level in a new source

	bool food_spot_active = false;	// Special food spot: Active
	int food_spot_row = 0;		// Special food spot: Initial row
	int food_spot_col = 0;		// Special food spot: Initial row
	int food_spot_size_rows = 0;	// Special food spot: Rows size
	int food_spot_size_cols = 0;	// Special food spot: Cols size
	float food_spot_density = 0.0f;	// Special food spot: Food density
	float food_spot_level = 0.0f;	// Special food spot: Food level

	unsigned short init_random_seq[3];	// Status of the init random sequence
	unsigned short food_random_seq[3];	// Status of the food random sequence
	unsigned short food_spot_random_seq[3];	// Status of the special food spot random sequence

	int	num_cells;		// Number of cells currently stored in the list
	Cell	*cells;			// List to store cells information


	// Statistics
	Statistics sim_stat;	
	sim_stat.history_total_cells = 0;
	sim_stat.history_dead_cells = 0;
	sim_stat.history_max_alive_cells = 0;
	sim_stat.history_max_new_cells = 0;
	sim_stat.history_max_dead_cells = 0;
	sim_stat.history_max_age = 0;
	sim_stat.history_max_food = 0.0f;

	/* 0. Initialize MPI */
	int rank;
	MPI_Init( &argc, &argv );
	MPI_Comm_rank( MPI_COMM_WORLD, &rank );

	/* 1. Read simulation arguments */
	/* 1.1. Check minimum number of arguments */
	if (argc < 11) {
		fprintf(stderr, "-- Error: Not enough arguments when reading configuration from the command line\n\n");
		show_usage( argv[0] );
		MPI_Abort( MPI_COMM_WORLD, EXIT_FAILURE );
	}

	/* 1.2. Read culture sizes, maximum number of iterations */
	rows = atoi( argv[1] );
	columns = atoi( argv[2] );
	max_iter = atoi( argv[3] );

	/* 1.3. Food data */
	max_food = atof( argv[4] );
	food_density = atof( argv[5] );
	food_level = atof( argv[6] );

	/* 1.4. Read random sequences initializer */
	for( i=0; i<3; i++ ) {
		init_random_seq[i] = (unsigned short)atoi( argv[7+i] );
	}

	/* 1.5. Read number of cells */
	num_cells = atoi( argv[10] );

	/* 1.6. Read special food spot */
	if (argc > 11 ) {
		if ( argc < 17 ) {
			fprintf(stderr, "-- Error in number of special-food-spot arguments in the command line\n\n");
			show_usage( argv[0] );
			MPI_Abort( MPI_COMM_WORLD, EXIT_FAILURE );
		}
		else {
			food_spot_active = true;
			food_spot_row = atoi( argv[11] );
			food_spot_col = atoi( argv[12] );
			food_spot_size_rows = atoi( argv[13] );
			food_spot_size_cols = atoi( argv[14] );
			food_spot_density = atof( argv[15] );
			food_spot_level = atof( argv[16] );

			// Check non-used trailing arguments
			if ( argc > 17 ) {
				fprintf(stderr, "-- Error: too many arguments in the command line\n\n");
				show_usage( argv[0] );
				MPI_Abort( MPI_COMM_WORLD, EXIT_FAILURE );
			}
		}
	}

#ifdef DEBUG
	/* 1.7. Print arguments */
	printf("Arguments, Rows: %d, Columns: %d, max_iter: %d\n", rows, columns, max_iter);
	printf("Arguments, Max.food: %f, Food density: %f, Food level: %f\n", max_food, food_density, food_level);
	printf("Arguments, Init Random Sequence: %hu,%hu,%hu\n", init_random_seq[0], init_random_seq[1], init_random_seq[2]);
	if ( food_spot_active ) {
		printf("Arguments, Food_spot, pos(%d,%d), size(%d,%d), Density: %f, Level: %f\n",
			food_spot_row, food_spot_col, food_spot_size_rows, food_spot_size_cols, food_spot_density, food_spot_level );
	}
	printf("Initial cells: %d\n", num_cells );
#endif // DEBUG


	/* 1.8. Initialize random sequences for food dropping */
	for( i=0; i<3; i++ ) {
		food_random_seq[i] = (unsigned short)nrand48( init_random_seq );
		food_spot_random_seq[i] = (unsigned short)nrand48( init_random_seq );
	}

	/* 1.9. Initialize random sequences of cells */
	cells = (Cell *)malloc( sizeof(Cell) * (size_t)num_cells );
	if ( cells == NULL ) {
		fprintf(stderr,"-- Error allocating: %d cells\n", num_cells );
		MPI_Abort( MPI_COMM_WORLD, EXIT_FAILURE );
	}
	for( i=0; i<num_cells; i++ ) {
		// Initialize the cell ramdom sequences
		for( j=0; j<3; j++ ) 
			cells[i].random_seq[j] = (unsigned short)nrand48( init_random_seq );
	}


#ifdef DEBUG
	/* 1.10. Print random seed of the initial cells */
	/*
	printf("Initial cells random seeds: %d\n", num_cells );
	for( i=0; i<num_cells; i++ )
		printf("\tCell %d, Random seq: %hu,%hu,%hu\n", i, cells[i].random_seq[0], cells[i].random_seq[1], cells[i].random_seq[2] );
	*/
#endif // DEBUG

	/* 2. Start global timer */
	MPI_Barrier( MPI_COMM_WORLD );
	double ttotal = cp_Wtime();

/*
 *
 * START HERE: DO NOT CHANGE THE CODE ABOVE THIS POINT
 *
 */
	//A cada procesador se le asigna un numero de filas
	
	int my_begin; //principio de cada procesador en el array general
	int my_end; //final de cada procesador en el array general
	int my_size; //tamaño de matriz asignado a cada procesador 
	int my_num_cells = 0; //numero de celulas de cada procesador
	int procs; //numero de procesadores
	int posicionMyCells = 0; //posicion real de la celula pasada a entero, si: row=21.2 col=22.4 -> posicionMyCells=21*22 pero usando my_cells para calcularlo, se usa en el bucle 4
	int rankTemp = 0; //rango del procesador al que se va a enviar la celula (usado en el bucle 4.3 para saber el procesador al que se le esta enviando una celula)
	int guarda = 0;
	int posicionMyCellsEnMiCulture = 0;
	int posicionMyNewCells = 0;

	MPI_Comm_size(MPI_COMM_WORLD,&procs);	

	//creacion estructura mpi celula
	MPI_Datatype MPI_CELULAS;
	    int lengths[9] = { 1, 1, 1, 1, 3, 1, 1, 3, 1};

	    //MPI_Aint disp[9] = { 0, sizeof(float), sizeof(float)2, sizeof(float)3, sizeof(float)4, sizeof(float)7, sizeof(float)8, sizeof(int) + sizeof(float)8, sizeof(unsigned short) * 3 + sizeof(int) + sizeof(float)*8};

	    MPI_Aint disp[9] ={
	        offsetof(Cell, pos_row),
	        offsetof(Cell, pos_col),
	        offsetof(Cell, mov_row),
	        offsetof(Cell, mov_col),
	        offsetof(Cell, choose_mov),
	        offsetof(Cell, storage),
	        offsetof(Cell, age),
	        offsetof(Cell, random_seq),
	        offsetof(Cell, alive),
	    };

	    MPI_Datatype types[9] = { MPI_FLOAT, MPI_FLOAT, MPI_FLOAT, MPI_FLOAT, MPI_FLOAT, MPI_FLOAT, MPI_INT, MPI_UNSIGNED_SHORT, MPI_C_BOOL};

	    MPI_Type_create_struct(9, lengths, disp, types, &MPI_CELULAS);
	    MPI_Aint lb, ext;
	    MPI_Datatype MPI_CELL_EXT;
	    MPI_Type_get_extent(MPI_CELULAS, &lb, &ext);
	    MPI_Type_create_resized(MPI_CELULAS, lb, ext, &MPI_CELL_EXT);
    MPI_Type_commit(&MPI_CELL_EXT);

	int size = rows * columns;
	//int my_rows = rows/procs;
	//int my_columns = columns/procs;

	//Si el numero de filas es divisible entre el numero de procs a cada proc se le asignan rows/size filas + 2 bordes
	//Si no, se hace lo mismo pero a la última se le asigna, ademas, el resto (modulo)

	//Declaramos inicio y fin de cada proceso y el tamaño del culture de cada proceso

	//El proceso 0 no va a hacer nada, solo recibe y envia cosas 
	/*
	if(rank!=0){
		if (size%(procs-1)==0) {
			my_size = size/(procs-1);
			my_begin = my_size*(rank-1);
			my_end = (my_size*(rank))-1;
		} else {
			if (rank < procs-1 ) {
				my_size = size/(procs-1);
				my_begin = my_size*(rank-1);
				my_end = (my_size*(rank))-1;
			} else {
				my_size = size - ((procs-2) * (size/(procs-1)));
				my_begin = (size - 1) - my_size;
				my_end = size - 1;
			}
		}
	}
	*/

	//Si el numero de filas es divisible entre el numero de procs a cada proc se le asignan rows/size filas + 2 bordes
	//Si no, se asigna a todos uno mas menos al ultimo que se le asignan los restantes

	//Declaramos inicio y fin de cada proceso y el tamaño del culture de cada proceso

	//El proceso 0 no va a hacer nada, solo recibe y envia cosas 
		if (size%procs==0) {
			my_size = size/procs;
			my_begin = my_size*rank;
			my_end = (my_size*(rank+1))-1;
		} else {
			if (rank < procs-1 ) {
				my_size = size/procs + 1;
				my_begin = my_size*rank;
				my_end = (my_size*(rank+1))-1;
			} else {
				my_size = size - ((procs-1) * ((size/procs)+1));
				my_begin = (size - 1) - my_size;
				my_end = size - 1;
			}
		}

	MPI_Status state;
	//printf("el rango %i tiene tamaño %i", rank, my_size);
	/* 3. Initialize culture surface and initial cells */
	culture = (float *)malloc( sizeof(float) * (size_t)my_size); //culture dividido


	//culture = (float *)malloc( sizeof(float) * (size_t)rows * (size_t)columns );


	//float *cultureAux = (float *)malloc( sizeof(float) * (size_t)rows * (size_t)columns );
	//culture_cells = (short *)malloc( sizeof(short) * (size_t)rows * (size_t)columns );

	culture_cells = (short *)malloc( sizeof(short) * (size_t)my_size); //culture_cells dividido


	if ( culture == NULL || culture_cells == NULL ) {
		fprintf(stderr,"-- Error allocating culture structures for size: %d x %d \n", rows, columns );
		exit( EXIT_FAILURE );
		MPI_Abort( MPI_COMM_WORLD, EXIT_FAILURE );
	}
	tiempoBucle1 = cp_Wtime();


	for( i=0; i<my_size; i++ ){
			culture_cells[i] = 0.0f;
			culture[i] = 0.0;
	}

	MPI_Request request, request2;

	tiempoBucle1 = cp_Wtime() - tiempoBucle1;
	tiempoTotalBucle1 = 	tiempoTotalBucle1 + tiempoBucle1;

	tiempoBucle2 = cp_Wtime();

	for( i=0; i<num_cells; i++ ) {
		cells[i].alive = true;
		// Initial age: Between 1 and 20
		cells[i].age = 1 + (int)(19 * erand48( cells[i].random_seq ));
		// Initial storage: Between 10 and 20 units
		cells[i].storage = (float)(10 + 10 * erand48( cells[i].random_seq ));
		// Initial position: Anywhere in the culture arena
		cells[i].pos_row = (float)(rows * erand48( cells[i].random_seq ));
		cells[i].pos_col = (float)(columns * erand48( cells[i].random_seq ));
		// Movement direction: Unity vector in a random direction
		float angle = (float)(2 * M_PI * erand48( cells[i].random_seq ));
		cells[i].mov_row = sinf( angle ); 
		cells[i].mov_col = cosf( angle );

		if ( cells[i].pos_row >= rows ) cells[i].pos_row -= rows;
		if ( cells[i].pos_col >= columns ) cells[i].pos_col -= columns;
		// Movement genes: Probabilities of advancing or changing direction: The sum should be 1.00

		cells[i].choose_mov[0] = 0.33f;
		cells[i].choose_mov[1] = 0.34f;
		cells[i].choose_mov[2] = 0.33f;
		posicionMyCells = (int) cells[i].pos_row * columns + (int) cells[i].pos_col;
		//printf("\n");
		//printf("%i", posicionMyCells);
		
		if(posicionMyCells <= my_end && posicionMyCells >= my_begin){ //si la celula pertenece al proceso incrementa la variable de numero de celulas de cada proceso
			my_num_cells++;
			//printf("incrementamos%i", my_num_cells);
		}
	}	

	
	Cell *my_cells = (Cell *)malloc( sizeof(Cell) * (size_t)my_num_cells ); //creamos un array con las celulas de cada proceso
	//printf("primero%i", my_num_cells);

	my_num_cells = 0; //reseteamos el valor
	posicionMyCells = 0;
	//repetimos el bucle de encima para rellenar my_cells con las celulas de cells que le correspondan a cada procesador

	for( i=0; i<num_cells; i++ ) {
		posicionMyCells = (int) cells[i].pos_row * columns + (int) cells[i].pos_col;
		if(posicionMyCells <= my_end && posicionMyCells >= my_begin){
			my_cells[my_num_cells] = cells[i]; //añadimos la celula al array de ese proceso
			my_num_cells++; //incrementamos su numero de celulas
		}
	}
	if(rank==1){
		//printf("segundo%i", my_num_cells);
	}


	posicionMyCells = 0; //reseteamos valor para el resto de la ejecucion
	//free(cells); //liberamos el espacio de cells, a partir de aqui solo se usara my_cells que sera el array de celulas de cada procesador
	
	tiempoBucle2 = cp_Wtime() - tiempoBucle2;
	tiempoTotalBucle2 = 	tiempoTotalBucle2 + tiempoBucle2;

	// Statistics: Initialize total number of cells, and max. alive
	sim_stat.history_total_cells = num_cells;
	sim_stat.history_max_alive_cells = num_cells;


#ifdef DEBUG
	/* Show initial cells data */
	printf("Initial cells data: %d\n", num_cells );
	for( i=0; i<num_cells; i++ ) {
		printf("\tCell %d, Pos(%f,%f), Mov(%f,%f), Choose_mov(%f,%f,%f), Storage: %f, Age: %d\n",
				i,
				cells[i].pos_row,
				cells[i].pos_col,
				cells[i].mov_row,
				cells[i].mov_col,
				cells[i].choose_mov[0],
				cells[i].choose_mov[1],
				cells[i].choose_mov[2],
				cells[i].storage,
				cells[i].age );
	}
#endif // DEBUG

	/* 4. Simulation */
	float current_max_food = 0.0f;
	int num_cells_alive = num_cells;
	int iter;

	for( iter=0; iter<max_iter && current_max_food <= max_food && num_cells_alive > 0; iter++ ) {
		int step_new_cells = 0;
		int step_dead_cells = 0;

		tiempoBucle9 = cp_Wtime();
		/* 4.1. Spreading new food */
		// Across the whole culture
		int row = 0;
		int col = 0;
		
		int posVecReducido; //guardara la posicion del culture pequeño donde se va a añadir la comida, es decir si era la 208 sera la 8 del pequeño del procesador 20

		int num_new_sources = (int)(rows * columns * food_density); //no cambiamos, el numero de recursos no tiene que ver con la particion
		//printf("\n");

		float *foodVec = (float *)malloc( sizeof(float) * num_new_sources );
		int *posVec = (int *)malloc( sizeof(int) * num_new_sources );
		float *foodVec2 = (float *)malloc( sizeof(float) * num_new_sources );
		int *posVec2 = (int *)malloc( sizeof(int) * num_new_sources );

		//printf("\n");
		//printf("\n");
		for (i=0; i<num_new_sources; i++) {
			//printf("\n");
		  	row = (int)(rows * erand48( food_random_seq ));
			col = (int)(columns * erand48( food_random_seq ));
			foodVec[i] = (float)( food_level * erand48( food_random_seq )); //cuanta comida a esa casilla tamañoArray = num_new_sources
			posVec[i] = row*columns+col; // a que casilla real se le da la comidatamañoArray = num_new_sources
		}

		for (i=0; i<num_new_sources; i++){
			if(posVec[i]>=my_begin && posVec[i]<=my_end){ //comprobamos en que array culture pequeño se encuentra, es decir que procesador guardara este alimento
				posVecReducido = posVec[i]-my_begin; //para saber la posicion a la que corresponderia en el culture pequeño la posicion en la que se le esta dando comida
				culture[posVecReducido] = culture[posVecReducido] + foodVec[i]; //en la matriz culture pequeña se añade el alimento nuevo foodVec[i] al antiguo que tenia esa posicion culture[posVecReducido]
			}
		}



		//lo mismo que arriba pero con el spot de comida activado

		if ( food_spot_active ) {
			num_new_sources = (int)(food_spot_size_rows * food_spot_size_cols * food_spot_density);
			for (i=0; i<num_new_sources; i++) {
				int row = food_spot_row + (int)(food_spot_size_rows * erand48( food_spot_random_seq ));
				int col = food_spot_col + (int)(food_spot_size_cols * erand48( food_spot_random_seq ));
				foodVec2[i] = (float)( food_spot_level * erand48( food_spot_random_seq ));
				posVec2[i] = row*columns+col;
			}
		}
		if ( food_spot_active ) {
			for (i=0; i<num_new_sources; i++) {
				if(posVec2[i]>=my_begin && posVec2[i]<=my_end){ //comprobamos en que array culture pequeño se encuentra, es decir que procesador guardara este alimento
					posVecReducido = posVec2[i]-my_begin;
					culture[posVecReducido] = culture[posVecReducido] + foodVec2[i];
				}				
			}	
		}


		for(i=0; i<my_size; i++){
			printf("%lf    rank    %i",culture[i], rank);
			printf("\n");
		}
		//printf("hlassad");
		tiempoBucle9 = cp_Wtime() - tiempoBucle9;
		tiempoTotalBucle9 = 	tiempoTotalBucle9 + tiempoBucle9;

		/* 4.2. Prepare ancillary data structures */
		/* 4.2.1. Clear ancillary structure of the culture to account alive cells in a position after movement */

		tiempoBucle10 = cp_Wtime();

 		/*4.2.2. Allocate ancillary structure to store the food level to be shared by cells in the same culture place */
		float *food_to_share = (float *)malloc( sizeof(float) * my_num_cells );
		if ( culture == NULL || culture_cells == NULL ) {
			fprintf(stderr,"-- Error allocating culture structures for size: %d x %d \n", rows, columns );
			exit( EXIT_FAILURE );
			MPI_Abort( MPI_COMM_WORLD, EXIT_FAILURE );
		}
		tiempoBucle10 = cp_Wtime() - tiempoBucle10;
		tiempoTotalBucle10 = 	tiempoTotalBucle10 + tiempoBucle10;
		Cell temporal;
		tiempoBucle3 = cp_Wtime();
		int history_max_age =  0;
		//printf("numero de celulas %i     rank %i",my_num_cells, rank);
		Cell *my_cellsEnviadas = (Cell *)calloc((size_t)my_num_cells , sizeof(Cell) ); //almacena las celulas que va a enviar, en caso de que no vaya a enviar ninguna se envia vacion (siempre se envia)
		//printf("bbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbb%i                  ", rank);
		//printf("\n");
		int *rank_Enviadas = (int *)calloc((size_t)my_num_cells , sizeof(int)); //almacena el rank al que se van a enviar las celulas
		//printf("CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC%i                  ", rank);
		//printf("\n");
		int *numero_Recibos = (int *)calloc((size_t)procs , sizeof(int)); //vector que guarda el numero de celulas que va a recibir cada proceso numero_Recibido([0],[2],[3]...)
		//printf("DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD%i                  ", rank);
		//printf("\n");
		//El proceso 0 va a recibir 0, el proceso 1 va a recibir 2, el proceso 2 va a recibir 3 ...
		int enviados = 0; //numero de celulas que envia cada proceso

		//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
		//En caso de que nos de problemas de memoria hacer dos veces el bucle de debajo, una para saber el tamaño exacto que hay que reservar en las 2 matrices de encima y otro para rellenarlas

		//Para hacer reduction en un array hay que marcar el rango del array que se quiere reducir, como en este caso es todo el array se marca con [:tam_array]E


		//printf ("hola");
		//printf("\n");
		//printf ("%i",my_num_cells);

		// 4.3 Movimiento celulas
		//printf("hola%i",my_num_cells);
		//printf("wdqwdwqdqwdqwdqwdqwdqwdqwdqw%i", iter);

		//printf("\n"); 
		/*
		for (i=0; i<num_cells; i++) { //recorre hasta el numero de celulas totales
			posicionMyCells = (int) my_cells[i].pos_row * columns + (int) my_cells[i].pos_col;
			if(posicionMyCells <= my_end && posicionMyCells >= my_begin){ //comprueba que este en el proceso
		*/
		for (i=0; i<my_num_cells; i++){
				//int mensaje = 0; //creamos mensaje para saber si hemos recibido objeto
				//MPI_Send(&mensaje, 1, MPI_INT, MPI_ANY_SOURCE, 1, MPI_COMM_WORLD, &state);
				//MPI_Recv(&mensaje, 1, MPI_INT, MPI_ANY_SOURCE, 1, MPI_COMM_WORLD, &state); //si hemos recibido un objeto el mensaje estara a 1 si no estara a 0
				// en caso de que no funcione con recv no bloqueante, enviar y recibir aqui el mensaje
				//MPI_Send(&mensaje, 1, MPI_INT, MPI_ANY_SOURCE, 1, MPI_COMM_WORLD);

				//printf("colANTES%i",(int) my_cells[i].pos_col);
				//printf("\n");

				my_cells[i].age ++;
				// Statistics: Max age of a cell in the simulation history
				if ( my_cells[i].age > history_max_age ) history_max_age = my_cells[i].age;

				/* 4.3.1. Check if the cell has the needed energy to move or keep alive */
				if ( my_cells[i].storage < 0.1f ) {
					// Cell has died
					my_cells[i].alive = false;
					step_dead_cells ++;
					//printf("hola");
					//printf("\n");
					continue;
				}
				else if ( my_cells[i].storage < 1.0f ) {
					// Almost dying cell, it cannot move, only if enough food is dropped here it will survive
					my_cells[i].storage -= 0.2f;
				}
				else {
					// Consume energy to move
					my_cells[i].storage -= 1.0f;

					/* 4.3.2. Choose movement direction */
					float prob = (float)erand48( my_cells[i].random_seq );
					if ( prob < my_cells[i].choose_mov[0] ) {
						// Turn left (90 degrees)
						float tmp = my_cells[i].mov_col;
						my_cells[i].mov_col = cells[i].mov_row;
						my_cells[i].mov_row = -tmp;
					}
					else if ( prob >= my_cells[i].choose_mov[0] + my_cells[i].choose_mov[1] ) {
						// Turn right (90 degrees)
						float tmp = my_cells[i].mov_row;
						my_cells[i].mov_row = my_cells[i].mov_col;
						my_cells[i].mov_col = -tmp;
					}
					// else do not change the direction

					/* 4.3.3. Update position moving in the choosen direction*/
					my_cells[i].pos_row += my_cells[i].mov_row;
					my_cells[i].pos_col += my_cells[i].mov_col;
					// Periodic arena: Left/Rigth edges are connected, Top/Bottom edges are connected
					if ( my_cells[i].pos_row < 0 ) my_cells[i].pos_row += rows;
					if ( my_cells[i].pos_row >= rows ) my_cells[i].pos_row -= rows;
					if ( my_cells[i].pos_col < 0 ) my_cells[i].pos_col += columns;
					if ( my_cells[i].pos_col >= columns ) my_cells[i].pos_col -= columns;
            	}	
            	//printf("colDESPUES%i",(int) my_cells[i].pos_col);
            	//printf("\n");


	            posicionMyCells = (int) my_cells[i].pos_row * columns + (int) my_cells[i].pos_col; //posicion nueva de la celula despues del movimiento



	            if(posicionMyCells<my_begin || posicionMyCells>my_end){ //Comprobamos si la celula se ha movido de procesador
	            	rankTemp = posicionMyCells/my_size; //Calculamos el procesador al que se ha movido hay que sumar 1 por que hemos quitado el proc 0
	            	if(rankTemp>procs-2){//El ultimo es el unico que puede dar mal
	            		rankTemp = procs-1;
	            	}
	            	printf("adiosssssssssssssssssssssssssssssssss%i",iter);
	            	printf("\n");
	            	rank_Enviadas[enviados] = rankTemp; //Guardamos el rank del proceso al que se va a enviar la celula
	            	my_cellsEnviadas[enviados] = my_cells[i]; //Guardamos la celula que vamos a enviar
	            	numero_Recibos[rankTemp] ++; //Aumentamos el numero de recibos de cada proceso

	            	my_cells[i].alive = false; //Matamos la celula para que descaparezca de este proceso;
	            	enviados++;	//Aumentamos el numero de celulas enviadas

	            }

		} // End cell movements

		int contadorCellsRecibidas = 0;


		int *numero_RecibosTotales = (int *)calloc((size_t)procs , sizeof(int));

		for(i=0; i<procs; i++){
			MPI_Allreduce(&numero_Recibos[i], &numero_RecibosTotales[i], 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD); //todos los procesos van a tener el array numero_RecibosTotales con el numero de celulas que van a recibir
		}	
		
		contadorCellsRecibidas = numero_RecibosTotales[rank];
		

		Cell *celulasRecibidas = (Cell *)malloc( sizeof(int) * (size_t)contadorCellsRecibidas ); //array en el que se van a almacenar las celulas recibidas TAM = celulas que va a recibir


		//todos los que no sean el proceso 0
		if(enviados!=0){ //si han enviado alguna celula

			for(int k = 0; k<enviados; k++){ //recorren el array de celulas enviadas y vamos enviando una a una sabiendo el rango al que se le envia y la celula que se envia
				//printf("h0");
			//falta el buffer************************************************************************************************************************************************************************
				//printf("enviamos");
				//printf("edadEnviar       %i quien envia         %i    a donde envia       %i    cuando enviar    %i", my_cellsEnviadas[k].age, rank, rank_Enviadas[k], iter);
				//printf("\n");
				//printf("\n");
				MPI_Send(&my_cellsEnviadas[k], 1, MPI_CELL_EXT, rank_Enviadas[k], 1, MPI_COMM_WORLD);	
				//printf("enviamowswwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwws                    %i", iter );
				//printf("\n");
				//printf("años%i",my_cellsEnviadas[k].age);					
			}
		}

					//printf("h0");

		MPI_Barrier(MPI_COMM_WORLD); //si no no tira

		if(contadorCellsRecibidas!=0){ //si ha recibido alguna celula			
			//printf("wssssssssssssssssssssssssssssssssssssssssswsswswswswswswswswswswswswsw");
			for(int k = 0; k<contadorCellsRecibidas; k++){ //recorre hasta el total que tiene que recibir
				//printf("h0                    %i", iter );
				//printf("\n");
				Cell celulaRecibida; //celula que se va a recibir
				MPI_Recv(&celulaRecibida, 1, MPI_CELL_EXT, MPI_ANY_SOURCE, 1, MPI_COMM_WORLD, &state); //recibe una celula
				//MPI_Wait(MPI_COMM_WORLD);
				//printf("edadRecibir     %i   lo recibe       %i      cuando recibe %i", celulaRecibida.age, rank, iter);
				//printf("\n");
				//printf("años%i",celulasRecibida.age);
				celulasRecibidas[k] = celulaRecibida; //almacena la celula recibida en el array de celulas totales recibidas
			}				
		}
		//printf("h11111111111111");

		//MPI_Barrier(MPI_COMM_WORLD); //si no no tira

		int my_num_cellsTemp = my_num_cells; 

		my_num_cells += contadorCellsRecibidas;


		my_cells = (Cell *)realloc(my_cells, sizeof(Cell) * my_num_cells );




		int indice = 0;

		//Fusionar las celulas recibidas a my_cells
		for(i=0; i<my_num_cells; i++){
			if(my_num_cellsTemp<=i){
				my_cells[i] = celulasRecibidas[indice];
				indice++;
			}
		}

		

		/*
		if(rank==1){
			printf("\n");
			printf("Esto si que no me lo esperaba%i", my_cells[3].age);
		}*/
		/*
		if(rank==0){
			for(i=0; i<my_size; i++){
				printf("antes de recbirlo%i", culture_cells[i]);
				printf("\n");
			}						
		}*/

		for (i=0; i<my_num_cells; i++) { 
			/*
			if(rank==0){
				printf("Estoy entrando en el bucle");
				printf("\n");
			}*/
			posicionMyCellsEnMiCulture = (int) my_cells[i].pos_row * columns + (int) my_cells[i].pos_col;
			posicionMyCellsEnMiCulture -= my_begin;
			//*******************************************************************************************************************************************************************************************
			//es lo mismo la linea de debajo y la comentada?
			/*
			if(rank==0){
				printf("Vamos a incrementar la posicion%i",posicionMyCells);
				printf("\n");
			}*/
			culture_cells[posicionMyCellsEnMiCulture]++;

			//acessMat( culture_cells, my_cells[i].pos_row, my_cells[i].pos_col ) ++; 
			food_to_share[i] = culture_cells[posicionMyCellsEnMiCulture];
		}
		/*
		if(rank==0){
			for(i=0; i<my_size; i++){
				printf("despues de recibirlo%i", culture_cells[i]);
				printf("\n");
			}						
		}*/



		//Hacer realloc con el nuevo tamaño de my_cells -> Sera contadorCellsRecibidas + my_num_cells****************************************************************************************************
		//Falta añadir las celulasRecibidas al array de my_cells de cada proceso*************************************************************************************************************************
		//Hacer el accesMat del bucle 3 de culture_cells y el food_to_share =**************************************************************************************************************************** 

		num_cells_alive -= step_dead_cells;
		
		if(rank==0){
			//printf("Historymaxage antes    %i", history_max_age);
			//printf("\n");
		}
		
		int peque = 0;	

		MPI_Reduce(&history_max_age, &peque, 1, MPI_INT, MPI_MAX, 0, MPI_COMM_WORLD);

		/*
		free (my_cellsEnviadas);
		free (rank_Enviadas);
		free (numero_Recibos);
		free (celulasRecibidas);
		free (numero_RecibosTotales);
		*/
		enviados = 0;
		/*
		my_cellsEnviadas = (Cell *)malloc( sizeof(int) * my_num_cells );
		rank_Enviadas = (int *)malloc( sizeof(int) * my_num_cells );
		numero_Recibos = (int *)malloc( sizeof(Cell) * procs );
		numero_RecibosTotales = (int *)malloc( sizeof(int) * procs );
		*/
		//printf("llega el 0");
		//printf("\n");

		if(sim_stat.history_max_age < peque){
			sim_stat.history_max_age =  peque;
		}

		///printf("llegamos00            %i",         rank);
		//printf("\n");

		if(rank==0){
			//printf("Historymaxage despues    %i", peque);
			//printf("\n");
		}

		//MPI_reduce(&sim_stat.history_max_age, &guarda, 1, MPI_FLOAT, MPI_MAX, 0, MPI_COMM_WORLD);
		//history_max_age = guarda;

		tiempoBucle3 = cp_Wtime() - tiempoBucle3;

		tiempoTotalBucle3 =	tiempoTotalBucle3 + tiempoBucle3;

		/* 4.4. Cell actions */
		// Space for the list of new cells (maximum number of new cells is num_cells)

		//printf("%i", my_num_cells);
		Cell *new_cells = (Cell *)malloc( sizeof(Cell) * my_num_cells );


		//Cell *my_new_cells = (Cell *)malloc( sizeof(Cell) * my_num_cells ); //QUEREMOS UN ARRAY NEW_CELLS PARA CADA PROCESADOR O NO LO QUEREMOS DIVIDIR

		if ( new_cells == NULL ) {
			fprintf(stderr,"-- Error allocating new cells structures for: %d cells\n", my_num_cells );
			exit( EXIT_FAILURE );
			MPI_Abort( MPI_COMM_WORLD, EXIT_FAILURE );
		}

		Cell *my_cellsEnviadas2 = (Cell *)calloc((size_t)my_num_cells , sizeof(Cell)); //almacena las celulas que va a enviar, en caso de que no vaya a enviar ninguna se envia vacion (siempre se envia)
		int *rank_Enviadas2 = (int *)calloc((size_t)my_num_cells , sizeof(int)); //almacena el rank al que se van a enviar las celulas
		int *numero_Recibos2 = (int *)calloc((size_t)procs ,  sizeof(int) ); //vector que guarda el numero de celulas que va a recibir cada proceso numero_Recibido([0],[2],[3]...)

		tiempoBucle4 = cp_Wtime();
		for (i=0; i<my_num_cells; i++) {
			if ( my_cells[i].alive ) {
				/* 4.4.1. Food harvesting */
				float food = food_to_share[i];
				posicionMyCellsEnMiCulture = (int) my_cells[i].pos_row * columns + (int) my_cells[i].pos_col;
				posicionMyCellsEnMiCulture -= my_begin;
				short count = culture_cells[posicionMyCellsEnMiCulture];
				float my_food = food / count;
				my_cells[i].storage += my_food;

				/* 4.4.2. Split cell if the conditions are met: Enough maturity and energy */
				if ( my_cells[i].age > 30 && my_cells[i].storage > 20 ) {

					int new = 0;

					new = step_new_cells ++;

					
					// New cell is a copy of parent cell
					new_cells[ new ] = my_cells[i];

					// Split energy stored and update age in both cells
					my_cells[i].age = 1;
					my_cells[i].storage /= 2.0f;

					new_cells[ new ].storage /= 2.0f;
					
					new_cells[ new ].age = 1;

					// Random seed for the new cell, obtained using the parent random sequence
					new_cells[ new ].random_seq[0] = (unsigned short)nrand48( cells[i].random_seq );
					new_cells[ new ].random_seq[1] = (unsigned short)nrand48( cells[i].random_seq );
					new_cells[ new ].random_seq[2] = (unsigned short)nrand48( cells[i].random_seq );
					// Both cells start in random directions
					
					float angle = (float)(2 * M_PI * erand48( my_cells[i].random_seq ));
					my_cells[i].mov_row = sinf( angle );
					my_cells[i].mov_col = cosf( angle );
					

					
					angle = (float)(2 * M_PI * erand48( new_cells[new].random_seq ));
					new_cells[new].mov_row = sinf( angle );
					new_cells[new].mov_col = cosf( angle );				

					cell_mutation( &cells[i] );

					cell_mutation( &new_cells[ new ] );

					posicionMyCells = (int) my_cells[i].pos_row * columns + (int) my_cells[i].pos_col;

					posicionMyNewCells = (int) new_cells[i].pos_row * columns + (int) new_cells[i].pos_col;

					if(posicionMyCells<my_begin || posicionMyCells>my_end){ //Comprobamos si la celula se ha movido de procesador
						//printf("sssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss%i", rank);
						//printf("\n");
	            		rankTemp = posicionMyCells/my_size; //Calculamos el procesador al que se ha movido hay que sumar 1 por que hemos quitado el proc 0
	            		if(rankTemp>procs-2){//El ultimo es el unico que puede dar mal
	            			rankTemp = procs-1;
	            		}
	            		rank_Enviadas2[enviados] = rankTemp; 
	            		my_cellsEnviadas2[enviados] = my_cells[i]; 
	            		numero_Recibos2[rankTemp] ++; 

	            		my_cells[i].alive = false; 
	            		enviados++;	
	            	}

	            	if(posicionMyNewCells<my_begin || posicionMyNewCells>my_end){ //Comprobamos si la celula se ha movido de procesador
	            		//printf("sssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss%i", rank);
	            		//printf("\n");
	            		rankTemp = posicionMyNewCells/my_size; //Calculamos el procesador al que se ha movido hay que sumar 1 por que hemos quitado el proc 0
	            		if(rankTemp>procs-2){//El ultimo es el unico que puede dar mal
	            			rankTemp = procs-1;
	            		}
	            		rank_Enviadas2[enviados] = rankTemp; 
	            		my_cellsEnviadas2[enviados] = my_cells[i]; 
	            		numero_Recibos2[rankTemp] ++; 

	            		my_cells[i].alive = false; 
	            		enviados++;	
	            	}
				}
			}
		}  // End cell actions
		contadorCellsRecibidas = 0;
		int *numero_RecibosTotales2 = (int *)calloc( (size_t)procs ,  sizeof(int) );

		for(i=0; i<procs; i++){
			MPI_Allreduce(&numero_Recibos2[i], &numero_RecibosTotales2[i], 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD); //todos los procesos van a tener el array numero_RecibosTotales con el numero de celulas que van a recibir
		}	
		
		contadorCellsRecibidas = numero_RecibosTotales2[rank];
		//printf("esta funcionando%i          rank            %i", contadorCellsRecibidas, rank);
		//printf("\n");

		Cell *celulasRecibidas2 = (Cell *)malloc( sizeof(int) * (size_t)contadorCellsRecibidas );
		
		if(enviados!=0){ //si han enviado alguna celula			
			//printf("\n");
			for(int k = 0; k<enviados; k++){ //recorren el array de celulas enviadas y vamos enviando una a una sabiendo el rango al que se le envia y la celula que se envia
			//falta el buffer************************************************************************************************************************************************************************
				//printf("enviamos");
				//printf("edadEnviar%i quien envia         %i    a donde envia       %i", my_cellsEnviadas[k].age, rank, rank_Enviadas[k]);
				//printf("\n");
				//printf("\n");
				MPI_Send(&my_cellsEnviadas2[k], 1, MPI_CELL_EXT, rank_Enviadas2[k], 1, MPI_COMM_WORLD);	
				//printf("años%i",my_cellsEnviadas[k].age);					
			}
		}

		//printf("holaaaaaaa%i", rank);
		//printf("\n");
		//printf("holaaaaaaadwddddddddddddd");
		//printf("\n");


		MPI_Barrier(MPI_COMM_WORLD); //si no no tira

		if(contadorCellsRecibidas!=0){ //si ha recibido alguna celula			
			//printf("\n");
			for(int k = 0; k<contadorCellsRecibidas; k++){ //recorre hasta el total que tiene que recibir
				Cell celulaRecibida; //celula que se va a recibir
				MPI_Recv(&celulaRecibida, 1, MPI_CELL_EXT, MPI_ANY_SOURCE, 1, MPI_COMM_WORLD, &state); //recibe una celula
				//MPI_Wait(MPI_COMM_WORLD);
				//printf("edadRecibir%i   lo recibe       %i", celulaRecibida.age, rank);
				//printf("\n");
				//printf("años%i",celulasRecibida.age);
				celulasRecibidas2[k] = celulaRecibida; //almacena la celula recibida en el array de celulas totales recibidas
			}				
		}

		my_num_cellsTemp = my_num_cells; 


		my_num_cells += contadorCellsRecibidas;

		my_cells = (Cell *)realloc( my_cells, sizeof(Cell) * my_num_cells );


		indice = 0;

		//Fusionar las celulas recibidas a my_cells
		for(i=0; i<my_num_cells; i++){
			if(my_num_cellsTemp<=i){
				my_cells[i] = celulasRecibidas2[indice];
				indice++;
			}
		}


		sim_stat.history_total_cells = sim_stat.history_total_cells + step_new_cells;
		num_cells_alive = num_cells_alive + step_new_cells;

		
		tiempoBucle4 = cp_Wtime() - tiempoBucle4;

		tiempoTotalBucle4 =	tiempoTotalBucle4 + tiempoBucle4;

		/* 4.5. Clean ancillary data structures */
		/* 4.5.1. Clean the food consumed by the cells in the culture data structure */
		tiempoBucle5 = cp_Wtime();
		int free_position = 0;
		int alive_in_main_list = 0;
		int celulasMuertasEnUnProcesador = 0;

		for (i=0; i<my_num_cells; i++) {
				if ( my_cells[i].alive ) {
					//accessMat( culture, cells[i].pos_row, cells[i].pos_col ) = 0.0f;
					posicionMyCellsEnMiCulture = (int) my_cells[i].pos_row * columns + (int) my_cells[i].pos_col;
					posicionMyCellsEnMiCulture -= my_begin;
					culture[posicionMyCellsEnMiCulture] = 0.0f; //actualizamos la casilla de solo ese procesador
					alive_in_main_list ++;
					if ( free_position != i ) {
						my_cells[free_position] = my_cells[i]; //Recolocamos vector de celulas de cada procesador
					}
					free_position ++;
				}
			//}
		}
		/*
		if(rank == 0){
				free_position = free_position_suma;
				alive_in_main_list = alive_in_main_list_suma;
		}
		*/
		//printf("\n");	
		//printf("hola buenos dias");
		tiempoBucle5 = cp_Wtime() - tiempoBucle5;
		tiempoTotalBucle5 =	tiempoTotalBucle5 + tiempoBucle5;
		/* 4.5.2. Free the ancillary data structure to store the food to be shared */
		//printf("DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD");

		tiempoBucle6 = cp_Wtime() - tiempoBucle6;

		tiempoTotalBucle6 =	tiempoTotalBucle6 + tiempoBucle6;
		// 4.6.2. Reduce the storage space of the list to the current number of cells
		my_num_cells = alive_in_main_list;
		//printf("\n");
		//printf("\n");

		my_cells = (Cell *)realloc( my_cells, sizeof(Cell) * my_num_cells ); //reasignamos espacio en funcion de las celulas vivas para cada array cells de cada procesador
		/* 4.7. Join cell lists: Old and new cells list */
		tiempoBucle7 = cp_Wtime();
		/*
		if ( step_new_cells > 0 ) {
			my_cells = (Cell *)realloc( cells, sizeof(Cell) * ( my_num_cells + step_new_cells ) ); //reasignamos espacio en funcion de las celulas vivas que habia y las nuevas que han nacido (step_new_cells)
			for (j=0; j<step_new_cells; j++)
				my_cells[ my_num_cells + j ] = new_cells[ j ]; //le asignamos la nuevas celulas almacenadas en el array new_cells
			my_num_cells += step_new_cells; //sumamos el numero de celulas de cada proceso a las que han nacido en ese proceso
		}*/
		//printf("\n");
		//printf("hola buenos 2dwddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddd");
		//printf("\n");
		tiempoBucle7 = cp_Wtime() - tiempoBucle7;

		tiempoTotalBucle7 =	tiempoTotalBucle7 + tiempoBucle7;

		/* 4.8. Decrease non-harvested food */
		tiempoBucle8 = cp_Wtime();
		current_max_food = 0.0f;

		//MPI_Scatter (culture, (rows*columns)/procs, MPI_FLOAT, cultureAux, (rows*columns)/procs, MPI_FLOAT, 0, MPI_COMM_WORLD);

		/*for( i=0; i<(rows * columns)/procs; i++ ){
			cultureAux[i] *= 0.95f;
			if ( cultureAux[i] > current_max_food )
				current_max_food = cultureAux[i];	
		}	*/


		//FALTA
		float peque2 = 0;
		//printf("hola buenos diasdwdwdwqdqdqwdqwdqwdqwd");
		for( i=0; i<my_size; i++ ){
			culture[i] *= 0.95f;
			culture_cells[i] = 0.0f;
			if ( culture[i] > current_max_food )
				current_max_food = culture[i];
		}
		MPI_Reduce(&current_max_food, &peque2, 1, MPI_FLOAT, MPI_MAX, 0, MPI_COMM_WORLD);
		current_max_food = peque2;

		//MPI_reduce(&current_max_food, &guarda, 1, MPI_FLOAT, MPI_MAX, 0, MPI_COMM_WORLD);
		//MPI_Allgather (cultureAux, (rows*columns)/procs, MPI_FLOAT, culture, (rows*columns)/procs, MPI_FLOAT, MPI_COMM_WORLD);

		//printf("guarda %lf", guarda);


		//free( new_cells );
		//free( food_to_share );



		/*

		if(rank==0)
		for( i=1; i<size+2; i++ )
			for( j=0; j<columns; j++ ) {
				accessMat( culture_cells, i, j ) *= 0.0f;
				accessMat( culture, i, j ) *= 0.95f; // Reduce 5%
				
			}

		else if(rank==procs-1)
		for( i=0; i<size+1; i++ )
			for( j=0; j<columns; j++ ) {
				accessMat( culture_cells, i, j ) *= 0.0f;
				accessMat( culture, i, j ) *= 0.95f; // Reduce 5%
				
			}

		else
		for( i=0; i<size+1; i++ )
			for( j=0; j<columns; j++ ) {
				accessMat( culture_cells, i, j ) *= 0.0f;
				accessMat( culture, i, j ) *= 0.95f; // Reduce 5%
				if ( accessMat( culture, i, j ) > current_max_food ) 
					current_max_food = accessMat( culture, i, j );
			}

		if(rank!=procs-1){
			MPI_Send(&accessMat( culture, size, 0 ), columns, MPI_FLOAT, rank+1, 1, MPI_COMM_WORLD);
			MPI_Recv(&accessMat( culture, size+1, 0 ), columns, MPI_FLOAT, rank+1, 1, MPI_COMM_WORLD,&state);
		}
		if(rank!=0){
			MPI_Send(&accessMat( culture, 1, 0 ), columns, MPI_FLOAT, rank-1, 1, MPI_COMM_WORLD);
			MPI_Recv(&accessMat( culture, 0, 0 ), columns, MPI_FLOAT, rank-1, 1, MPI_COMM_WORLD,&state);
		}
		
		*/
		tiempoBucle8 = cp_Wtime() - tiempoBucle8;

		tiempoTotalBucle8 =	tiempoTotalBucle8 + tiempoBucle8;

		/* 4.9. Statistics */
		// Statistics: Max food
		if ( current_max_food > sim_stat.history_max_food ) sim_stat.history_max_food = current_max_food;
		// Statistics: Max new cells per step
		if ( step_new_cells > sim_stat.history_max_new_cells ) sim_stat.history_max_new_cells = step_new_cells;
		// Statistics: Accumulated dead and Max dead cells per step
		sim_stat.history_dead_cells += step_dead_cells;
		if ( step_dead_cells > sim_stat.history_max_dead_cells ) sim_stat.history_max_dead_cells = step_dead_cells;
		// Statistics: Max alive cells per step
		if ( num_cells_alive > sim_stat.history_max_alive_cells ) sim_stat.history_max_alive_cells = num_cells_alive;

		MPI_Barrier(MPI_COMM_WORLD);



#ifdef DEBUG
		/* 4.10. DEBUG: Print the current state of the simulation at the end of each iteration */
		print_status( iter, rows, columns, culture, num_cells, cells, num_cells_alive, sim_stat );
#endif // DEBUG
	}


/*
 *
 * STOP HERE: DO NOT CHANGE THE CODE BELOW THIS POINT
 *
 */
	/* 5. Stop global time */
	MPI_Barrier( MPI_COMM_WORLD );
	ttotal = cp_Wtime() - ttotal;

#ifdef DEBUG
	printf("List of cells at the end of the simulation: %d\n\n", num_cells );
	for( i=0; i<num_cells; i++ ) {
		printf("Cell %d, Alive: %d, Pos(%f,%f), Mov(%f,%f), Choose_mov(%f,%f,%f), Storage: %f, Age: %d\n",
				i,
				cells[i].alive,
				cells[i].pos_row, 
				cells[i].pos_col, 
				cells[i].mov_row, 
				cells[i].mov_col, 
				cells[i].choose_mov[0], 
				cells[i].choose_mov[1], 
				cells[i].choose_mov[2], 
				cells[i].storage,
				cells[i].age );
	}
#endif // DEBUG

	/* 6. Output for leaderboard */
	if ( rank == 0 ) {
		printf("\n");
		/* 6.1. Total computation time */
		printf("Time: %lf\n", ttotal );

		printf("Inicialización 1 3: %lf", tiempoTotalBucle1);
		printf("\n");
		printf("Inicialización 2 3: %lf", tiempoTotalBucle2);
		printf("\n");
		printf("Tirar nueva comida 4.1: %lf", tiempoTotalBucle9);
		printf("\n");
		printf("Preparar estructuras 4.2: %lf", tiempoTotalBucle10);
		printf("\n");
		printf("Movimiento celulas 4.3: %lf", tiempoTotalBucle3);
		printf("\n");
		printf("Accion celulas 4.4: %lf", tiempoTotalBucle4);
		printf("\n");
		printf("Limpieza comida consumida 4.5: %lf", tiempoTotalBucle5);
		printf("\n");
		printf("Limpieza celulas muertas 4.6: %lf", tiempoTotalBucle6);
		printf("\n");
		printf("Nuevas celulas 4.7: %lf", tiempoTotalBucle7);
		printf("\n");
		printf("Decremento comida que no se consume 4.8: %lf", tiempoTotalBucle8);
		printf("\n");

		/* 6.2. Results: Number of iterations and other statistics */
		printf("Result: %d, ", iter);
		printf("%d, %d, %d, %d, %d, %d, %d, %f\n", 
			num_cells_alive, 
			sim_stat.history_total_cells, 
			sim_stat.history_dead_cells, 
			sim_stat.history_max_alive_cells, 
			sim_stat.history_max_new_cells, 
			sim_stat.history_max_dead_cells, 
			sim_stat.history_max_age,
			sim_stat.history_max_food
		);
	}

	/* 7. Free resources */	
	free( culture );
	free( culture_cells );
	free( cells );

	/* 8. End */
	MPI_Finalize();
	return 0;
}
