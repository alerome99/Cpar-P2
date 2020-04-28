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
void cell_new_direction( Cell *cell ) {
	float angle = (float)(2 * M_PI * erand48( cell->random_seq ));
	cell->mov_row = sinf( angle );
	cell->mov_col = cosf( angle );
}

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

    int resto = size%procs;

    my_size = size/procs ;

    if(rank < resto ){
        my_size=my_size+1;
        my_begin = my_size * rank ;
        my_end = my_begin+(my_size-1);
    }else{
        my_begin = resto + my_size * rank;
        my_end = my_begin+(my_size-1) ;
    }



	MPI_Status state;
	culture = (float *)malloc( sizeof(float) * (size_t)my_size); //culture dividido

	culture_cells = (short *)malloc( sizeof(short) * (size_t)my_size); //culture_cells dividido


	if ( culture == NULL || culture_cells == NULL ) {
		fprintf(stderr,"-- Error allocating culture structures for size: %d x %d \n", rows, columns );
		exit( EXIT_FAILURE );
		MPI_Abort( MPI_COMM_WORLD, EXIT_FAILURE );
	}

	for( i=0; i<my_size; i++ ){
			culture_cells[i] = 0.0f;
			culture[i] = 0.0;
	}

	for( i=0; i<num_cells; i++ ) {
		cells[i].alive = true;
		// Initial age: Between 1 and 20
		cells[i].age = 1 + (int)(19 * erand48( cells[i].random_seq ));
		// Initial storage: Between 10 and 20 units
		cells[i].storage = (float)(10 + 10 * erand48( cells[i].random_seq ));
		//printf("storage   %lf    ",cells[i].storage);
		//printf("\n");
		// Initial position: Anywhere in the culture arena
		cells[i].pos_row = (float)(rows * erand48( cells[i].random_seq ));
		cells[i].pos_col = (float)(columns * erand48( cells[i].random_seq ));
		// Movement direction: Unity vector in a random direction
		posicionMyCells = (int) cells[i].pos_row * columns + (int) cells[i].pos_col;
		
		if(posicionMyCells <= my_end && posicionMyCells >= my_begin){ //si la celula pertenece al proceso incrementa la variable de numero de celulas de cada proceso
			float angle = (float)(2 * M_PI * erand48( cells[i].random_seq ));
			cells[i].mov_row = sinf( angle ); 
			cells[i].mov_col = cosf( angle );

			if ( cells[i].pos_row >= rows ) cells[i].pos_row -= rows;
			if ( cells[i].pos_col >= columns ) cells[i].pos_col -= columns;
			// Movement genes: Probabilities of advancing or changing direction: The sum should be 1.00

			cells[i].choose_mov[0] = 0.33f;
			cells[i].choose_mov[1] = 0.34f;
			cells[i].choose_mov[2] = 0.33f;
		}
		else{
			cells[i].alive = false;
		}
	}

	int free_position = 0;
	int alive_in_main_list = 0;

	for( i=0; i<num_cells; i++ ) {
		if ( cells[i].alive ) {
			alive_in_main_list ++;
			if ( free_position != i ) {
				cells[free_position] = cells[i];
			}
			free_position ++;
		}
	}

	my_num_cells = free_position ;

	cells = (Cell *)realloc(cells, sizeof(Cell) * my_num_cells );	

	posicionMyCells = 0; //reseteamos valor para el resto de la ejecucion
	 //liberamos el espacio de cells, a partir de aqui solo se usara my_cells que sera el array de celulas de cada procesador
	

	// Statistics: Initialize total number of cells, and max. alive
	sim_stat.history_total_cells = my_num_cells;
	sim_stat.history_max_alive_cells = num_cells;
	//float *food_to_share = (float *)calloc((size_t)my_num_cells , sizeof(float));



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
	//posicionMyNewCells=0;
	float current_max_food = 0.0f;
	int num_cells_alive = num_cells;
	int iter;
	int posVecReducido;
	int enviados; 
	int step_new_cells;
	int row ;
	int col;
	int posicion;
	int step_dead_cells; 
	int num_new_sources=0;
	float prob;
	float tmp;
	int enviadosTotal = 0;
    int offset = 0;
    int my_num_cellsTemp; 
    int indiceTemp;
    int indice;
    float reduceMaxFood;
	int reduceCellsAlive;
	int reduceCellsMuertas;
	int reduceStepNewCells;
	int cociente;
    int aux;
    float food;
    float my_food;
    short count;
	for( iter=0; iter<max_iter && current_max_food <= max_food && num_cells_alive > 0; iter++ ) {

		step_new_cells = 0;
		step_dead_cells = 0;

		/* 4.1. Spreading new food */
		// Across the whole culture
		row = 0;
		col = 0;
		posicion = 0;
		
		posVecReducido; //guardara la posicion del culture pequeño donde se va a añadir la comida, es decir si era la 208 sera la 8 del pequeño del procesador 20

		num_new_sources = (int)(rows * columns * food_density); //no cambiamos, el numero de recursos no tiene que ver con la particion

		float *foodVec = (float *)malloc( sizeof(float) * num_new_sources );
		int *posVec = (int *)malloc( sizeof(int) * num_new_sources );
		float *foodVec2 = (float *)malloc( sizeof(float) * num_new_sources );
		int *posVec2 = (int *)malloc( sizeof(int) * num_new_sources );


		for (i=0; i<num_new_sources; i++) {
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


 		/* 4.2.2. Allocate ancillary structure to store the food level to be shared by cells in the same culture place */
		if ( culture == NULL || culture_cells == NULL ) {
			fprintf(stderr,"-- Error allocating culture structures for size: %d x %d \n", rows, columns );
			MPI_Abort( MPI_COMM_WORLD, EXIT_FAILURE );
		}

		Cell *my_cellsEnviadas = (Cell *)calloc((size_t)my_num_cells , sizeof(Cell) ); //almacena las celulas que va a enviar, en caso de que no vaya a enviar ninguna se envia vacion (siempre se envia)
		int *numero_Envios = (int *)calloc((size_t)procs , sizeof(int)); //vector que guarda el numero de celulas que va a recibir cada proceso numero_Recibido([0],[2],[3]...)

		enviados = 0; //numero de celulas que envia cada proceso
		/* 4.3. Cell movements */

		for (i=0; i<my_num_cells; i++) {
				cells[i].age ++;
				// Statistics: Max age of a cell in the simulation history
				if ( cells[i].age > sim_stat.history_max_age ) sim_stat.history_max_age = cells[i].age;

				/* 4.3.1. Check if the cell has the needed energy to move or keep alive */
				if ( cells[i].storage < 0.1f ) {
					cells[i].alive = false;


					num_cells_alive --;
					step_dead_cells ++;
					continue;
				}
				if ( cells[i].storage < 1.0f ) {

					// Almost dying cell, it cannot move, only if enough food is dropped here it will survive
					cells[i].storage -= 0.2f;
				}
				else {
					
					// Consume energy to move
					cells[i].storage -= 1.0f;
					//printf("%lf",cells[i].storage);
					//printf("\n");


					/* 4.3.2. Choose movement direction */
					prob = (float)erand48( cells[i].random_seq );
					if ( prob < cells[i].choose_mov[0] ) {
						// Turn left (90 degrees)
						tmp = cells[i].mov_col;
						cells[i].mov_col = cells[i].mov_row;
						cells[i].mov_row = -tmp;
					}

					else if ( prob >= cells[i].choose_mov[0] + cells[i].choose_mov[1] ) {
						// Turn right (90 degrees)
						tmp = cells[i].mov_row;
						cells[i].mov_row = cells[i].mov_col;
						cells[i].mov_col = -tmp;
					}

					// else do not change the direction
					
					/* 4.3.3. Update position moving in the choosen direction*/
					cells[i].pos_row += cells[i].mov_row;
					cells[i].pos_col += cells[i].mov_col;
					// Periodic arena: Left/Rigth edges are connected, Top/Bottom edges are connected
					if ( cells[i].pos_row < 0 ) cells[i].pos_row += rows;
					if ( cells[i].pos_row >= rows ) cells[i].pos_row -= rows;
					if ( cells[i].pos_col < 0 ) cells[i].pos_col += columns;
					if ( cells[i].pos_col >= columns ) cells[i].pos_col -= columns;
				}

				posicionMyCells = (int) cells[i].pos_row * columns + (int) cells[i].pos_col; //posicion nueva de la celula despues del movimiento

	            if(posicionMyCells<my_begin || posicionMyCells>my_end){ //Comprobamos si la celula se ha movido de procesador
	            	cociente = size/procs;
           			rankTemp = posicionMyCells/cociente;
           			if(rankTemp<resto){
           				aux = rankTemp;
           			}else{
           				aux = resto;
           			}
           			if((posicionMyCells - rankTemp*cociente< aux) && (rankTemp != 0)){
            			rankTemp --;
           			}
	            	my_cellsEnviadas[enviados] = cells[i]; //Guardamos la celula que vamos a enviar
	            	numero_Envios[rankTemp] ++; //Aumentamos el numero de recibos de cada proceso
	            	cells[i].alive = false; //Matamos la celula para que descaparezca de este proceso;
	            	enviados++;	//Aumentamos el numero de celulas enviadas
	            }
		} // End cell movements

		free_position = 0;
		alive_in_main_list = 0;

		for( i=0; i<my_num_cells; i++ ) {
			if ( cells[i].alive ) {
				alive_in_main_list ++;
				if ( free_position != i ) {
					cells[free_position] = cells[i];
				}
				free_position ++;
			}
		}

		my_num_cells = free_position ;

		cells = (Cell *)realloc(cells, sizeof(Cell) * my_num_cells );

		int *numero_Recibos = (int *)malloc(sizeof(int) * (size_t)procs );

		MPI_Alltoall(numero_Envios, 1, MPI_INT, numero_Recibos, 1, MPI_INT, MPI_COMM_WORLD);

		int *dispEnvio = (int *)calloc(sizeof(int), (size_t)procs);

	    int *dispRecibo = (int *)calloc(sizeof(int), (size_t)procs);

	    int *aux2 = (int *)calloc(sizeof(int),(size_t)procs);

        int cellsRecibidasPorProceso = 0;


        for(i=1; i<procs; i++){

            dispEnvio[i] = dispEnvio[i-1] + numero_Envios[i-1];

            dispRecibo[i] = dispRecibo[i-1] + numero_Recibos[i-1];

            cellsRecibidasPorProceso += numero_Recibos[i];
        }

        cellsRecibidasPorProceso += numero_Recibos[0];
        

        Cell *numero_CelulasRecibidasCadaUno = (Cell *)malloc(sizeof(Cell ) * (size_t)cellsRecibidasPorProceso);

        Cell *celulasEnviadas = (Cell *)malloc(sizeof(Cell ) * (size_t)enviados);
		
        for(i=0; i<enviados; i++){
        	posicionMyCells = (int) my_cellsEnviadas[i].pos_row * columns + (int) my_cellsEnviadas[i].pos_col; //posicion nueva de la celula despues del movimiento

			cociente = size/procs;
           	rankTemp = posicionMyCells/cociente;
           			if(rankTemp<resto){
           				aux = rankTemp;
           			}else{
           				aux = resto;
           			}
           			if((posicionMyCells - rankTemp*cociente< aux) && (rankTemp != 0)){
            			rankTemp --;
           			}

            celulasEnviadas[dispEnvio[rankTemp] + aux2[rankTemp]++] = my_cellsEnviadas[i];
        }


        MPI_Alltoallv(celulasEnviadas, numero_Envios, dispEnvio, MPI_CELL_EXT, numero_CelulasRecibidasCadaUno,numero_Recibos
			, dispRecibo, MPI_CELL_EXT, MPI_COMM_WORLD);

        my_num_cellsTemp = my_num_cells; 

        my_num_cells += cellsRecibidasPorProceso;

        cells = (Cell *)realloc(cells, sizeof(Cell) * my_num_cells );

		indice = 0;

		for (i=my_num_cellsTemp; i<my_num_cells; i++){
			cells[i] = numero_CelulasRecibidasCadaUno[indice];
			indice ++;
		}

		
		for (i=0; i<my_num_cells; i++) { 		
			posicionMyCellsEnMiCulture = (int) cells[i].pos_row * columns + (int) cells[i].pos_col;

			posicionMyCellsEnMiCulture -= my_begin;	

			culture_cells[posicionMyCellsEnMiCulture]++;			
		}


		/* 4.4. Cell actions */
		// Space for the list of new cells (maximum number of new cells is num_cells)
		Cell *new_cells = (Cell *)malloc( sizeof(Cell) * my_num_cells );
		if ( new_cells == NULL ) {
			fprintf(stderr,"-- Error allocating new cells structures for: %d cells\n", num_cells );
			MPI_Abort( MPI_COMM_WORLD, EXIT_FAILURE );
		}

		posicionMyCellsEnMiCulture = 0;

		//int *medidorDeStepNewCells = (int *)calloc((size_t)my_num_cells , sizeof(int));

		for (i=0; i<my_num_cells; i++) {
				/* 4.4.1. Food harvesting */
				posicionMyCellsEnMiCulture = (int) cells[i].pos_row * columns + (int) cells[i].pos_col;
				posicionMyCellsEnMiCulture -= my_begin;
				food = culture[posicionMyCellsEnMiCulture];
				count = culture_cells[posicionMyCellsEnMiCulture];
				my_food = 0;
				if(count!=0){
					my_food = food / count;
				}else{
					my_food = food;
				}

				cells[i].storage += my_food;
				
				/* 4.4.2. Split cell if the conditions are met: Enough maturity and energy */
				if ( cells[i].age > 30 && cells[i].storage > 20 ) {


					num_cells_alive ++; //hace un reduce MPI_SUM ****************************************************************************************************************************************
					sim_stat.history_total_cells ++; //hace un reduce MPI_SUM ****************************************************************************************************************************************
	
					step_new_cells ++;

					// New cell is a copy of parent cell
					new_cells[ step_new_cells-1 ] = cells[i];

					// Split energy stored and update age in both cells
					cells[i].storage /= 2.0f;
					new_cells[ step_new_cells-1 ].storage /= 2.0f;
					cells[i].age = 1;
					new_cells[ step_new_cells-1 ].age = 1;

					// Random seed for the new cell, obtained using the parent random sequence
					new_cells[ step_new_cells-1 ].random_seq[0] = (unsigned short)nrand48( cells[i].random_seq );
					new_cells[ step_new_cells-1 ].random_seq[1] = (unsigned short)nrand48( cells[i].random_seq );
					new_cells[ step_new_cells-1 ].random_seq[2] = (unsigned short)nrand48( cells[i].random_seq );

					// Both cells start in random directions
					float angle = (float)(2 * M_PI * erand48( cells[i].random_seq ));
					cells[i].mov_row = sinf( angle );
					cells[i].mov_col = cosf( angle );
					

					
					
					angle = (float)(2 * M_PI * erand48( new_cells[step_new_cells-1].random_seq ));
					new_cells[step_new_cells-1].mov_row = sinf( angle );
					new_cells[step_new_cells-1].mov_col = cosf( angle );
					
				
					// Mutations of the movement genes in both cells
					cell_mutation( &cells[i] );
					cell_mutation( &new_cells[ step_new_cells-1 ] );
				}				
		} // End cell actions	
		
		


		/* 4.5. Clean ancillary data structures */
		/* 4.5.1. Clean the food consumed by the cells in the culture data structure */

		for (i=0; i<my_num_cells; i++) {
			if ( cells[i].alive ) {
				posicionMyCellsEnMiCulture = (int) cells[i].pos_row * columns + (int) cells[i].pos_col;
				posicionMyCellsEnMiCulture -= my_begin;
				culture[posicionMyCellsEnMiCulture] = 0.0f; //actualizamos la casilla de solo ese procesador
			}
		}



		if ( step_new_cells > 0 ) {
			cells = (Cell *)realloc( cells, sizeof(Cell) * ( my_num_cells + step_new_cells ) );

			for (j=0; j<step_new_cells; j++){

				cells[ my_num_cells + j ] = new_cells[ j ];
				}
			my_num_cells += step_new_cells;
		}


		
		free( new_cells );
		free( my_cellsEnviadas );
		free( numero_Envios );
		free( numero_Recibos );
		free( numero_CelulasRecibidasCadaUno );

		/* 4.8. Decrease non-harvested food */
		current_max_food = 0.0f;
		for( i=0; i<my_size; i++ ){
			culture[i] *= 0.95f;
			culture_cells[i] = 0.0f;
			if ( culture[i] > current_max_food )
				current_max_food = culture[i];
		}



		reduceMaxFood = 0;
		reduceCellsAlive = 0;
		reduceCellsMuertas = 0;
		reduceStepNewCells = 0;



		MPI_Allreduce(&current_max_food, &reduceMaxFood, 1, MPI_FLOAT, MPI_MAX, MPI_COMM_WORLD);
		current_max_food = reduceMaxFood;


		MPI_Allreduce(&my_num_cells, &reduceCellsAlive, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
		num_cells_alive = reduceCellsAlive;


		MPI_Reduce(&step_dead_cells, &reduceCellsMuertas, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
		step_dead_cells = reduceCellsMuertas;


		MPI_Reduce(&step_new_cells, &reduceStepNewCells, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
		step_new_cells = reduceStepNewCells;

		
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



#ifdef DEBUG
		/* 4.10. DEBUG: Print the current state of the simulation at the end of each iteration */
		print_status( iter, rows, columns, culture, num_cells, cells, num_cells_alive, sim_stat );
#endif // DEBUG
	}

	int reduceHistoryTotalCells = 0;
	MPI_Reduce(&sim_stat.history_total_cells, &reduceHistoryTotalCells, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
	if(rank==0){
		sim_stat.history_total_cells = reduceHistoryTotalCells;
	}

	int reduceMaxAge = 0;
	MPI_Reduce(&sim_stat.history_max_age, &reduceMaxAge, 1, MPI_INT, MPI_MAX, 0, MPI_COMM_WORLD);
	if(rank==0){
		sim_stat.history_max_age = reduceMaxAge;
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
