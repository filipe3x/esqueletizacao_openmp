#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include <papi.h>

#include "skeletonize.h"

#define ERROR_RETURN(retval) { fprintf(stderr, "Error PAPI %d %s:line %d: \n", retval,__FILE__,__LINE__);  exit(retval); }

typedef struct struct_name_events {
	char* name;
	int valor;
}*Name_Events, Valor_Name_Events;

static int create_events(Name_Events ***name_events, int *length_events, int **length,
		int **EventSet, long long *** values) ;

Name_Events **name_events;
int length_events;
int *length;
int *EventSet;// = PAPI_NULL;
long long **values;

void papi_print(void) {

	for (int i = 0; i < length_events; i++) {
		for (int j = 0; j < length[i]; j++) {
			fprintf(stderr,"%s:%lld\n", name_events[i][j]->name, values[i][j]);
		}
	}
}

int papi_init (void) {
  int retval;
	
  //init PAPI
  if ((retval = PAPI_library_init(PAPI_VER_CURRENT)) != PAPI_VER_CURRENT) {
    printf("Library initialization error! \n");
    return (0);
  }

  // create PAPI events
  return create_events(&name_events, &length_events, &length, &EventSet, &values);
}

void papi_finalize(void) {
  PAPI_shutdown();
}

void papi_start_event (int i) {
  int retval;

  EventSet[i] = PAPI_NULL;
  if ((retval = PAPI_create_eventset(&(EventSet[i]))) != PAPI_OK)
	ERROR_RETURN(retval);

  for (int j = 0; j < length[i]; j++) {
	if ((retval = PAPI_add_event(EventSet[i],
		  name_events[i][j]->valor)) != PAPI_OK)
		ERROR_RETURN(retval);
	}

  if ((retval = PAPI_start(EventSet[i])) != PAPI_OK)
	ERROR_RETURN(retval);
}

void papi_stop_event (int i) {
  int retval;
  
  if ((retval = PAPI_stop(EventSet[i], values[i])) != PAPI_OK)
	ERROR_RETURN(retval);
}


// creates the events one wants to measure
// returns the number of such events, or better, the nbr of times 
// the code being measured has to be repeated to measure all events
static int create_events(Name_Events ***name_events, int *length_events, int **length,
		int **EventSet, long long *** values) {

	int *aux;
	Name_Events **aux_name_events;

	*length_events = 1; //number of measurements???
	//*length_events = 4;

	*length = (int*) malloc(sizeof(int) * *length_events);
	*EventSet = (int*) malloc(sizeof(int) * *length_events);

	(*length)[0] = 1; //number of counters - arch dependent
	//(*length)[1] = 1;
	//(*length)[2] = 1;
	//(*length)[3] = 1;

	//Allocate arrays
	aux_name_events = (Name_Events**) malloc(sizeof(Name_Events*)
			* *length_events);
	*values = (long long **) malloc(sizeof(long long *) * *length_events);
	for (int i = 0; i < (*length_events); i++) {
		(*values)[i] = (long long *) malloc(sizeof(long long) * (*length)[i]);
		aux_name_events[i] = (Name_Events*) malloc(sizeof(Name_Events)
				* (*length)[i]);
		for (int j = 0; j < (*length)[i]; j++) {
			aux_name_events[i][j] = (Name_Events) malloc(
					sizeof(Valor_Name_Events));
		}
	}

	aux_name_events[0][0]->name = strdup("PAPI_L2_TCM");
	aux_name_events[0][0]->valor = PAPI_L2_TCM;
	//aux_name_events[0][1]->name = strdup("PAPI_L1_TCM");
	//aux_name_events[0][1]->valor = PAPI_L1_TCM;
	
	// Event PAPI_XXX
	//aux_name_events[1][0]->name = strdup("PAPI_XXX");
	//aux_name_events[1][0]->valor = PAPI_XXX;

	//aux_name_events[2][0]->name = strdup("PAPI_XXX");
	//aux_name_events[2][0]->valor = PAPI_XXX;
	
	//aux_name_events[3][0]->name = strdup("PAPI_XXX");
	//aux_name_events[3][0]->valor = PAPI_XXX;
	
	*name_events = aux_name_events;

	return (*length_events);
}

