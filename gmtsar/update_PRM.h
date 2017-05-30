#include <stdio.h>
#include "gmtsar.h"

void init_lookup_table();

int update_PRM(char* prmfile, char* entry, char* value);

int setvalue(struct PRM *prm, char* name, char* value);


char* get_PRM_string(char * prmfile, char * valuename);

int get_PRM_int(char * prmfile, char * valuename);

double get_PRM_double(char * prmfile, char * valuename);
