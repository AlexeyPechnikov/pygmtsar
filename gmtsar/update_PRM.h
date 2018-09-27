#ifndef UPDATE_PRM_H
#define UPDATE_PRM_H

#include "gmtsar.h"
#include <stdio.h>

void init_lookup_table();

int update_PRM_sub(char *prmfile, char *entry, char *value);

int setvalue(struct PRM *prm, char *name, char *value);

char *get_PRM_string(char *prmfile, char *valuename);

int get_PRM_int(char *prmfile, char *valuename);

double get_PRM_double(char *prmfile, char *valuename);

char *get_PRM(char *prmfile, char *valuename);

#endif /* UPDATE_PRM_H */
