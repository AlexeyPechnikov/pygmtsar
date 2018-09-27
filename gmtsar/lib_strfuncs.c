#include "lib_functions.h"
#include <ctype.h>
#include <string.h>

/* Removes leading and trailing whitespaces from a string */
char *trimwhitespace(char *str) {
	char *end;

	// Trim leading space
	while (isspace((unsigned char)*str))
		str++;

	if (*str == 0) // All spaces?
		return str;

	// Trim trailing space
	end = str + strlen(str) - 1;
	while (end > str && isspace((unsigned char)*end))
		end--;

	// Write new null terminator
	*(end + 1) = 0;

	return str;
}
