#include "update_PRM.h"
#include <stdio.h>

int main(int argc, char **argv) {

	char *szString = NULL;

	if (argc != 3 && argc != 4) {
		printf("Usage: get_PRM file.PRM param \n\nExample: get_PRM "
		       "IMG-HH-ALPSRP049040660-H1.0__A.PRM rshift\n");
		return 1;
	}

	szString = get_PRM(argv[1], argv[2]);
	printf("%s\n", szString);

	if (szString)
		free(szString);

	return 0;
}
