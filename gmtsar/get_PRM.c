#include <stdio.h>
#include "update_PRM.h"


int main(int argc, char **argv)
{

	char* szString=NULL;
	int nDecimals = 0;

	if(argc != 3 && argc != 4){
		printf("Usage: get_PRM file.PRM param \n\nExample: get_PRM IMG-HH-ALPSRP049040660-H1.0__A.PRM rshift\n");
		return 1;
	}
	
	if(argc == 4)
		nDecimals = atoi(argv[3]);
	else
		nDecimals = 2;


	szString = get_PRM(argv[1], argv[2], nDecimals);
	printf("%s\n",szString);

	if(szString)
		free(szString);
	
	return 0;
}
