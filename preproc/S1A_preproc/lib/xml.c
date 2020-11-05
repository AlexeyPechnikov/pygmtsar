/***************************************************************************
 * Creator:  Xiaohua(Eric) XU                                              *
 *           (Scripps Institution of Oceanography)                         *
 * Date   :  04/06/2015                                                    *
 ***************************************************************************/

/***************************************************************************
 * Modification history:                                                   *
 *                                                                         *
 * Date   :  DTS made the set of xml functions into a library              *
 * Date   :  EX made changed the way of getting tree (more robust)         *
 ***************************************************************************/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "lib_defs.h"

/* global variables needed for the xml library */
int N = 0;
int MAX_TREE_SIZE = 600000; // size of the tree in maximum
int MAX_CHAR_SIZE = 60000;  // size of char arrays in maximum
char STR[4000][60000];

int search_tree(tree *list, char *str, char *s_out, int type, int loc, int num) {
	/***************************************************************************
	    list     - pointer to the tree structure
	    str_out  - pointer to the output string(everything will be output as
	 string). str      - is the string used for searching in the shape of
	 /lev1_name/lev2_name/.../levN_name/ Note all the names should be seperated by
	 '/', including the start and end with '/'. type     - choice of return where
	               1 string to str_out
	               2 converted date to str_out with the format yyyyddd.frac_of_day
	               3 values in the names, for example <orbitLists 'count = 31'>
	 will return 31

	    For names occurred in multiple locations, use loc and num to return the
	 one at correct location. loc      - the term to be repeatedly found, num -
	 index of array. For example, to search
	 /product/OrbitInformation/OrbitList/orbit/position/x/ and there are multiple
	 state vectors under Orbitlist, set loc = 4, the 4th value of x

	  After reading parameters in string format, use str2double to
	  convert to double, use str2ints to convert to integer array, or use
	  str2dbs to convert to double array.
	 ***************************************************************************/
	// search the num-th target at loc in str in the tree
	long int ct = 0;
	int n = 1, i, j, j1, j2;
	char s_name[200], *tmp_name, tmp_num[200];
	while (strlocate(str, '/', n) != -1) {
		n++;
	}
	n = n - 2;

	tmp_name = (char *)malloc(MAX_CHAR_SIZE * sizeof(char));

	for (i = 0; i < n; i++) {
		j1 = strlocate(str, '/', i + 1);
		j2 = strlocate(str, '/', i + 2);
		strasign(s_name, str, j1 + 1, j2 - 1);

		if (strncmp(list[ct].name, "OutOfSpace", 10) == 0) {
			cat_nums(tmp_num, list[ct].name);
			strcpy(tmp_name, STR[(int)str2double(tmp_num)]);
		}
		else {
			strcpy(tmp_name, list[ct].name);
		}

		while (ct < MAX_TREE_SIZE && ct >= 0 && strncmp(tmp_name, s_name, strlen(s_name)) != 0) {
			ct = list[ct].sibr;
			if (strncmp(list[ct].name, "OutOfSpace", 10) == 0) {
				cat_nums(tmp_num, list[ct].name);
				strcpy(tmp_name, STR[(int)str2double(tmp_num)]);
			}
			else {
				strcpy(tmp_name, list[ct].name);
			}
		}

		if (ct >= MAX_TREE_SIZE || ct < 0 || list[ct].firstchild == -1) {
			fprintf(stderr, "Unable to find designated string...after %ld searches..\n", ct);
			return (-1);
		}

		if (loc == i + 1) {
			if (num > 1) {
				for (j = 1; j < num; j++) {
					if (list[ct].sibr != -1) {
						ct = list[ct].sibr;
					}
					else if (list[ct].sibr == -1) {
						fprintf(stderr, "Unable to find designated string...\n");
						return (-1);
					}
				}
			}
		}
		ct = list[ct].firstchild;
	}

	if (type == 1) {
		if (strncmp(list[ct].name, "OutOfSpace", 10) == 0) {
			cat_nums(s_name, list[ct].name);
			strcpy(s_out, STR[(int)str2double(s_name)]);
		}
		else {
			strcpy(s_out, list[ct].name);
		}
	}
	else if (type == 2) {
		cat_nums(s_name, list[ct].name);
		str_date2JD(s_out, s_name);
	}
	else if (type == 3) {
		cat_nums(s_out, list[list[ct].parent].name);
	}
	free(tmp_name);
	return (list[ct].parent);
}

int get_tree(FILE *fp, tree *list, int num_parse) {

	/*
	  fp        - file pointer to the xml file
	  list      - structure to hold the entire xml file
	  num_parse - number of header lines to skip

	  The xml file is required to have the names in pairs, either as
	      <name_of_par>value_of_par</name_of_par>
	  or as
	      <name_of_par>
	        ...
	      </name_of_par>
	  The current code can also deal with some not well shaped conditions like
	      <name_of_par.../name_of_par>
	  or lines not indented correctly like
	      <name_of_par>
	        ...
	        <name_of_par>
	  For questions, please send to xix016@ucsd.edu
	*/

	char *buffer;
	char tmp_char[200], tmp_s[200], *tmp_c;
	int i1, i2, j1, j2, have_slash;
	long int count = 0;
	// int *num_space;
	long int level[100] = {-1}, lev_ct = 0;
	char lev_rec[100][200];

	buffer = (char *)malloc(MAX_CHAR_SIZE * sizeof(char));
	tmp_c = (char *)malloc(MAX_CHAR_SIZE * sizeof(char));

	for (i1 = 0; i1 < 100; i1++) {
		strcpy(lev_rec[i1], "CLOSED");
	}

	// num_space = (int *)malloc(aprox_size*5*sizeof(int));

	// fprintf(stderr," %d \n",sizeof(buffer));

	for (i1 = 0; i1 < num_parse; i1++) {
		fgets(buffer, MAX_CHAR_SIZE * sizeof(char), fp);
	}

	while (fgets(buffer, MAX_CHAR_SIZE * sizeof(char), fp) != NULL) {
		// num_space[count] = space_count(buffer);
		i1 = strlocate(buffer, '<', 1);
		j1 = strlocate(buffer, '>', 1);
		i2 = strlocate(buffer, '<', 2);
		j2 = strlocate(buffer, '>', 2);

		if (i1 < 0 || j1 < 0) {
			fprintf(stderr, "Not an well formatted XML file...%d, %d\n", i1, j1);
			return (-1);
		}
		else if (buffer[i1 + 1] == '/') {
			have_slash = 1;
			strasign(tmp_char, buffer, i1 + 2, j1 - 1);
		}
		else {
			have_slash = 0;
			strasign(tmp_char, buffer, i1 + 1, j1 - 1);
		}

		// initinate the numbers
		list[count].sibr = -1;
		list[count].sibl = -1;
		list[count].parent = -1;
		list[count].firstchild = -1;

		// create tree
		// first node
		if (count == 0) {
			strcpy(list[count].name, tmp_char);
			level[lev_ct] = count;
			strcpy(lev_rec[lev_ct], tmp_char);
			if (i2 != -1 || (buffer[j1 - 1] == '/')) {
				list[count].firstchild = count + 1;
				create_child(list, buffer, j1 + 1, i2 - 1, count++);
				strcpy(lev_rec[lev_ct], "CLOSED");
			}
		}
		// child
		else if (count != 0 && have_slash == 0 && strncmp(lev_rec[lev_ct], "CLOSED", 6) != 0) {
			lev_ct++;
			strcpy(list[count].name, tmp_char);
			strcpy(lev_rec[lev_ct], tmp_char);
			list[count].parent = level[lev_ct - 1];
			list[count].sibl = -1;
			list[count].sibr = -1;
			list[count].firstchild = -1;
			list[level[lev_ct - 1]].firstchild = count;
			level[lev_ct] = count;
			if (i2 != -1 || (buffer[j1 - 1] == '/')) {
				list[count].firstchild = count + 1;
				create_child(list, buffer, j1 + 1, i2 - 1, count++);
				// fprintf(stderr,"CHILD: %s,LVL: %ld\n",tmp_char,lev_ct);
				strcpy(lev_rec[lev_ct], "CLOSED");
			}
		}
		// sibling
		else if (count != 0 && have_slash == 0 && strncmp(lev_rec[lev_ct], "CLOSED", 6) == 0) {
			strcpy(list[count].name, tmp_char);
			strcpy(lev_rec[lev_ct], tmp_char);
			list[count].sibl = level[lev_ct];
			list[count].sibr = -1;
			list[count].firstchild = -1;
			list[count].parent = list[level[lev_ct]].parent;
			list[level[lev_ct]].sibr = count;
			level[lev_ct] = count;
			if (i2 != -1 || (buffer[j1 - 1] == '/')) {
				list[count].firstchild = count + 1;
				create_child(list, buffer, j1 + 1, i2 - 1, count++);
				strcpy(lev_rec[lev_ct], "CLOSED");
			}
		}
		// go to parent level

		else if (count != 0 && have_slash == 1) {
			// fprintf(stderr,"%s\n",tmp_char);
			if (strncmp(lev_rec[lev_ct - 1], "OutOfSpace", 10) == 0) {
				cat_nums(tmp_s, lev_rec[lev_ct - 1]);
				strcpy(tmp_c, STR[(int)str2double(tmp_s)]);
			}
			else {
				strcpy(tmp_c, lev_rec[lev_ct - 1]);
			}
			if (strncmp(tmp_char, tmp_c, strlen(tmp_char)) == 0) {
				strcpy(lev_rec[lev_ct - 1], "CLOSED");
				lev_ct--;
			}
		}

		// printf("%s",buffer);
		// printf("%d %d %d\n",num_space[count],count,lev_ct);
		count++;
	}
	free(buffer);
	free(tmp_c);
	// fclose(fp);
	// free(num_space);
	return (1);
}

int space_count(char *str) {
	// count the number of spaces in front of each line
	int i, j;
	j = 0;
	for (i = 0; i < strlen(str); i++) {
		if (str[i] == ' ') {
			j = j + 1;
			continue;
		}
		else if (str[i] == '\t') {
			j = j + 2;
			continue;
		}
		else {
			break;
		}
	}
	return (j);
}

int itoa_xml(int d, char *buf, int base) {
	char *p = buf;
	char *p1, *p2;
	unsigned long ud = d;
	int divisor = 10;

	/* If %d is specified and D is minus, put `-' in the head.  */
	if (base == 'd' && d < 0) {
		*p++ = '-';
		buf++;
		ud = -d;
	}
	else if (base == 'x') {
		divisor = 16;
	}

	/* Divide UD by DIVISOR until UD == 0.  */
	do {
		int remainder = ud % divisor;

		*p++ = (remainder < 10) ? remainder + '0' : remainder + 'a' - 10;
	} while (ud /= divisor);

	/* Terminate BUF.  */
	*p = 0;

	/* Reverse BUF.  */
	p1 = buf;
	p2 = p - 1;
	while (p1 < p2) {
		char tmp = *p1;
		*p1 = *p2;
		*p2 = tmp;
		p1++;
		p2--;
	}
	return (1);
}

int strasign(char *str_out, char *str, int n1, int n2) {
	// asign n1-n2 of str to str_out
	int i;
	if (n1 > n2) {
		return (-1);
	}

	if (n2 - n1 > 199) {
		// fprintf(stderr,"OutOfSpace %d, n1: %d, n2 %d\n",N,n1,n2);
		strcpy(str_out, "OutOfSpace");
		char c[100];
		itoa_xml(N, c, 'd');
		strcat(str_out, c);
		// strcpy(STR[N],str);
		// fprintf(stderr,"%s\n",str_out);

		// STR[N] = (char *)malloc(MAX_CHAR_SIZE*sizeof(char));

		for (i = n1; i <= n2; i++) {
			STR[N][i - n1] = str[i];
		}
		STR[N][i - n1] = '\0';
		// fprintf(stderr,"%s\n",STR[N]);
		N++;
		return (1);
	}

	for (i = n1; i <= n2; i++) {
		str_out[i - n1] = str[i];
	}
	str_out[n2 - n1 + 1] = '\0';
	return (1);
}

int strlocate(char *str, int c, int n) {
	// locate the n-th c in str
	int i, j = 0;
	for (i = 0; i < strlen(str); i++) {
		if (str[i] == (char)c) {
			j++;
			if (j == n) {
				return (i);
			}
		}
	}
	return (-1);
}

int create_child(tree *T, char *str, int i, int j, int ct) {
	// create the firstchild of a tree-element
	if (j >= i) {
		strasign(T[ct + 1].name, str, i, j);
	}
	else {
		T[ct + 1].name[0] = '\0';
	}
	T[ct + 1].sibr = -1;
	T[ct + 1].sibl = -1;
	T[ct + 1].parent = ct;
	T[ct + 1].firstchild = -1;
	// num_space[ct+1] = num_space[ct];
	return (1);
}

int show_tree(tree *T, int ct, int lvl) {
	// print out the tree
	if (ct == 0) {
		printf("In the brackets the numbers are (count, child, sibling)\n");
	}
	int i;
	for (i = 0; i < 2 * lvl; i++) {
		putchar(' ');
	}
	if (T[ct].sibr == -1 && T[ct].firstchild == -1) {
		putchar('=');
	}

	if (strncmp(T[ct].name, "OutOfSpace", 10) == 0) {
		char c[100];
		cat_nums(c, T[ct].name);
		printf("%s   (%d,%d,%d)\n", STR[(int)str2double(c)], ct, T[ct].firstchild, T[ct].sibr);
	}
	else {
		printf("%s   (%d,%d,%d)\n", T[ct].name, ct, T[ct].firstchild, T[ct].sibr);
	}

	if (T[ct].firstchild != -1) {
		show_tree(T, T[ct].firstchild, lvl + 1);
	}
	if (T[ct].sibr != -1) {
		show_tree(T, T[ct].sibr, lvl);
	}
	return (1);
}

int cat_nums(char *str_out, char *str) {
	// cat out the numbers in str to str_out
	int i = 0, j = 0;
	while (str[i] != '\0') {
		if (str[i] >= '0' && str[i] <= '9') {
			str_out[j++] = str[i];
		}
		i++;
	}
	str_out[j] = '\0';
	return (j);
}

double date2MJD(int yr, int mo, int day, int hr, int min, double sec) {
	/* convert to date to MJD */
	double part1, part2;
	double MJD;

	part1 = 367 * ((double)yr) - floor(7 * ((double)yr + floor(((double)mo + 9) / 12.0)) / 4.0) + floor(275 * (double)mo / 9.0) +
	        (double)day;
	part2 = -678987 + ((sec / 60.0 + (double)min) / 60.0 + (double)hr) / 24.0;

	MJD = part1 + part2;

	return MJD;
}

int str_date2JD(char *str_JD, char *str_date) {
	double str2double(char *);
	int yr, mo, day, hr, min, doy;
	double sec = 0, MJDday, MJDyr, MJDfrac;
	char tmp[30];

	strasign(tmp, str_date, 0, 3);
	yr = (int)str2double(tmp);
	strasign(tmp, str_date, 4, 5);
	mo = (int)str2double(tmp);
	strasign(tmp, str_date, 6, 7);
	day = (int)str2double(tmp);
	strasign(tmp, str_date, 8, 9);
	hr = (int)str2double(tmp);
	strasign(tmp, str_date, 10, 11);
	min = (int)str2double(tmp);
	strasign(tmp, str_date, 12, 13);
	sec = sec + str2double(tmp);
	strasign(tmp, str_date, 14, 19);
	sec = sec + str2double(tmp) / 1000000.0;
	// printf("%d %d %d %d %d %.10f\n",yr,mo,day,hr,min,sec);
	MJDyr = date2MJD(yr, 1, 1, 0, 0, 0);
	MJDday = date2MJD(yr, mo, day, 0, 0, 0);
	MJDfrac = (((hr * 60.) + min) * 60 + sec) / 86400;
	doy = (int)(MJDday - MJDyr + .1);
	sprintf(str_JD, "%.12f", doy + MJDfrac);

	return (1);
}

double str2double(char *str) {
	int i, n, m;
	double value = 0.0, value1 = 0.0, value2 = 0.0, sgn = 1.0;
	char tmp1[100], tmp2[100], tmp[100], str_tmp[100];

	strasign(str_tmp, str, 0, strlen(str));

	// decide the sign
	if (str_tmp[0] == '-' || str_tmp[0] == '+') {
		if (str_tmp[0] == '-') {
			sgn = -1.0;
		}
		strasign(tmp, str_tmp, 1, strlen(str_tmp));
		strasign(str_tmp, tmp, 0, strlen(tmp));
	}

	// decide where it is sci form
	if (strlocate(str_tmp, 'e', 1) != -1 || strlocate(str_tmp, 'E', 1) != -1) {
		n = strlocate(str_tmp, 'e', 1);
		if (n == -1) {
			n = strlocate(str_tmp, 'E', 1);
		}
		strasign(tmp2, str_tmp, n + 1, strlen(str)); // exponential part
		strasign(tmp1, str_tmp, 0, n - 1);           // digits part
	}
	else {
		strasign(tmp1, str_tmp, 0, strlen(str_tmp)); // digits part
	}
	// decide whether it has fraction
	n = strlocate(tmp1, '.', 1);
	if (n != -1) {
		strasign(tmp, tmp1, 0, n - 1);
		m = strlen(tmp);
		for (i = 0; i < m; i++) {
			value1 = value1 + (double)((int)tmp[i] - 48) * pow(10.0, (double)(m - i - 1));
		}
		m = strlen(tmp1);
		strasign(tmp, tmp1, n + 1, m);
		m = strlen(tmp);
		for (i = 0; i < m; i++) {
			value2 = value2 + (double)((int)tmp[i] - 48) * pow(10.0, (double)(-i - 1));
		}
		// puts(tmp2);
		// fprintf(stderr,"%.12f    %.12f %.12f\n",value1,value2,str2double(tmp2));
		value = value1 + value2;
	}
	else {
		m = strlen(tmp1);
		for (i = 0; i < m; i++) {
			value = value + (double)((int)tmp1[i] - 48) * pow(10.0, (double)(m - i - 1));
		}
	}

	if (strlocate(str_tmp, 'e', 1) != -1 || strlocate(str_tmp, 'E', 1) != -1) {
		value = value * pow(10.0, str2double(tmp2));
	}

	return (value * sgn);
}

int str2ints(int *a, char *c) {
	int i = 0, n1 = 0, n2 = 0, len;
	char tmp_c[100];

	len = strlen(c);

	for (n2 = 0; n2 < len; n2++) {
		if (strncmp(&c[n2], " ", 1) == 0) {
			strasign(tmp_c, c, n1, n2 - 1);
			a[i] = (int)str2double(tmp_c);
			i++;
			n1 = ++n2;
		}
	}

	n2 = len;

	if (n2 > n1 + 1) {
		strasign(tmp_c, c, n1, n2);
		a[i] = (int)str2double(tmp_c);
		i++;
	}
	/*
	 while((n2 = strlocate(c,' ',i+1)) != -1){
	 strasign(tmp_c,c,n1,n2-1);
	 a[i] = (int)str2double(tmp_c);
	 i++;
	 n1 = n2+1;
	 }
	 */
	return (i);
}

int str2dbs(double *a, char *c) {
	int i = 0, n1 = 0, n2 = 0, len;
	char tmp_c[100];

	len = strlen(c);
	for (n2 = 0; n2 < len; n2++) {
		if (c[n2] == ' ' || c[n2] == '\0' || c[n2] == '\n') {
			strasign(tmp_c, c, n1, n2 - 1);
			a[i] = str2double(tmp_c);
			i++;
			n1 = ++n2;
		}
	}

	n2 = len;

	if (n2 > n1 + 1) {
		strasign(tmp_c, c, n1, n2);
		a[i] = str2double(tmp_c);
		i++;
	}
	/*
	 while((n2 = strlocate(c,' ',i+1)) != -1){
	 strasign(tmp_c,c,n1,n2-1);
	 a[i] = (int)str2double(tmp_c);
	 i++;
	 n1 = n2+1;
	 }
	 */
	return (i);
}

int null_MEM_STR() {
	int i;
	for (i = 0; i < 100; i++) {
		STR[i][0] = '\0';
	}
	N = 0;
	return (1);
}

int prefix_str(char *a, char *b) {
	// put b infront of a
	char *str;
	str = (char *)malloc((strlen(a) + strlen(b)) * 2 * sizeof(char));

	strcpy(str, b);
	strcat(str, a);
	strcpy(a, str);

	free(str);
	return (1);
}

int print_space(int j, FILE *fp) {
	int i;

	for (i = 0; i < j; i++)
		fprintf(fp, "%s", " ");

	return (1);
}

/*
int find_assembly(struct tree **T, int ct, int ii) {

    char str[2000],str_out[60000];
    int jj;

    str[0] = '\0';
    jj = ct;
    while(T[ii][jj].parent != -1) {
        prefix_str(str,"/");
        sscanf(T[ii][jj].name,"%s ",str_out);
        prefix_str(str,str_out);
        jj = T[ii][jj].parent;
    }
    prefix_str(str,"/");
    sscanf(T[ii][jj].name,"%s ",str_out);
    prefix_str(str,str_out);
    prefix_str(str,"/");

    //fprintf(stderr,"%s\n",str);


    jj = search_tree(T[ii],str,str_out,1,0,1);
    //search_tree(tree *list, char *str, char *s_out, int type, int loc, int
num)
    //fprintf(stderr,"%d found...\n",jj);
    return(jj);

}
*/

int print_tree(struct tree *T, int ct, int mode, FILE *fp) {

	char str[200];

	if (strncmp(T[ct].name, "OutOfSpace", 10) == 0) {
		char c[100];
		cat_nums(c, T[ct].name);
		if (mode == 1) {
			fprintf(fp, "<%s>\n", STR[(int)str2double(c)]);
		}
		else if (mode == 2) {
			sscanf(STR[(int)str2double(c)], "%s ", str);
			fprintf(fp, "</%s>\n", str);
		}
		else {
			sscanf(STR[(int)str2double(c)], "%s ", str);
			if (strncmp(T[T[ct].firstchild].name, "OutOfSpace", 10) == 0) {
				char c2[100];
				cat_nums(c2, T[T[ct].firstchild].name);
				fprintf(fp, "<%s>%s</%s>\n", STR[(int)str2double(c)], STR[(int)str2double(c2)], str);
			}
			else {
				fprintf(fp, "<%s>%s</%s>\n", STR[(int)str2double(c)], T[T[ct].firstchild].name, str);
			}
		}
	}
	else {
		if (mode == 1) {
			fprintf(fp, "<%s>\n", T[ct].name);
		}
		else if (mode == 2) {
			sscanf(T[ct].name, "%s ", str);
			fprintf(fp, "</%s>\n", str);
		}
		else {
			sscanf(T[ct].name, "%s ", str);
			if (strncmp(T[T[ct].firstchild].name, "OutOfSpace", 10) == 0) {
				char c2[100];
				cat_nums(c2, T[T[ct].firstchild].name);
				fprintf(fp, "<%s>%s</%s>\n", T[ct].name, STR[(int)str2double(c2)], str);
			}
			else {
				fprintf(fp, "<%s>%s</%s>\n", T[ct].name, T[T[ct].firstchild].name, str);
			}
		}
	}

	return (1);
}

int assemble_trees(int nfiles, struct tree **T, int ct, int lvl, FILE *fp) {

	int j, k;

	if (ct == 0)
		fprintf(fp, "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n");

	print_space(2 * lvl, fp);
	/*
	    if(T[0][ct].sibr == -1 && T[0][ct].firstchild == -1) {
	        putchar('=');
	    }
	*/

	if (T[0][ct].firstchild != -1 && T[0][T[0][ct].firstchild].firstchild == -1) {
		print_tree(T[0], ct, 3, fp);
		/*
		        for (i=1;i<nfiles;i++) {
		            j = find_assembly(T,ct,i);
		            if (j>=0) {
		                if
		   (strncmp(T[0][ct].name,T[i][j].name,strlen(T[0][ct].name)-1) == 0 &&
		   strncmp(T[0][T[0][ct].firstchild].name,T[i][T[i][j].firstchild].name,strlen(T[0][T[0][ct].firstchild].name)-1)
		   == 0) mark = 0; else mark = 1; if (mark == 1) { print_space(2*lvl,fp);
		                    print_tree(T[i],j,3,fp);
		                }
		            }
		        }
		*/

		if (T[0][ct].sibr != -1) {
			assemble_trees(nfiles, T, T[0][ct].sibr, lvl, fp);
		}
		else {
			j = ct;
			k = 1;
			while (T[0][j].sibr == -1 && j != 0) {
				print_space(2 * (lvl - k), fp);
				print_tree(T[0], T[0][j].parent, 2, fp);
				j = T[0][j].parent;
				k++;
			}
		}
	}
	else {
		print_tree(T[0], ct, 1, fp);
		/*
		        for (i=1;i<nfiles;i++) {
		            j = find_assembly(T,ct,i);
		            if (j>=0) {
		                if
		   (strncmp(T[0][ct].name,T[i][j].name,strlen(T[0][ct].name)-1) == 0) mark =
		   0; else mark = 1; if (mark == 1) { print_space(2*lvl,fp);
		                    print_tree(T[i],j,1,fp);
		                }
		            }
		        }
		*/

		if (T[0][ct].firstchild != -1) {
			assemble_trees(nfiles, T, T[0][ct].firstchild, lvl + 1, fp);
		}

		if (T[0][ct].sibr != -1) {
			assemble_trees(nfiles, T, T[0][ct].sibr, lvl, fp);
		}
		if (T[0][ct].firstchild == -1 && T[0][ct].sibr == -1) {
			j = ct;
			k = 1;
			while (T[0][j].sibr == -1 && j != 0) {
				print_space(2 * (lvl - k), fp);
				print_tree(T[0], T[0][j].parent, 2, fp);
				j = T[0][j].parent;
				k++;
			}
		}
	}

	return (1);
}

/*
add_branch(int n, struct tree **T, str, ct) {
    
    int i;
    char str_out[60000];
    i = search_tree(T[0],str,str_out,1,0,1);


    return(1);

}
*/
/*
int add_index(struct tree *T, int ct, int ii) {

    int i;

    //printf("hahaha %d\n",ii);
    for (i=0;i<ct;i++){
        if(T[i].parent != -1) T[i].parent = T[i].parent+ii;
        if(T[i].sibl != -1) T[i].sibl = T[i].sibl+ii;
        if(T[i].sibr != -1) T[i].sibr = T[i].sibr+ii;
        if(T[i].firstchild != -1) T[i].firstchild = T[i].firstchild+ii;
    }

    //printf("hahaha\n");

    return(1);
}

double time2double(char *str){
    double k;
    char s_name[200],s_out[200];

    cat_nums(s_name,str);
    str_date2JD(s_out, s_name);

    k = str2double(s_out);

    return(k);

}
*/
/*
int edit_tree(int nfiles,int nlmx,struct tree **T){

    // note this is only editing firstchild and sibr, thus its not really
getting a good tree structure

    int ii1,jj1,kk1,ct,nn1,ii2,jj2,kk2,kk3,nn2,qq;
    char tmp_c[60000],s_name[200],s_out[200];
    double t1,t2;
  
    ct = 0;
    while (T[0][ct].sibr != -1 || T[0][ct].firstchild != -1) {
        if (T[0][ct].sibr != -1) ct = T[0][ct].sibr;
        else ct = T[0][ct].firstchild;
    }
    printf("Original(first) xml has %d tree elements.\n",ct);

    for (qq=1;qq<nfiles;qq++) {
        ct = 0;
        while (T[qq][ct].sibr != -1 || T[qq][ct].firstchild != -1) {
            if (T[qq][ct].sibr != -1) ct = T[qq][ct].sibr;
            else ct = T[qq][ct].firstchild;
        }
        add_index(&T[0][nlmx*qq*5],ct+1,nlmx*qq*5);
    }

    // start editing the useful information to create one xml_tree
    for (qq=1;qq<nfiles;qq++) {

        // burst List
        ii1 = search_tree(T[0],"/product/swathTiming/burstList/",tmp_c,3,0,1);
        nn1 = (int)str2double(tmp_c);
        kk1 = T[0][ii1].firstchild;
        for (jj1 = 1;jj1<nn1;jj1++) kk1 = T[0][kk1].sibr;
        ii2 = search_tree(T[qq],"/product/swathTiming/burstList/",tmp_c,3,0,1);
        T[0][kk1].sibr = T[0][ii2+nlmx*qq*5].firstchild;
        nn2 = (int)str2double(tmp_c);
        kk2 = T[0][ii2+nlmx*qq*5].firstchild;
        for (jj2 = 1;jj2<nn2;jj2++) {
            T[0][kk2].parent = ii1;
            kk2 = T[0][kk2].sibr;
        }
        T[0][kk2].parent = ii1;
        sprintf(tmp_c,"burstList count=\"%d\"",nn1+nn2);
        strcpy(T[0][ii1].name,tmp_c);

        // orbit List
        ii1 =
search_tree(T[0],"/product/generalAnnotation/orbitList/",tmp_c,3,0,1); nn1 =
(int)str2double(tmp_c); kk1 = T[0][ii1].firstchild; for (jj1 = 1;jj1<nn1;jj1++)
kk1 = T[0][kk1].sibr; t1 =
time2double(T[0][T[0][T[0][kk1].firstchild].firstchild].name); ii2 =
search_tree(T[qq],"/product/generalAnnotation/orbitList/",tmp_c,3,0,1); nn2 =
(int)str2double(tmp_c); printf("%d %lf\n",nn1,t1); t2 = t1-1;ct = 0; kk2 =
T[0][ii2+nlmx*qq*5].firstchild;
        t2=time2double(T[0][T[0][T[0][kk2].firstchild].firstchild].name);
        while (t2<t1) {
            kk2 = T[0][kk2].sibr;
            t2=time2double(T[0][T[0][T[0][kk2].firstchild].firstchild].name);
            ct++;
        }
        T[0][kk1].sibr = kk2;
        kk2 = T[0][ii2+nlmx*qq*5].firstchild;
        for (jj2 = 1;jj2<nn2-ct;jj2++) {
            T[0][kk2].parent = ii1;
            kk2 = T[0][kk2].sibr;
        }
        T[0][kk2].parent = ii1;
        sprintf(tmp_c,"orbitList count=\"%d\"",nn1+nn2-ct);
        strcpy(T[0][ii1].name,tmp_c);
    

    }

    return(1);

}


*/
