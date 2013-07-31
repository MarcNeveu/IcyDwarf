/*
 * CHNOSZinfo.h
 *
 *  Created on: Apr 4, 2013
 *      Author: Marc Neveu (mneveu@asu.edu)
 *
 *      Equivalent of the info() function in CHNOSZ. Used to read the database OBIGT-2
 *      used by CHNOSZ (see CHNOSZ user manual at www.chnosz.net).
 *
 *      This function is modified from http://www.daniweb.com/software-development/c/
 *      threads/97843/parsing-a-csv-file-in-c#.
 *
 *      Note: this won't work in case of quoted entries containing commas:
 *      e.g., "1,2-diiodobutane".
 */

#ifndef CSVPARSE_H_
#define CSVPARSE_H_

#define MAXFLDSR 3500                         // > Number of rows in OBIGT.csv
#define MAXFLDSC 20                           // Number of columns in OBIGT-2.csv
#define MAXFLDSIZE 50                         // Character string field size

/* Note (MN 4/5/2013): it looks like the product MAXFLDSR*MAXFLDSC*MAXFLDSIZE cannot
 * be too large. I could not open OBIGT.csv simply because it had to many lines
 * (MAXFLDSR too large) unless I decreased MAXFLDSIZE from 200 to 50.
 */

void parse(char *record, char *delim, char arr[][MAXFLDSIZE], int *fldcnt);
struct OBIGTentry OBIGT(char *filename, int recordcnt);
int ID(char *filename, char *species, char *state);
void info(char *filename, int record);

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

struct OBIGTentry {
   	 char* name;
   	 char* abbrv;
   	 char* formula;
   	 char* state;
   	 char* ref1;
   	 char* ref2;
   	 char* date;
   	 float G;
   	 float H;
   	 float S;
   	 float Cp;
   	 float V;
   	 float a1_a;
   	 float a2_b;
   	 float a3_c;
   	 float a4_d;
   	 float c1_e;
   	 float c2_f;
   	 float omega_lambda;
   	 float z_T;
    };

/* Parse any line of a .csv file. Note: this won't work in case of quoted entries
 * containing commas: e.g., "1,2-diiodobutane". */

void parse(char *record, char *delim, char arr[][MAXFLDSIZE], int *fldcnt){
    int fld = 0;
	char *p = strtok(record,delim);

    while(p){
        strcpy(*(arr+fld),p);
		fld++;
		p = strtok('\0',delim);
	}

	*fldcnt = fld;
}

/* Open and read the entire OBIGT-2 database (304 records). Return only the
 * information for the record specified. This information is contained in an
 * OBIGTentry structure. Ideally, the entire OBIGT-2 database would be returned and
 * manipulated, but returning arrays by functions in C is a nightmare.
 */

struct OBIGTentry OBIGT(char *filename, int record){
	char tmp[1024] = {0x0};
	int fldcnt = 0;
	int recordcnt = 0;
	char arr[MAXFLDSR][MAXFLDSC][MAXFLDSIZE] = {{{0x0}}};
	struct OBIGTentry entry = {};

	FILE *in = fopen(filename,"r");

	if (in == NULL){
		printf("File open error\n");
		exit(0);
	}

	while(fgets(tmp,sizeof(tmp),in) != 0){    // Read a record
		parse(tmp,",",arr[recordcnt],&fldcnt); // Fill table arr with read records
	    recordcnt++;
	}

    fclose(in);

    entry.name = arr[record][0];
    entry.abbrv = arr[record][1];
    entry.formula = arr[record][2];
    entry.state = arr[record][3];
    entry.ref1 = arr[record][4];
    entry.ref2 = arr[record][5];
    entry.date = arr[record][6];
    sscanf(arr[record][7],"%f",&entry.G);
    sscanf(arr[record][8],"%f",&entry.H);
    sscanf(arr[record][9],"%f",&entry.S);
    sscanf(arr[record][10],"%f",&entry.Cp);
    sscanf(arr[record][11],"%f",&entry.V);
    sscanf(arr[record][12],"%f",&entry.a1_a);
    sscanf(arr[record][13],"%f",&entry.a2_b);
    sscanf(arr[record][14],"%f",&entry.a3_c);
    sscanf(arr[record][15],"%f",&entry.a4_d);
    sscanf(arr[record][16],"%f",&entry.c1_e);
    sscanf(arr[record][17],"%f",&entry.c2_f);
    sscanf(arr[record][18],"%f",&entry.omega_lambda);
    sscanf(arr[record][19],"%f",&entry.z_T);

   return entry;
}

/* Retrieve the line in the database corresponding to the species of interest.
 * Enter the name of the database, the formula for the species as it appears in the
 * "formula" column of the OBIGT databases, and the state ("cr", "cr2", "liq", "aq",
 * or "g").
 */

int ID(char *filename, char *species, char *state){
	char tmp[1024] = {0x0};
	int fldcnt = 0;
	int recordcnt = 0;
	char arr[MAXFLDSR][MAXFLDSC][MAXFLDSIZE] = {{{0x0}}};
	int record;
	char *state2;
	int i = 0;

	if (strcmp(state,"g") == 0)
		state2 = "gas";
	else
		state2 = state;

	FILE *in = fopen(filename,"r");

	if (in == NULL){
		printf("File open error\n");
		exit(0);
	}

	while(fgets(tmp,sizeof(tmp),in) != 0){    // Read a record
		parse(tmp,",",arr[recordcnt],&fldcnt); // Fill table arr with read records
	    recordcnt++;
	}

    fclose(in);

    for (i=0 ; i<MAXFLDSR ; i++) {
    	if (strcmp(*(*(arr+i)+2),species) == 0 && strcmp(*(*(arr+i)+3),state2) == 0){
    		record = i;
    		break;
    	}
    }
    return record;
}

/* Equivalent to the info() command in CHNOSZ */

void info(char *filename, int record){
	struct OBIGTentry entry = {};

	/* Optional: get current working directory to locate OBIGT.csv */

    char cwd[1024];
    if (getcwd(cwd, sizeof(cwd)) != NULL)
        printf("Current working directory: %s\n", cwd);
    else
        printf("getcwd() error\n");

    /* Retrieve info */

	entry = OBIGT(filename,record);
	char *aq = "aq";                             // To avoid an annoying warning

	/* Display info */

	printf("-----------------\n");
	printf("Info for %s (%s)\n",entry.name,entry.state);
	printf("name: %s\n",entry.name);
	printf("abbrv: %s\n",entry.abbrv);
	printf("formula: %s\n",entry.formula);
	printf("state: %s\n",entry.state);
	printf("ref1: %s\n",entry.ref1);
	printf("ref2: %s\n",entry.ref2);
	printf("date: %s\n",entry.date);
	printf("G: %f\n",entry.G);
	printf("H: %f\n",entry.H);
	printf("S: %f\n",entry.S);
	printf("Cp: %f\n",entry.Cp);
	printf("V: %f\n",entry.V);
	if (strcmp(entry.state,aq) == 0) {        // Solutes have HKF parameters
		printf("a1: %f\n",entry.a1_a);
		printf("a2: %f\n",entry.a2_b);
		printf("a3: %f\n",entry.a3_c);
		printf("a4: %f\n",entry.a4_d);
		printf("c1: %f\n",entry.c1_e);
		printf("c2: %f\n",entry.c2_f);
		printf("omega: %f\n",entry.omega_lambda);
		printf("Z: %f\n",entry.z_T);
	}
	else {                                    // Solid, liquid, or gas have Maier-
		printf("a: %f\n",entry.a1_a);         // Kelley parameters
		printf("b: %f\n",entry.a2_b);
		printf("c: %f\n",entry.a3_c);
		printf("d: %f\n",entry.a4_d);
		printf("e: %f\n",entry.c1_e);
		printf("f: %f\n",entry.c2_f);
		printf("lambda: %f\n",entry.omega_lambda);
		printf("T: %f\n",entry.z_T);
	}
	printf("-----------------\n");
}

#endif /* CSVPARSE_H_ */
