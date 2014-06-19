/*
 * WaterRock.h
 *
 *  Created on: Jun 19, 2014
 *      Author: Marc Neveu (mneveu@asu.edu)
 *
 *      Header file in development. Calls a geochemical code (PHREEQC at the moment).
 */

#ifndef WATERROCK_H_
#define WATERROCK_H_

#include <stdio.h>
#include <stdlib.h>

#include "../IcyDwarf.h"

int WaterRock(char path[1024]);

int WaterRock(char path[1024]) {

	char *bin = (char*)malloc(1024);
	char *infile = (char*)malloc(1024);
	char *outfile = (char*)malloc(1024);
	char *dbase = (char*)malloc(1024);

	bin[0] ='\0';
	if (release == 1) {
		bin = path;
		strncat(bin,path,strlen(path)-38);
	}
	else if (cmdline == 1) {
		bin = path;
		strncat(bin,path,strlen(path)-40);
	}
	else strncat(bin,path,strlen(path)-38);
	strcat(bin,"Codes/PHREEQC/phreeqc-2.18.3/bin/phreeqc ");

	infile[0] = '\0';
	strncat(infile,bin,strlen(bin)-12);
	strcat(infile,"examples/ex1 ");

	outfile[0] = '\0';
	strncat(outfile,bin,strlen(bin)-12);
	strcat(outfile,"examples/ex1.out ");

	dbase[0] = '\0';
	strncat(dbase,bin,strlen(bin)-12);
	strcat(dbase,"database/phreeqc.dat");

	strcat(bin, infile);
	strcat(bin, outfile);
	strcat(bin, dbase);

	system(bin);

	free(bin);
	free(infile);
	free(outfile);
	free(dbase);

	return 0;
}

#endif /* WATERROCK_H_ */
