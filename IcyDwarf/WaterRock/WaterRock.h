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

	char *idi = (char*)malloc(1024);
	idi[0] ='\0';
	if (release == 1) {
		idi = path;
		strncat(idi,path,strlen(path)-38);
	}
	else if (cmdline == 1) {
		idi = path;
		strncat(idi,path,strlen(path)-40);
	}
	else strncat(idi,path,strlen(path)-38);
	strcat(idi,"Codes/PHREEQC/phreeqc-2.18.3/bin/phreeqc");

	printf("%s\n",idi);
	system(idi);

	free(idi);
	return 0;
}

#endif /* WATERROCK_H_ */
