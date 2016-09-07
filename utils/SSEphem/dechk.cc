//This file is part of SSEphem.

//  SSEphem is a collection of C/C++ source files that provide an interface
//  to the JPL Solar System Ephemeris.
//  Copyright (C) 2009 Associated Universities, Inc. Washington DC, USA

//  SSEphem is free software: you can redistribute it and/or modify
//  it under the terms of the GNU Lesser General Public License as published
//  by the Free Software Foundation, either version 3 of the License, or
//  (at your option) any later version.

//  SSEphem is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU Lesser General Public License for more details.

//  You should have received a copy of the GNU Lesser General Public License
//  along with SSEphem.  If not, see <http://www.gnu.org/licenses/>.

#include <stdio.h>
#include <string.h>
#include "jpl_eph.h"

char *delete_trailing_blanks(char *str)
{
    for (int i = strlen(str) - 1; i >= 0; i--) {
        if (str[i] == ' ') {
            str[i] = '\0';
        } else {
            break;
        }
    }
    return str;
}

void print_header(de_headerT hdr)
{
    printf ("%s\n", delete_trailing_blanks(hdr.title1));
    printf ("%s\n", delete_trailing_blanks(hdr.title2));
    printf ("%s\n", delete_trailing_blanks(hdr.title3));
    printf ("Constant names 1, 2, 199, and 200\n");
    printf ("%s %s .... %s %s\n", hdr.constant_name[0], hdr.constant_name[1],
                             hdr.constant_name[198], hdr.constant_name[199]);
    printf ("%5.3f  %5.3f  %5.3f\n", hdr.epoch1, hdr.epoch2, hdr.interval);
    printf ("%d constants\n", hdr.num_constants);
    printf ("AU = %24.16e\n", hdr.au);
    printf ("EMRAT = %24.16e\n", hdr.emrat);
    printf ("ipt values:\n");
    for (int j = 0; j < 3; j++) {
        for (int i = 0; i < 12; i++) {
            printf ("%6d", hdr.ipt[j][i]);
        }
        printf ("\n");
    }
    printf ("DE %d\n", hdr.DEnumber);
    printf ("lpt = %d %d %d\n", hdr.lpt[0], hdr.lpt[1], hdr.lpt[2]);
    printf ("Constants 1 & 2:     %24.16e  %24.16e\n", hdr.constant_value[0],
                                       hdr.constant_value[1]);
    printf ("Constants 199 & 200: %24.16e  %24.16e\n", hdr.constant_value[198],
                                       hdr.constant_value[199]);
    printf ("%5d coefficients\n", hdr.num_coeff);
}

void main(int argc, char *argv[])
{
    if (argc < 2) {
        printf ("Usage: dechk file_name\n");
        return;
    }
    FILE *strm = fopen(argv[1], "rb");
    de_headerT hdr;
    fread(&hdr, sizeof(hdr), 1, strm);
    print_header(hdr);
    double *buf = new double[hdr.num_coeff];
    int num_words;
    while ((num_words = fread(buf, sizeof(double), hdr.num_coeff, strm)) ==
                                                            hdr.num_coeff) {
        printf ("%24.16e %24.16e %24.16e\n", buf[0], buf[1], buf[2]);
    }
//    printf ("num_words = %d\n", num_words);
}

