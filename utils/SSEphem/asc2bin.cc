// This file is part of SSEphem.
//
//  SSEphem is a collection of C/C++ source files that provide an interface
//  to the JPL Solar System Ephemeris.
//  Copyright (C) 2009 Associated Universities, Inc. Washington DC, USA
//
//  SSEphem is free software: you can redistribute it and/or modify
//  it under the terms of the GNU Lesser General Public License as published
//  by the Free Software Foundation, either version 3 of the License, or
//  (at your option) any later version.
//
//  SSEphem is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU Lesser General Public License for more details.
//
//  You should have received a copy of the GNU Lesser General Public License
//  along with SSEphem.  If not, see <http://www.gnu.org/licenses/>.

//   ASC2JPL creates a direct access JPL Planetary Ephemeris file from
//   one or more ascii text files.

//   This program, 'asc2jpl', requires (via standard input) an ascii
//   header file ('header.XXX'), followed by one or more ascii ephemeris 
//   data files ('ascSYYYY.XXX').  All files must have the same ephemeris
//   number, XXX.  Further, the data files must be consecutive in time
//   with no gaps between them. 

//   By default, the output ephemeris will span the same interval as the input
//   text file(s).  If you are interested in only a portion of data, set the
//   below T1 and T2 to the begin and end times of the span you desire.  T1
//   and T2 must be specified in  Julian Ephemeris Days (ET).

//   A sample commnd and sequence of files might be:

//    cat header.403 asc+1980.403 asc+2000.403 | asc2jpl 2450000 2460000

//   (The data files for DE200 and DE403 contain 20 years each;
//   for DE404, 50 years)

//  By default, the output ephemeris will span the same interval as the
//  input ascii data file(s).  The user may reset these to other JED's.

#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <stdlib.h>
#include "jpl_eph.h"

void fix_double_notation(char *str);
int  next_group();
void print_header(de_headerT hdr);

FILE *inst;

int main(int argc, char *argv[])
{
    de_headerT hdr;

    if (argc < 2) {
	printf ("Usage: asc2bin ascii_file [JD beginning [JD end]]\n");
	return 1;
    }
    inst = fopen(argv[1], "rt");
    if (inst == NULL) {
        printf ("Cannot open file %s\n", argv[1]);
        return 1;
    }
    double t1 = 0.0;
    double t2 = 9999999.0;
    if (argc > 2) {
        t1 = atof(argv[2]);
    }
    if (argc > 3) {
        t2 = atof(argv[3]);
    }

    // Write a fingerprint to the screen.
    printf ("\nC code version of the JPL Solar System Ephemerides\n");
    printf ("   ASCII-to-Binary program last modified 24-Dec-1995\n");

    char line_str[85];
    if (fgets(line_str, 84, inst) == NULL) {
        printf ("Unexpected EOF reading number of coefficients\n");
        return 1;
    }
    char label1[20], label2[20];
    int  ksize, ncoeff;
    sscanf(line_str, "%s %d %s %d", label1, &ksize, label2, &ncoeff);
    if (strcmp(label1, "KSIZE=")) {
        printf ("KSIZE label not recognized (%s)\n", line_str);
        return 1;
    }
    if (strcmp(label2, "NCOEFF=")) {
        printf ("NCOEFF label not recognized (%s)\n", line_str);
        return 1;
    }
    if (ncoeff > 0) {
        hdr.num_coeff = ncoeff;
    } else {
        printf ("Number of coefficients too small (%d)\n", ncoeff);
        return 1;
    }

    // Now for the alphameric heading records (GROUP 1010)
    if (next_group() != 1010) {
        printf ("Group 1010 not found\n");
        return 1;
    }

    // get the three title lines
    if (fgets(line_str, 84, inst) == NULL) {
        printf ("Unexpected EOF reading title line 1\n");
        return 1;
    }
    strcpy(hdr.title1, line_str);
    if (fgets(line_str, 84, inst) == NULL) {
        printf ("Unexpected EOF reading title line 2\n");
        return 1;
    }
    strcpy(hdr.title2, line_str);
     if (fgets(line_str, 84, inst) == NULL) {
        printf ("Unexpected EOF reading title line 3\n");
        return 1;
    }
    strcpy(hdr.title3, line_str);

    // Read start, end and record span  (GROUP 1030)
    if (next_group() != 1030) {
        printf ("Group 1030 not found\n");
        return 1;
    }
     if (fgets(line_str, 84, inst) == NULL) {
        printf ("Unexpected EOF reading group 1030\n");
        return 1;
    }
    fix_double_notation(line_str);
    sscanf(line_str, "%lf %lf %lf", &hdr.epoch1, &hdr.epoch2, &hdr.interval);

    // Read number of constants and names of constants (GROUP 1040/4).
    if (next_group() != 1040) {
        printf ("Group 1040 not found\n");
        return 1;
    }
    if (fgets(line_str, 84, inst) == NULL) {
        printf ("Unexpected EOF reading number of constants\n");
        return 1;
    }
    sscanf(line_str, "%d", &hdr.num_constants);
    int i;
    for (i = 0; i < hdr.num_constants; i += 10) {
        if (fgets(line_str, 84, inst) == NULL) {
            printf ("Unexpected EOF reading constants names\n");
            return 1;
        }
        sscanf(line_str, "%s %s %s %s %s %s %s %s %s %s",
               hdr.constant_name[i],     hdr.constant_name[i + 1],
               hdr.constant_name[i + 2], hdr.constant_name[i + 3],
               hdr.constant_name[i + 4], hdr.constant_name[i + 5],
               hdr.constant_name[i + 6], hdr.constant_name[i + 7],
               hdr.constant_name[i + 8], hdr.constant_name[i + 9]);
    }
    // Read number of values and values (GROUP 1041/4)
    if (next_group() != 1041) {
        printf ("Group 1041 not found\n");
        return 1;
    }
    if (fgets(line_str, 84, inst) == NULL) {
        printf ("Unexpected EOF reading number of constant values\n");
        return 1;
    }
    int num_const_values;
    sscanf(line_str, "%d", &num_const_values);
    if (num_const_values != hdr.num_constants) {
        printf ("Number of const labels (%d) and const values (%d) disagree\n",
                 num_const_values, hdr.num_constants);
    }
    for (i = 0; i < num_const_values; i += 3) {
        if (fgets(line_str, 84, inst) == NULL) {
            printf ("Unexpected EOF reading constants values\n");
            return 1;
        }
        fix_double_notation(line_str);
        sscanf(line_str, "%lf %lf %lf", &hdr.constant_value[i],
                     &hdr.constant_value[i + 1],  &hdr.constant_value[i + 2]);
    }
    int found = 0;
    for (i = 0; i < hdr.num_constants; i++) {
        if (!strcmp(hdr.constant_name[i], "AU")) {
            found++;
            hdr.au = hdr.constant_value[i];
            break;
        }
    }
    if (found == 0) {
        printf ("AU value not found\n");
    }
    found = 0;
    for (i = 0; i < hdr.num_constants; i++) {
        if (!strcmp(hdr.constant_name[i], "EMRAT")) {
            found++;
            hdr.emrat = hdr.constant_value[i];
            break;
        }
    }
    if (found == 0) {
        printf ("EMRAT value not found\n");
    }
    found = 0;
    for (i = 0; i < hdr.num_constants; i++) {
        if (!strcmp(hdr.constant_name[i], "DENUM")) {
            found++;
            hdr.DEnumber = (int)hdr.constant_value[i];
            break;
        }
    }
    if (found == 0) {
        printf ("DENUM value not found\n");
    }
    // Read pointers needed by INTERP (GROUP 1050)
    if (next_group() != 1050) {
        printf ("Group 1050 not found\n");
        return 1;
    }
    int *offset = new int[numObj];
    int *numcof = new int[numObj];
    int *numsub = new int[numObj];
    if (fgets(line_str, 84, inst) == NULL) {
        printf ("Unexpected EOF reading group 1050\n");
        return 1;
    }
    sscanf(line_str, "%d %d %d %d %d %d %d %d %d %d %d %d %d",
               &offset[0],  &offset[1],  &offset[2], &offset[3], &offset[4],
               &offset[5],  &offset[6],  &offset[7], &offset[8], &offset[9],
               &offset[10], &offset[11], &offset[12]);
    // reduce these offsets by 1 to conform to C indexing
    int k;
    for (k = 0; k < numObj; k++) {
        offset[k]--;
    }
    if (fgets(line_str, 84, inst) == NULL) {
        printf ("Unexpected EOF reading group 1050\n");
        return 1;
    }
    sscanf(line_str, "%d %d %d %d %d %d %d %d %d %d %d %d %d",
               &numcof[0],  &numcof[1],  &numcof[2], &numcof[3], &numcof[4],
               &numcof[5],  &numcof[6],  &numcof[7], &numcof[8], &numcof[9],
               &numcof[10], &numcof[11], &numcof[12]);
    if (fgets(line_str, 84, inst) == NULL) {
        printf ("Unexpected EOF reading group 1050\n");
        return 1;
    }
    sscanf(line_str, "%d %d %d %d %d %d %d %d %d %d %d %d %d",
               &numsub[0],  &numsub[1],  &numsub[2], &numsub[3], &numsub[4],
               &numsub[5],  &numsub[6],  &numsub[7], &numsub[8], &numsub[9],
               &numsub[10], &numsub[11], &numsub[12]);
    for (i = 0; i < numObj; i++) {
        hdr.obj[i].offset = offset[i];
        hdr.obj[i].num_coeff = numcof[i];
        hdr.obj[i].num_subinterv = numsub[i];
        hdr.obj[i].num_components = 3;
    }
    hdr.obj[Nutation].num_components = 2;
    // Open output file ('DEcXXX')
    char file_name[10];
    sprintf(file_name, "DEc%d", hdr.DEnumber);
    FILE *out_stream = fopen(file_name, "wb");
    if (out_stream == NULL) {
        printf ("Could not open output file %s\n", file_name);
        return 1;
    }
    //  print_header(hdr);
    fwrite(&hdr, sizeof(hdr), 1, out_stream);

    if (next_group() != 1070) {
        printf ("Group 1070 not found\n");
        return 1;
    }

    // The header file has been exhausted.  The following read is from the
    //  first data file.
    if (fgets(line_str, 84, inst) == NULL) {
        printf ("Unexpected EOF reading coefficients file\n");
        return 1;
    }
    int record_num, num_coef;
    sscanf(line_str, "%d %d", &record_num, &num_coef);
    if (num_coef != hdr.num_coeff) {
        printf (
           "Number of coeff (%d) in record %d disagrees with header (%d)\n",
           num_coef, record_num, hdr.num_coeff);
        fclose(out_stream);
        return 1;
    }
    int in_rn = 0;
    int out_rn = 0;
    double last_time = 1.0e50;
    double *ptr = new double[num_coef + 3];
    double *buffer = ptr;
    buffer[1] = 0.0;
    while (buffer[1] < t2) {
        int break_flag = 0;
        for (i = 0; i < num_coef + 2; i += 3) {
            if (fgets(line_str, 84, inst) == NULL) {
                printf ("EOF reading coefficients, after record %d, line %d\n",
                         in_rn, i);
                break_flag++;
                break;
            }
            fix_double_notation(line_str);
            sscanf(line_str, "%lf %lf %lf", ptr, ptr + 1, ptr + 2);
            ptr += 3;
        }
        if (break_flag) break;
        in_rn++;
        ptr = buffer;
        if ((buffer[0] > last_time) && (buffer[1] > t1)) {
            printf ("Coefficients time interval gap, record %d\n", in_rn);
            break;
        }
        if (buffer[1] < t1) {
            if (((in_rn - 1) % 10) == 0) {
                printf (
                   "Searching for requested start epoch, input record %d\n",
                     in_rn);
            }
        }
        last_time = buffer[1];
        if ((buffer[1] >= t1) && (buffer[0] <= t2)) {
            out_rn++;
            if (((out_rn - 1) % 10) == 0) {
                printf ("Converting record (in %d out %d), %5.1f %5.1f\n",
                                        in_rn, out_rn, buffer[0], buffer[1]);
            }
            if (out_rn == 1) {
                hdr.epoch1 = buffer[0];
                hdr.interval = buffer[1] - buffer[0];
            }
            hdr.epoch2 = buffer[1];
            fwrite(buffer, sizeof(double), num_coef, out_stream);
        }
        if (fgets(line_str, 84, inst) == NULL) {
            printf ("Coefficients file(s) finished.\n");
            break;
        }
        sscanf(line_str, "%d %d", &record_num, &num_coef);
        if (num_coef != hdr.num_coeff) {
            printf (
           "Number of coeff (%d) in record %d disagrees with header (%d)\n",
           num_coef, in_rn, hdr.num_coeff);
            break;
        }
    }
    // Update header and close
    rewind(out_stream);
    fwrite(&hdr, sizeof(hdr), 1, out_stream);
    fclose(out_stream);
    return 0;
}

//  Read the GROUP header, discard the blank line after the header, and
//  return the group number.

int next_group()
{
    char line_str[85];
    int char_flag = 0;
    while (char_flag == 0) {
        if (fgets(line_str, 84, inst) == NULL) {
            return -1;
        }
	int i;
        for (i = strlen(line_str) - 1; i >= 0; i--) {
            if ((line_str[i] == ' ') || (line_str[i] == '\n')) {
                line_str[i] = '\0';
            } else {
                char_flag++;
                break;
            }
        }
    }
    char group[85];
    int group_num;
    sscanf(line_str, "%s %d", group, &group_num);
    if (strcmp(group, "GROUP")) {
        printf("Expected group line not seen (%s)\n", line_str);
        return -1;
    }
//    Found the header.  Read the blank line so we can get at the data.
    fgets(line_str, 84, inst);
    return group_num;
}
// Change the FORTRAN double precision notation to C notation.

void fix_double_notation(char *str)
{
    int i;
    for (i = 1; i < strlen(str) - 1; i++) {
        if ((str[i] == 'd') || (str[i] == 'D')) {
            if (isdigit(str[i - 1]) && (isdigit(str[i + 1]) ||
                               (str[i + 1] == '-') || (str[i + 1] == '+'))) {
                str[i] = 'e';
            }
        }
    }
}

char *delete_trailing_blanks(char *str)
{
    int i;
    for (i = strlen(str) - 1; i >= 0; i--) {
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
    printf ("Object offset   \n");
    int i;
    for (i = 0; i < 13; i++) {
        printf ("%5d", hdr.obj[i].offset);
    }
    printf ("No. subintervals\n");
    for (i = 0; i < 13; i++) {
        printf ("%5d", hdr.obj[i].num_subinterv);
    }
    printf ("No. coefficients\n");
    for (i = 0; i < 13; i++) {
        printf ("%5d", hdr.obj[i].num_coeff);
    }
    printf ("\n");
    printf ("DE %d\n", hdr.DEnumber);
    printf ("Constants 1 & 2:     %24.16e  %24.16e\n", hdr.constant_value[0],
                                       hdr.constant_value[1]);
    printf ("Constants 199 & 200: %24.16e  %24.16e\n", hdr.constant_value[198],
                                       hdr.constant_value[199]);
    printf ("%5d coefficients\n", hdr.num_coeff);
}
