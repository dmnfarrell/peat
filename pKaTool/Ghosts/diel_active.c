/*
# pKaTool - scripts for analysing chemical shift perturbations
# Copyright (C) 2010 Predrag Kukic & Jens Erik Nielsen
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
# Contact information:
# Email: Jens.Nielsen@gmail.com
# Normal mail:
# Jens Nielsen
# SBBS, Conway Institute
# University College Dublin
# Dublin 4, Ireland
#
*/

#include "apbscfg.h"
#include "apbs/apbs.h"
#include "stdio.h"

#define IJK(i,j,k)  (((k)*(nx)*(ny))+((j)*(nx))+(i))
#define ERRRC 2



int diel_active(Vgrid *grid, double xcoor, double ycoor, double zcoor, double radius, double diel);

int usage(int rc) {

    char *usage = "\n\n\
    ----------------------------------------------------------------------\n\
    This driver program convolves grid data with various filters.  It is\n\
    invoked as:\n\
      smooth <args>\n\n\
    where <args> is a list of the following arguments:\n\
      REQUIRED GENERAL ARGUMENTS:\n\
      --format=<format>  where <format> specifies the data format and is one\n\
                          of the following: dx (OpenDX)\n\
      --input=<file>     where <file> is the input mesh data in the\n\
                         specified format\n\
      --output=<file>    where <file> is the output mesh data in the\n\
                         specified format\n\
      --diel=<diel>  	where <diel> is the new value used in the \n\
	  					protein active site\n\
      --xcoor=<xcoor>  	where <xcoor> is the xcoordinate of the \n\
	  					point around which you define the protein active site\n\
      --ycoor=<ycoor>  	where <ycoor> is the xcoordinate of the \n\
	  					point around which you define the protein active site\n\
      --zcoor=<zcoor>  	where <zcoor> is the xcoordinate of the \n\
	  					point around which you define the protein active site\n\
      --radius=<radius>  where <radius> is the radius used in specifying the \n\
	  					protein active site\n\
    ----------------------------------------------------------------------\n\n";

    Vnm_print(2, usage);

    exit(rc);
 
    return 0;
}

int main(int argc, char **argv) {

    /* *************** VARIABLES ******************* */
    int i;
    Vgrid *grid = VNULL;
    /* Input parameters */
    Vdata_Format format; int gotFormat = 0;
    char inPath[VMAX_BUFSIZE]; int gotInPath = 0;
    char outPath[VMAX_BUFSIZE]; int gotOutPath = 0;
	double diel; int gotDiel = 0;
    double xcoor; int gotXcoor = 0;
    double ycoor; int gotYcoor = 0;
	double zcoor; int gotZcoor = 0;
	double radius; int gotRadius = 0;
	

    char *header = "\n\n\
    ----------------------------------------------------------------------\n\
    ----------------------------------------------------------------------\n\
    \n\n";

    /* *************** CHECK INVOCATION ******************* */
    Vio_start();
    Vnm_redirect(1);
    Vnm_print(1, "%s", header);
    for (i=1; i<argc; i++) {
        Vnm_print(1, "Parsing: %s...\n", argv[i]);
        if (strstr(argv[i], "--format") != NULL) {
            if (strstr(argv[i], "dx") != NULL) {
                gotFormat = 1;
                format = VDF_DX;
            } else {
                Vnm_print(2, "Error:  %s\n", argv[i]);
                usage(2);
            }
        } else if (strstr(argv[i], "--input") != NULL) {
            if (sscanf(argv[i], "--input=%s", inPath) == 1) gotInPath = 1;
            else {
                Vnm_print(2, "Error:  %s\n", argv[i]);
                usage(2);
            }
        } else if (strstr(argv[i], "--output") != NULL) {
            if (sscanf(argv[i], "--output=%s", outPath) == 1) gotOutPath = 1;
            else {
                Vnm_print(2, "Error:  %s\n", argv[i]);
                usage(2);
            }
		} else if (strstr(argv[i], "--diel") != NULL) {
            if (sscanf(argv[i], "--diel=%lf", &diel) == 1) gotDiel = 1;
            else {
                Vnm_print(2, "Error:  %s\n", argv[i]);
                usage(2);
            }
        } else if (strstr(argv[i], "--xcoor") != NULL) {
            if (sscanf(argv[i], "--xcoor=%lf", &xcoor) == 1) gotXcoor = 1;
            else {
                Vnm_print(2, "Error:  %s\n", argv[i]);
                usage(2);
            }
		} else if (strstr(argv[i], "--ycoor") != NULL) {
            if (sscanf(argv[i], "--ycoor=%lf", &ycoor) == 1) gotYcoor = 1;
            else {
                Vnm_print(2, "Error:  %s\n", argv[i]);
                usage(2);
            }
		} else if (strstr(argv[i], "--zcoor") != NULL) {
            if (sscanf(argv[i], "--zcoor=%lf", &zcoor) == 1) gotZcoor = 1;
            else {
                Vnm_print(2, "Error:  %s\n", argv[i]);
                usage(2);
            }
        } else if (strstr(argv[i], "--radius") != NULL) {
            if (sscanf(argv[i], "--radius=%lf", &radius) == 1) 
                gotRadius = 1;
            else {
                Vnm_print(2, "Error:  %s\n", argv[i]);
                usage(2);
            }
        } else {
            Vnm_print(2, "Error:  %s\n", argv[i]);
            usage(2);
        }
    }
    if (!gotFormat) {
        Vnm_print(2, "Error:  --format not specified!\n");
        usage(2);
    } 
    if (!gotInPath) {
        Vnm_print(2, "Error:  --input not specified!\n");
        usage(2);
    }
    if (!gotOutPath) {
        Vnm_print(2, "Error:  --output not specified!\n");
        usage(2);
    }
	if (!gotDiel) {
        Vnm_print(2, "Error:  --diel new dielectric value not specified!\n");
        usage(2);
    }
    if (!gotRadius) {
        Vnm_print(2, "Error:  --radius not specified!\n");
        usage(2);
    }
	if (!(gotXcoor && gotYcoor && gotZcoor)) {
        Vnm_print(2, "Error:  --xcoor --ycoor --zcoor not specified!\n");
        usage(2);
    }

    /* *************** READ DATA ******************* */
    Vnm_print(1, "main:  Reading data from %s...\n", inPath);
    grid = Vgrid_ctor(0, 0, 0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, VNULL);
    if (format == VDF_DX) {
        if (!Vgrid_readDX(grid, "FILE", "ASC", VNULL, inPath)) {
            Vnm_print(2, "main:  Problem reading OpenDX-format grid from %s\n",
              inPath);
            return ERRRC;
        }
    }

    /* *************** SMOOTH ******************* */
	Vnm_print(1, "Changing the dielectric constant of the active site...\n");
	diel_active(grid, xcoor, ycoor, zcoor, radius, diel);

    /* *************** READ DATA ******************* */
    Vnm_print(1, "main:  Writing data to %s...\n", outPath);
    if (format == VDF_DX) 
      Vgrid_writeDX(grid, "FILE", "ASC", VNULL, outPath, "Changed active site dielectric constant",
        VNULL);

    return 0;
}

int diel_active(Vgrid *grid, double xcoor, double ycoor, double zcoor, double radius, double diel) {

    int nx, ny, nz, iband, jband, kband, i, j, k, ii, jj, kk;
    int kkmin, jjmin, iimin, kkmax, jjmax, iimax;
    double hx, hy, hzed, xmin, ymin, zmin, *newData, *oldData;
    double dist2, radius2;

    Vnm_print(1, "Changing the diel constant around the point: xcoord. = %g A, ycoord = %g A, zcoord = %g A using radius=%g A.\n",
      xcoor, ycoor, zcoor, radius);

    nx = grid->nx; ny = grid->ny; nz = grid->nz;
    hx = grid->hx; hy = grid->hy; hzed = grid->hzed;
    xmin = grid->xmin; ymin = grid->ymin; zmin = grid->zmin;
    Vnm_print(1, "Grid:  %d x %d x %d points\n", nx, ny, nz);
    Vnm_print(1, "Grid:  %g, %g, %g A spacing\n", hx, hy, hzed);
    Vnm_print(1, "Grid:  (%g, %g, %g) A origin\n", xmin, ymin, zmin);

    /* Convert HALF bandwidth to grid units */
    iband = (int)(radius/hx);
    jband = (int)(radius/hy);
    kband = (int)(radius/hzed);
    Vnm_print(1, "The region converted to %d x %d x %d grid units.\n",
      iband, jband, kband);
    Vnm_print(1, "This means water dielectric constant within (%g, %g, %g)\n", 
      (iband+1)*hx, (jband+1)*hy, (kband+1)*hzed);
    Vnm_print(1, "around the given point will be given the new dielctric constant value.\n");

    /* Get data */
    oldData = grid->data;
    newData = Vmem_malloc(VNULL, (nx*ny*nz), sizeof(double));
    radius2 = VSQR(radius);

	/* Copy over old data.  All but the boundary values will be replaced in the next step so this is more copying than is strictly
	   necessary... */
	for (i=0; i<(nx*ny*nz); i++) newData[i] = oldData[i];

    /* Apply filter */
    for (k=1; k<(nz-1); k++) {
        kkmin = VMAX2(1, (k - kband));
        kkmax = VMIN2((nz-1), (k + kband));
        for (j=1; j<(ny-1); j++) {
            jjmin = VMAX2(1, (j - jband));
            jjmax = VMIN2((ny-1), (j + jband));
            for (i=1; i<(nx-1); i++) {
                iimin = VMAX2(1, (i - iband));
                iimax = VMIN2((nx-1), (i + iband));
                dist2 = VSQR(xmin+hx*i-xcoor) + VSQR(ymin+hy*j-ycoor) + VSQR(zmin+hzed*k-zcoor);
                if (dist2 < radius2) {
                    if (newData[IJK(i,j,k)] > 78.0) newData[IJK(i,j,k)] = diel;
				}
            } /* i loop */
        } /* j loop */
    } /* k loop */

    /* Replace data */
    for (i=0; i<(nx*ny*nz); i++) grid->data[i] = newData[i];
    Vmem_free(VNULL, (nx*ny*nz), sizeof(double), (void **)&newData);

    return 0;
}
