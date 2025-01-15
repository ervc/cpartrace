#include <stdlib.h>
#include <stdio.h>

// #include "constants.h"

typedef struct Domain {
    char fargodir[100];
    size_t nx;
    size_t ny;
    size_t nz;

    // domain
    double* phiCenters;
    double* rCenters;
    double* thetaCenters;
    double* phiEdges;
    double* rEdges;
    double *thetaEdges;
    // grids
    // MeshField* phiGrid;
    // MeshField* rGrid;
    // MeshField* thetaGrid;
} Domain;

void read_domfile(double* edges, double* centers, size_t nx, char* domfile, int ghostCells, int rescale);

Domain *init_Domain(char* fargodir, size_t nx, size_t ny, size_t nz) {
    Domain *domain = (Domain*)malloc(sizeof(*domain));
    // string copy directory
    snprintf(domain->fargodir,100,"%s",fargodir);
    domain->nx = nx;
    domain->ny = ny;
    domain->nz = nz;

    char xfile[100];
    char yfile[100];
    char zfile[100];
    int cx;
    cx = snprintf(xfile,100,"%s/domain_x.dat",fargodir);
    cx = snprintf(yfile,100,"%s/domain_y.dat",fargodir);
    cx = snprintf(zfile,100,"%s/domain_z.dat",fargodir);
    if (cx>100) {
        perror("Fargodir is too long!");
        exit(1);
    }

    // allocate memory for domain arrays
    domain->phiCenters = (double*)malloc(sizeof(double)*nx);
    domain->rCenters = (double*)malloc(sizeof(double)*ny);
    domain->thetaCenters = (double*)malloc(sizeof(double)*nz);
    domain->phiEdges = (double*)malloc(sizeof(double)*(nx+1));
    domain->rEdges = (double*)malloc(sizeof(double)*(ny+1));
    domain->thetaEdges = (double*)malloc(sizeof(double)*(nz+1));

    read_domfile(domain->phiEdges,domain->phiCenters,nx,xfile,0,0);
    read_domfile(domain->rEdges,domain->rCenters,ny,yfile,3,1);
    read_domfile(domain->thetaEdges,domain->thetaCenters,nz,zfile,3,0);

    // create_grids(domain);

    return domain;
}

Domain *init_Jupiter_Domain(char* fargodir, size_t nx, size_t ny, size_t nz) {
    Domain *domain = (Domain*)malloc(sizeof(*domain));
    // string copy directory
    snprintf(domain->fargodir,100,"%s",fargodir);
    domain->nx = nx;
    domain->ny = ny;
    domain->nz = nz;

    char xfile[100];
    char yfile[100];
    char zfile[100];
    int cx;
    cx = snprintf(xfile,100,"%s/domain_x.dat",fargodir);
    cx = snprintf(yfile,100,"%s/domain_y.dat",fargodir);
    cx = snprintf(zfile,100,"%s/domain_z.dat",fargodir);
    if (cx>100) {
        perror("Fargodir is too long!");
        exit(1);
    }

    // allocate memory for domain arrays
    domain->phiCenters = (double*)malloc(sizeof(double)*nx);
    domain->rCenters = (double*)malloc(sizeof(double)*ny);
    domain->thetaCenters = (double*)malloc(sizeof(double)*nz);
    domain->phiEdges = (double*)malloc(sizeof(double)*(nx+1));
    domain->rEdges = (double*)malloc(sizeof(double)*(ny+1));
    domain->thetaEdges = (double*)malloc(sizeof(double)*(nz+1));

    // read in the domain files, no ghost cells or rescaling
    read_domfile(domain->phiEdges,domain->phiCenters,nx,xfile,0,0);
    read_domfile(domain->rEdges,domain->rCenters,ny,yfile,0,0);
    read_domfile(domain->thetaEdges,domain->thetaCenters,nz,zfile,0,0);

    // create_grids(domain);

    return domain;
}

void free_Domain(Domain *domain) {
    free (domain->phiCenters);
    free (domain->rCenters);
    free (domain->thetaCenters);
    free (domain->phiEdges);
    free (domain->rEdges);
    free (domain->thetaEdges);
    // free_MeshField(domain->phiGrid);
    // free_MeshField(domain->rGrid);
    // free_MeshField(domain->thetaGrid);
    free (domain);
}

void read_domfile(double *edges, double *centers, size_t nx, char* domfile, int ghostCells, int rescale) {
    /**
     * Reads the domain file EDGES into an array of CENTER values. 
     * ghostCells is the number of ghost cells to ignore when reading in 
     * cell edges. rescale is a bool to rescale the value or not.
     * 
     */
    // printf("Reading domain from %s\n",domfile);
    FILE *file;
    file = fopen(domfile,"r");
    if (file==NULL) {
        printf("Cannot open domfile: %s\n",domfile);
        perror("Error!");
        exit(1);
    }
    // left and right cell edges
    double ledge = 0.0;
    double redge = 0.0; 
    // rescale lengths to cm, dont rescale angles though!
    double scale = 1.0;
    if (rescale) {
        scale = R0;
    }

    for (int i=0; i<nx+ghostCells+1; i++) {
        fscanf(file,"%lf",&redge);
        if (i<=ghostCells) {
            ledge=redge;
            continue;
        }
        double center = (ledge+redge)/2;
        centers[i-ghostCells-1] = center*scale;
        edges[i-ghostCells-1] = ledge*scale;
        ledge = redge;
    }
    edges[nx] = redge*scale;
}

// void create_grids(Domain* domain) {
//     size_t nx = domain->nx;
//     size_t ny = domain->ny;
//     size_t nz = domain->nz;
//     //allocate the memory
//     // domain->phiGrid   = init_MeshField(nx,ny,nz);
//     // domain->rGrid     = init_MeshField(nx,ny,nz);
//     // domain->thetaGrid = init_MeshField(nx,ny,nz);

//     size_t idx = 0;
//     double phi,r,theta;
//     for (size_t k=0; k<nz; k++) {
//         theta = domain->thetaCenters[k];
//         for (size_t j=0; j<ny; j++) {
//             r = domain->rCenters[j];
//             for (size_t i=0; i<nx; i++) {
//                 phi = domain->phiCenters[i];

//                 // domain->phiGrid->data[idx] = phi;
//                 // domain->rGrid->data[idx] = r;
//                 // domain->thetaGrid->data[idx] = theta;

//                 idx++;
//             }
//         }
//     }
// }