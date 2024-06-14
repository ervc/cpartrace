#include <stdlib.h>
#include <stdio.h>

// #include "constants.h"

typedef struct MeshField {
    /**
     * @brief Struct to store and access data in arrays for FARGO outputs
     * 
     */
    size_t nx;
    size_t ny;
    size_t nz;
    double *data;
} MeshField;

// functions
MeshField *init_MeshField(size_t nx, size_t ny, size_t nz);
MeshField *init_MeshField_fromFile(char* fname, size_t nx, size_t ny, size_t nz, int rescale);
size_t get_idx(MeshField *mesh, size_t k, size_t j, size_t i);
void read_datfile(MeshField *mesh, char* fname, int rescale);
double get_data(MeshField *mesh, size_t k, size_t j, size_t i);
void set_data(MeshField *mesh, size_t k, size_t j, size_t i, double value);


// definitions
MeshField *init_MeshField(size_t nx, size_t ny, size_t nz) {
    // initialize an empty MeshField struct
    MeshField *mesh = (MeshField*)malloc(sizeof(*mesh));
    if (!mesh) {
        perror("Malloc Failed");
        exit(1);
    }
    mesh->nx = nx;
    mesh->ny = ny;
    mesh->nz = nz;
    mesh->data = (double*)malloc(sizeof(double)*nx*ny*nz);
    if (!mesh->data) {
        perror("Malloc Failed");
        exit(1);
    }
    return mesh;
}

MeshField *init_MeshField_fromFile(
    char* fname, size_t nx, size_t ny, size_t nz, int rescale) {
    // Initialize meshfield from a file
    MeshField *mesh = init_MeshField(nx,ny,nz);
    read_datfile(mesh,fname,rescale);
    return mesh;
}

void free_MeshField(MeshField *mesh) {
    // Free the memory for the meshfield
    free(mesh->data);
    free(mesh);
}

// If you want to modify the struct, pass a pointer. If you just want to
// read from the struct, pass the object

size_t get_idx(MeshField *mesh, size_t k, size_t j, size_t i) {
    return k*mesh->nx*mesh->ny + j*mesh->nx + i;
}

void set_data(MeshField *mesh, size_t k, size_t j, size_t i, double value) {
    size_t idx = get_idx(mesh,k,j,i);
    mesh->data[idx] = value;
}

double get_data(MeshField *mesh, size_t k, size_t j, size_t i) {
    size_t idx = get_idx(mesh,k,j,i);
    return mesh->data[idx];
}

void read_datfile(MeshField *mesh, char* fname, int rescale) {
    // init a MeshField struct with data from a file
    // printf("Reading from %s\n",fname);
    int nx=mesh->nx;
    int ny=mesh->ny;
    int nz=mesh->nz;
    // read the data in
    // fargo data is stored as sequence of 8 byte doubles
    FILE* file;
    file = fopen(fname,"rb");
    // for scaling the data to cgs
    double scale = 1.0;
    switch (rescale) {
        case RHO:
            scale = MSUN/R0/R0/R0; //g cm-3
            break;
        case VPHI:
        case VR:
        case VTHETA:
            scale = R0/TIME;
            break;
        default:
            scale = 1.0;
            break;
    }
    size_t idx = 0;
    for (int k=0; k<nz; k++) {
        for (int j=0; j<ny; j++) {
            for (int i=0; i<nx; i++) {
                char buffer[8];
                fread (buffer,8,1,file);
                double value = *((double*)buffer);
                // rescale the data here to cgs
                mesh->data[idx] = value*scale;
                idx++;
                // set_data(mesh,k,j,i,value*scale);
                }
            }
        }
    fclose(file);
}

