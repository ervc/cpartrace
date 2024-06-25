#include "partrace.h"
// this technically makes this a c++ file
// #include "mlinterp.hpp"

// defaults
#define NX 2048
#define NY 256
#define NZ 32

int main(int argc, char **argv) {
    printf("*** CPARTRACE VERSION %s ***\n",VERSION);
    
    char infile[100];
    Inputs *inputs = init_Inputs();
    if (argc <=1 ) {
        printf("No input supplied, using defaults\n");
    } else {
        strcpy(infile, argv[1]);
        printf("Reading input file: %s\n", infile);
        inputs = read_inputs(infile);
    }
    print_Inputs(inputs);
    
    const size_t nx = NX;
    const size_t ny = NY;
    const size_t nz = NZ;
    
    printf("making Model...\n");
    Model *model = init_Model(inputs->fargodir,inputs->nout,nx,ny,nz);
    printf("Model initialized!\n");

    // seed the random number generator
    srand(time(NULL));

    double rmin = 7*AU;
    double rmax = 12*AU;
    double phimin = -M_PI;
    double phimax = M_PI;
    double z0 = 0*AU;
    double size0 = inputs->partsize;
    int np = inputs->nparts;
    double sizes[np];
    double xs[np];
    double ys[np];
    double zs[np];
    for (int i=0; i<np; i++) {
        sizes[i] = size0; // /( (double)pow(10.0,i) );
        double phi = random_range(phimin,phimax);
        double r = random_range(rmin,rmax);
        double z = z0;
        xs[i] = r*cos(phi);
        ys[i] = r*sin(phi);
        zs[i] = z;
    }

    double t0 = inputs->t0;
    double tf = inputs->tf;
    double dtout = 1*YR;
    int final_status=0;
    int all_final[np];

    if ( makedir(inputs->outputdir) < 0 ) { exit(1); }

    char resFilename[100];
    if (inputs->residenceTimes) {
        sprintf(resFilename, "%s/residenceTimes.dat",inputs->outputdir);
    } else {
        strcpy(resFilename,"NULL");
    }

    char velFilename[100];
    if (inputs->velocities) {
        sprintf(velFilename, "%s/velocities.dat",inputs->outputdir);
    } else {
        strcpy(velFilename,"NULL");
    }

    char crossFilename[100];
    if (inputs->crossings) {
        sprintf(crossFilename, "%s/partCrossings.txt",inputs->outputdir);
    } else {
        strcpy(crossFilename,"NULL");
    }

    // TODO: Parallelize this loop
    for (int i=0; i<np; i++) {
        printf("Starting loop\n");
        char filename[100];
        // save every 10th output
        if ((i%1) == 0) {
            sprintf(filename, "%s/diffout%d.txt",inputs->outputdir,i);
        } else {
            strcpy(filename,"NULL");
        }
        printf("Starting number: %d\n",i);
        if (strcmp(filename,"NULL") != 0) {
            printf("Saving output to %s\n",filename);
        }
        Particle *p = init_Particle(model, sizes[i], xs[i], ys[i], zs[i]);
        printf("Integrating...\n");
        final_status = integrate(p, t0, tf, dtout, filename,
                                 resFilename, velFilename, crossFilename);
        all_final[i] = final_status;
        free_Particle(p);
    }

    printf("All statuses: ");
    for (int i=0; i<np; i++) {
        printf("%d, ",all_final[i]);
    }
    printf("\n");

    free_Inputs(inputs);
    free_Model(model);
    return 0;
}