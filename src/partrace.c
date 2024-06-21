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

    double r0 = 8*AU;
    double phi0 = -M_PI;
    double z0 = 0*AU;
    double size0 = inputs->partsize;
    int np = inputs->nparts;
    double dphi = 2*M_PI/np;
    double sizes[np];
    double xs[np];
    double ys[np];
    double zs[np];
    for (int i=0; i<np; i++) {
        sizes[i] = size0; // /( (double)pow(10.0,i) );
        double phi = phi0 + dphi*i;
        xs[i] = r0*cos(phi);
        ys[i] = r0*sin(phi);
        zs[i] = z0;
    }

    double t0 = inputs->t0;
    double tf = inputs->tf;
    double dtout = 1*YR;
    int final_status=0;
    int all_final[np];
    

    if (DIFFUSION) {
        printf("Seeding rng for diffusion\n");
        // seed the random number generator
        srand(time(NULL));
    }

    if ( makedir(inputs->outputdir) < 0 ) { exit(1); }

    for (int i=0; i<np; i++) {
        printf("Starting loop\n");
        char filename[100];
        // save every 10th outputß
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
        final_status = integrate(p, t0, tf, dtout, filename);
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