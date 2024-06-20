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

    double x0 = 5.0 * AU;
    double y0 = 5.0 * AU;
    double z0 = -0.05 * AU;
    double size0 = inputs->partsize;
    int np = inputs->nparts;
    double sizes[np];
    double xs[np];
    double ys[np];
    double zs[np];
    for (int i=0; i<np; i++) {
        sizes[i] = size0; // /( (double)pow(10.0,i) );
        xs[i] = x0;
        ys[i] = y0;
        zs[i] = z0 + (0.1*AU/np)*i;
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
        // if ((i%10) == 0) {
        //     sprintf(filename, "%s/diffout%d.txt",inputs->outputdir,i);
        //     printf("%s\n",filename);
        // } else {
        //     strcpy(filename,"NULL");
        // }
        sprintf(filename, "%s/diffout%d.txt",inputs->outputdir,i);
        printf("Starting number: %d\nSaving output to %s\n",i,filename);
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