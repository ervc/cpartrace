#include "partrace.h"
// this technically makes this a c++ file
// #include "mlinterp.hpp"

// defaults
#define NX 2048
#define NY 256
#define NZ 32

int main() {
    if (DIFFUSION) {
        // seed the random number generator
        srand(time(NULL));
    }
    
    const size_t nx = NX;
    const size_t ny = NY;
    const size_t nz = NZ;
    printf("Hello, World\n");
    int nout = 50;
    char fargodir[100] = "/Users/ericvc/fargo/outputs/alpha3_mplan0/";
    printf("making Model...\n");
    Model *model = init_Model(fargodir,nout,nx,ny,nz);
    printf("Model initialized!\n");

    double x0 = 5.0 * AU;
    double y0 = 5.0 * AU;
    double z0 = -0.05 * AU;
    double size0 = 1.0;
    int np = 1;
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

    double tf = 1.0e4*YR;
    double dtout = 1*YR;
    int final_status=0;
    int all_final[np];
    char filename[50];
    if (!DIFFUSION) {
        printf("No ");
    } else {
        printf("Yes ");
    }
    printf("Diffusion\n");
    for (int i=0; i<np; i++) {
        sprintf(filename, "outputs/diffout%d.txt",i);
        printf("Starting number: %d\nSaving output to %s\n",i,filename);
        Particle *p = init_Particle(model, sizes[i], xs[i], ys[i], zs[i]);
        printf("Integrating...\n");
        final_status = integrate(p, tf, dtout, filename);
        all_final[i] = final_status;
        free_Particle(p);
    }

    printf("All statuses: ");
    for (int i=0; i<np; i++) {
        printf("%d, ",all_final[i]);
    }
    printf("\n");

    free_Model(model);
    return 0;
}