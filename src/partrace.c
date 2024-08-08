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
    if ( makedir(inputs->outputdir) < 0 ) { exit(1); }
    FILE *fin;
    char inputout[100];
    sprintf(inputout,"%s/inputs.in",inputs->outputdir);
    fin = fopen(inputout,"w");
    if (fin==NULL) {
        fprintf(stderr,"Cannot create input file in output directory!");
        exit(1);
    }
    fprintf_Inputs(fin,inputs);
    fclose(fin);
    
    const size_t nx = NX;
    const size_t ny = NY;
    const size_t nz = NZ;
    
    printf("making Model...\n");
    Model *model = init_Model(inputs->fargodir,inputs->nout,nx,ny,nz);
    printf("Model initialized!\n");

    // seed the random number generator
    srand(time(NULL));

    double rmin = 8*AU;
    double rmax = 10*AU;
    double phimin = -M_PI;
    double phimax = M_PI;
    double zmin = -0.5; //scaleheight
    double zmax =  0.5; //scaleheight
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
        //double r = rmin+ i*(rmax-rmin)/np;
        double z = r*0.05*random_range(zmin,zmax);
        xs[i] = r*cos(phi);
        ys[i] = r*sin(phi);
        zs[i] = z;
    }

    double t0 = inputs->t0;
    double tf = inputs->tf;
    double dtout = inputs->dtout;
    int final_status=0;
    int all_final[np];

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

    char allpartsFilename[100];
    sprintf(allpartsFilename, "%s/allparts.txt",inputs->outputdir);
    FILE *allpartsf;
    allpartsf = fopen(allpartsFilename,"w+");
    fprintf(allpartsf,"x0\ty0\tz0\txf\tyf\tzf\n");
    fclose(allpartsf);

    // TODO: Parallelize this loop
    for (int i=0; i<np; i++) {
        printf("Starting loop\n");
        char filename[100];
        // save every 10th output
        if ((i%1) == 0) {
            sprintf(filename, "%s/particle%d.txt",inputs->outputdir,i);
        } else {
            strcpy(filename,"NULL");
        }
        printf("Starting number: %d\n",i);
        if (strcmp(filename,"NULL") != 0) {
            printf("Saving output to %s\n",filename);
        }
        Particle *p = init_Particle(model, sizes[i], xs[i], ys[i], zs[i]);
        printf("Integrating...\n");
        final_status = integrate(p, t0, tf, dtout, inputs->diffusion,
                                 filename, resFilename, velFilename, crossFilename);
        // save to the allparts file
        allpartsf = fopen(allpartsFilename,"a");
        fprintf(allpartsf, "%f\t%f\t%f\t%f\t%f\t%f\n",xs[i],ys[i],zs[i],p->x,p->y,p->z);
        fclose(allpartsf);
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

int run_partrace(char *inputfile) {
    char *argv[2];
    strcpy(argv[0],"./partrace");
    strcpy(argv[1],inputfile);
    return main(2,argv);
}