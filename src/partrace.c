#include "partrace.h"
#include <mpi.h>
// this technically makes this a c++ file
// #include "mlinterp.hpp"

// defaults
// my fargo models:
// NX	2048 or 1024
// NY	256 or 128
// NZ	32 or 36
// Felipe's radmc models:
// NX   512
// NY   256
// NZ   32
#define NX 2048
#define NY 256
#define NZ 32
#define NLVL 5

void init_random_particles(Inputs *inputs, double *sizes, double *xs, double *ys, double *zs);
void read_partfile(Inputs *inputs, double *sizes, double *xs, double *ys, double *zs);

int main(int argc, char **argv) {
    int rank = 0, nprocs = 1;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

    if (rank == 0) {
        printf("*** CPARTRACE VERSION %s ***\n", VERSION);
        printf("Parallel implementation. Running on %d procs\n", nprocs);
    }
    printf("  Rank %d / %d running\n", rank, nprocs);

    // read inputs
    char infile[100];
    Inputs *inputs = init_Inputs();
    if (argc <= 1) {
        if (rank == 0) printf("No input supplied, using defaults\n");
    } else {
        strcpy(infile, argv[1]);
        if (rank == 0) printf("Reading input file: %s\n", infile);
        inputs = read_inputs(infile);
    }
    if (rank == 0) {
        if (makedir(inputs->outputdir) < 0) { 
            MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
        }
    }
    MPI_Barrier(MPI_COMM_WORLD);

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

    // make the model
    if (rank==0) printf("making Model...\n");
    int nlvl = 1;
    if (inputs->modeltype==JUPITER_MODEL) {nlvl = 5;}
    else {nlvl = 1;}
    Model *models[nlvl];
    if (inputs->modeltype==JUPITER_MODEL) {
        // size_t nxs[] = {680, 120, 120, 120, 120};
        // size_t nys[] = {215, 120, 120, 120, 120};
        // size_t nzs[] = {20, 34, 62, 86, 86};
        for (int i=0; i<NLVL; i++) {
            char leveldir[100];
            sprintf(leveldir,"%s/fargolev%d/",inputs->fargodir,i);
            models[i] = init_Model(inputs->modeltype,leveldir,"0", rank);
        }
        // default model level0
        if (rank == 0) printf("Models initialized\n");
    } else {
        nlvl = 1;
        models[0] = init_Model(inputs->modeltype,inputs->fargodir,inputs->nout, rank);
        if (rank == 0) printf("Model initialized!\n");
    }

    // seed the random number generator
    // make sure seed is different on different ranks!
    srand(time(NULL) * (rank+1));

    int np = inputs->nparts;
    double *sizes = malloc(np * sizeof(double));
    double *xs = malloc(np * sizeof(double));
    double *ys = malloc(np * sizeof(double));
    double *zs = malloc(np * sizeof(double));
    if (sizes == NULL || xs == NULL || ys == NULL || zs == NULL) {
        fprintf(stderr, "Rank %d: failed to allocate particle arrays\n", rank);
        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    }

    if (rank == 0) {
        if (strcmp(inputs->partfile, "NULL") == 0) {
            init_random_particles(inputs, sizes, xs, ys, zs);
        } else {
            read_partfile(inputs, sizes, xs, ys, zs);
        }
    }

    MPI_Bcast(sizes, np, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(xs, np, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(ys, np, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(zs, np, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    double t0 = inputs->t0;
    double tf = inputs->tf;
    double dtout = inputs->dtout;
    int BACKWARDS = 0; // is the integration backwards in time?
    if ( tf < t0 ) {
        printf("Integrating backwards in time\n");
        BACKWARDS = 1;
        if (inputs->diffusion) {
            printf("!!! ERROR: Integrating backwards with diffusion is not physically correct and is not recommended !!!\n");
            return 1;
        }
    }
    if ( ((tf-t0)*dtout)<0 ) {
        if (BACKWARDS) {
            printf("Integration is backwards in time but dtout is postive. Setting dtout = -dtout\n");
        }
        else {
            printf("Integration is forward in time but dtout is negative. Seting dtout = -dtout\n");
        }
        dtout = -dtout;
    }
    Intout result;
    result.status = 0;
    result.tf = 0.0;
    int *all_final = malloc(np * sizeof(int));
    if (all_final == NULL) {
        fprintf(stderr, "Rank %d: failed to allocate all_final\n", rank);
        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    }

    char resFilename[128];
    if (inputs->residenceTimes) {
        if (nlvl > 1) {
            if (rank == 0) printf("Cannot currently track residence times with multilevel model\n");
            MPI_Finalize();
            return 1;
        }
        if (BACKWARDS) {
            if (rank == 0) printf("Cannot track residence times with backwards integration\n");
            MPI_Finalize();
            return 1;
        }
        Model* model = models[0];
        if (nprocs > 1) {
            sprintf(resFilename, "%s/residenceTimes_rank%d.dat", inputs->outputdir, rank);
        } else {
            sprintf(resFilename, "%s/residenceTimes.dat", inputs->outputdir);
        }
        if (inputs->reset) {
            if (rank == 0) {
                printf("!!! Resetting Residence Times !!!\n");
            }
            FILE *resFile;
            resFile = fopen(resFilename,"wb");
            size_t bigSize = 2*model->nz*model->ny*model->nx;
            double zero = 0;
            // bigSize+1 because we also include the number of particles as the first double
            for (int i=0; i<bigSize+1; i++) {
                fwrite(&zero, sizeof(double), 1, resFile);
            }
            fclose(resFile);
        }
    } else {
        strcpy(resFilename,"NULL");
    }

    char velFilename[128];
    if (inputs->velocities) {
        if (nprocs > 1) {
            sprintf(velFilename, "%s/velocities_rank%d.dat", inputs->outputdir, rank);
        } else {
            sprintf(velFilename, "%s/velocities.dat", inputs->outputdir);
        }
    } else {
        strcpy(velFilename, "NULL");
    }

    char crossFilename[128];
    if (inputs->crossings) {
        if (nprocs > 1) {
            sprintf(crossFilename, "%s/partCrossings_rank%d.txt", inputs->outputdir, rank);
        } else {
            sprintf(crossFilename, "%s/partCrossings.txt", inputs->outputdir);
        }
        // if reset option then write a blank file
        if (inputs->reset) {
            FILE *crossFile;
            crossFile = fopen(crossFilename, "w");
            fclose(crossFile);
        }
    } else {
        strcpy(crossFilename, "NULL");
    }

    char allpartsFilename[128];
    if (nprocs > 1) {
        sprintf(allpartsFilename, "%s/allparts_rank%d.txt", inputs->outputdir, rank);
    } else {
        sprintf(allpartsFilename, "%s/allparts.txt", inputs->outputdir);
    }
    FILE *allpartsf;
    // if the file doesn't exist yet or reset is picked, create it and write the header
    if(!fileExists(allpartsFilename) || inputs->reset) {
        allpartsf = fopen(allpartsFilename,"w");
        fprintf(allpartsf,"tf\tx0\ty0\tz0\txf\tyf\tzf\tstatus\n");
        fclose(allpartsf);
    }

    // parallel loop over particles
    for (int i = rank; i < np; i += nprocs) {
        printf("[rank %d] Starting particle %d\n", rank, i);
        char filename[128];
        // save every dsave-th output
        if ((i % inputs->dsave) == 0) {
            sprintf(filename, "%s/particle%d.txt", inputs->outputdir, i + inputs->nstart);
        } else {
            strcpy(filename, "NULL");
        }
        if (strcmp(filename, "NULL") != 0) {
            printf("[rank %d] Saving output to %s\n", rank, filename);
        }
        Particle *p = init_Particle(models, nlvl, sizes[i], xs[i], ys[i], zs[i]);
        printf("[rank %d] Integrating...\n", rank);
        result = integrate(p, t0, tf, dtout, inputs->diffusion,
                                 filename, resFilename, velFilename, crossFilename);
        // save to the allparts file
        allpartsf = fopen(allpartsFilename, "a");
        fprintf(allpartsf, "%f\t%f\t%f\t%f\t%f\t%f\t%f\t%d\n", result.tf, xs[i], ys[i], zs[i], p->x, p->y, p->z, result.status);
        fclose(allpartsf);
        all_final[i] = result.status;
        free_Particle(p);
    }

    printf("[rank %d] statuses: ", rank);
    for (int i = rank; i < np; i += nprocs) {
        printf("%d, ", all_final[i]);
    }
    printf("\n");

    free_Inputs(inputs);
    free_Models(models, nlvl);
    free(sizes);
    free(xs);
    free(ys);
    free(zs);
    free(all_final);

    MPI_Finalize();
    return 0;
}

void init_random_particles(Inputs *inputs, double *sizes, double *xs, double *ys, double *zs) {
    double rmin = inputs->rmin;
    double rmax = inputs->rmax;
    double phimin = inputs->phimin;
    double phimax = inputs->phimax;
    double thetamin = inputs->thetamin;
    double thetamax = inputs->thetamax;
    for (int i=0; i<inputs->nparts; i++) {
        sizes[i] = inputs->partsize;
        double phi = random_range(phimin,phimax);
        double r = random_range(rmin,rmax);
        double theta = random_range(thetamin,thetamax);
        xs[i] = r*cos(phi)*sin(theta);
        ys[i] = r*sin(phi)*sin(theta);
        zs[i] = r*cos(theta);
    }
}

void read_partfile(Inputs *inputs, double *sizes, double *xs, double *ys, double *zs) {
    FILE *file;
    file = fopen(inputs->partfile, "r");
    if (file==NULL) {
        printf("Cannot open partfile: %s\n", inputs->partfile);
        exit(EXIT_FAILURE);
    }
    double s, x, y, z;
    int nline = 0;
    while ( fscanf(file, "%lf %lf %lf %lf", &s, &x, &y, &z) == 4 ) {
        if ( nline>inputs->nparts ) {
            fprintf(stdout, "TOO MANY LINES IN PARTFILE!\n");
            exit(EXIT_FAILURE);
        }
        sizes[nline] = s;
        xs[nline] = x;
        ys[nline] = y;
        zs[nline] = z;
        nline++;
    }
    fclose(file);
}