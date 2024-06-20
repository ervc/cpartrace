/**
 * We need to read in:
 * strings:
 *  fargodir, outputdir, nout, partfile
 * doubles:
 *  t0, tf, partsize, partdens
 * integers:
 *  nparts
 * booleans:
 *  diffusion, residence times, velocities
 */

#include <errno.h>
#include <sys/stat.h>

typedef struct Inputs {
    char *fargodir;
    char *outputdir;
    char *nout;
    char *partfile;

    double t0;
    double tf;
    double partsize;
    double partdens;

    int nparts;

    int diffusion;
    int residenceTimes;
    int velocities;
} Inputs;

Inputs *init_Inputs() {
    Inputs *in = malloc(sizeof(*in));
    // assign defaults
    in->fargodir  = malloc(sizeof(char)*100);
    strcpy(in->fargodir,"/Users/ericvc/fargo/outputs/alpha3_mplan100/");
    in->outputdir = malloc(sizeof(char)*100);
    strcpy(in->outputdir,"outputs/");
    in->nout      = malloc(sizeof(char)*5);
    strcpy(in->nout,"50");
    in->partfile  = malloc(sizeof(char)*100);
    strcpy(in->partfile,"0");
    in->t0       = 0.0;
    in->tf       = 1.0e5 * YR;
    in->partsize = 1.0; // cm
    in->partdens = PARTDENSITY;
    in->nparts = 100;
    in->diffusion = 1;
    in->residenceTimes = 0;
    in->velocities = 0;
    return in;
}

void free_Inputs(Inputs *in) {
    free (in->fargodir);
    free (in->outputdir);
    free (in->nout);
    free (in->partfile);
    free (in);
}

int read_bool(const char *bin) {
    if (strcmp(bin,"1")==0) {
        return 1;
    } else if (strcmp(bin,"0")==0) {
        return 0;
    }

    if (strcmp(bin,"TRUE")==0) {
        return 1;
    } else if (strcmp(bin,"FALSE")==0) {
        return 0;
    }
    if (bin[0]=='Y') {
        return 1;
    } else if (bin[1]=='N') {
        return 0;
    }
    printf("Cannot parse boolean: %s\n",bin);
    exit(1);
}

Inputs *read_inputs(const char* infile) {
    Inputs *in = init_Inputs();
    FILE *file;
    file = fopen(infile, "r");
    if (file == NULL) {
        perror("Cannot open input file");
        exit(1);
    }
     char* line = NULL;
    size_t len=0;
    ssize_t read=0;
    while ((read = getline(&line, &len, file)) != -1) {
        // printf("Read line : %s",line);
        char *split_str;
        split_str = strtok(line, " \t\n");
        char *key;
        char *val_s;
        size_t idx = 0;
        while (split_str != NULL) {
            if (idx == 0) {
                key = split_str;
            } else if (idx == 1) {
                val_s = split_str;
            }
            idx++;
            split_str = strtok(NULL," \t\n");
        }
        if (strcmp(key, "FARGODIR") == 0) {
            strcpy(in->fargodir,val_s);
        } else if (strcmp(key, "OUTPUTDIR") == 0) {
            strcpy(in->outputdir,val_s);
        } else if (strcmp(key,"NOUT") == 0) {
            strcpy(in->nout,val_s);
        } else if (strcmp(key,"PARTFILE") == 0) {
            strcpy(in->partfile,val_s);
        } else if (strcmp(key,"T0") == 0) {
            in->t0 = atof(val_s);
        } else if (strcmp(key,"TF") == 0) {
            in->tf = atof(val_s);
        } else if (strcmp(key,"PARTSIZE") == 0) {
            in->partsize = atof(val_s);
        } else if (strcmp(key,"PARTDENS") == 0) {
            in->partdens = atof(val_s);
        } else if (strcmp(key,"NPARTS") == 0) {
            in->nparts = atoi(val_s);
        } else if (strcmp(key,"DIFFUSION") == 0) {
            in->diffusion = read_bool(val_s);
        } else if (strcmp(key,"RESIDENCETIMES") == 0) {
            in->residenceTimes = read_bool(val_s);
        } else if (strcmp(key,"VELOCITIES") == 0) {
            in->velocities = read_bool(val_s);
        }
    }
    return in;
}

void print_Inputs(Inputs *in) {
    printf("Strings: \n");
    printf("  FARGODIR  : %s\n",in->fargodir);
    printf("  OUTPUTDIR : %s\n",in->outputdir);
    printf("  NOUT      : %s\n",in->nout);
    printf("  PARTFILE  : %s\n",in->partfile);
    printf("Doubles: \n");
    printf("  T0       : %f\n",in->t0);
    printf("  TF       : %f\n",in->tf);
    printf("  PARTSIZE : %f\n",in->partsize);
    printf("  PARTDENS : %f\n",in->partdens);
    printf("Integers: \n");
    printf("  NPARTS : %d\n",in->nparts);
    printf("Booleans: \n");
    printf("  DIFFUSION      : %d\n",in->diffusion);
    printf("  RESIDENCETIMES : %d\n",in->residenceTimes);
    printf("  VELOCITIES     : %d\n",in->velocities);
}

/**
 * @brief Make directory, copy pasted from Jinku Hu @ delftstack.com
 * 
 * @param name 
 * @return int 
 */
int makedir(const char *name) {
  errno = 0;
  int ret = mkdir(name, S_IRWXU);
  if (ret == -1) {
    switch (errno) {
      case EACCES:
        printf("the parent directory does not allow write\n");
        return -1;
      case EEXIST:
        printf("pathname already exists\n");
        return 1;
      case ENAMETOOLONG:
        printf("pathname is too long\n");
        return -1;
      default:
        perror("Error in mkdir\n");
        return -1;
    }
  }

  return 0;
}