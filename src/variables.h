#include <stdlib.h>
#include <stdio.h>
#include <string.h>

// Maximum number of items that can be saved in the dictionary
#define MAX_VAR_SIZE 100
// maixmum number of characters in the dictionary
#define MAX_KEY_SIZE  100
// number of keys we want to read
#define NKEYS 5

/**
 * @brief Dictionary object to store variables from the fargo model.
 * 
 * Shout out to Risha on GeeksforGeeks for the dictionary boilerplate I based this off of.
 * 
 */
typedef struct Variables {
    int size;
    char **keys;
    double *values;
} Variables;

/**
 * @brief Initialize empty dictionary for variables
 * 
 * @return Variables* 
 */
Variables *init_Variables() {
    Variables *var = malloc(sizeof(*var));
    var->size = 0;
    var->keys = malloc(sizeof(char*) * MAX_VAR_SIZE);
    var->values = malloc(sizeof(double*) * MAX_VAR_SIZE);
    return var;
}

/**
 * @brief free memory associated with Variables
 * 
 * @param var 
 */
void free_Variables(Variables *var) {
    free(var->keys);
    free(var->values);
    free(var);
}

/**
 * @brief add a variable to the dictionary
 * 
 * @param var 
 * @param key 
 * @param value 
 * @return int 
 */
int add_variable(Variables *var, char* key, double value) {
    if (var->size==MAX_KEY_SIZE-1) {
        perror("Variables dictionary is full!");
        exit(1);
    }
    var->keys[var->size] = malloc(sizeof(char)*MAX_KEY_SIZE);
    strcpy(var->keys[var->size], key);
    var->values[var->size] = value;
    var->size++;
    return 0;
}

/**
 * @brief Get the index of a key from the dictionary
 * 
 * @param var 
 * @param key 
 * @return int 
 */
int get_index(Variables *var, char* key) {
    for (int i=0; i<var->size; i++) {
        if (strcmp(var->keys[i],key)==0) {
            return i;
        }
    }
    // key not found
    return -1;
}

/**
 * @brief Get the value of a given key in the dictionary
 * 
 * @param var 
 * @param key 
 * @return double 
 */
double get_value(Variables *var, char* key) {
    int index = get_index(var, key);
    if (index==-1) {
        printf("Key not found: %s\n", key);
        exit(1);
    }
    return var->values[index];
}

/**
 * @brief print the keys and values from the dictionary
 * 
 * @param var 
 */
void print_variables(Variables *var) {
    printf("[\n");
    for (int i=0; i<var->size; i++) {
        printf("  %s : %f",var->keys[i],var->values[i]);
        if (i<var->size-1) {
            printf(",\n");
        }
    }
    printf("\n]\n");
}

Variables *init_Variables_fromFile(char *fname) {
    Variables *var = init_Variables();
    FILE *file;
    file = fopen(fname, "r");
    if (file == NULL) {
        perror("Cannot open variables file");
        exit(1);
    }
    // list of keys that we want to find
    char *key_search[NKEYS] = {
        "ALPHA", "ASPECTRATIO", "FLARINGINDEX", "PLANETMASS", "OMEGAFRAME"
    };
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
        for (int i=0; i<NKEYS; i++) {
            if (strcmp(key,key_search[i])==0) {
                add_variable(var,key,atof(val_s));
            }
        }
    }

    return var;
}