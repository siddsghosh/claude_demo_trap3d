#ifndef YAML_READER_H
#define YAML_READER_H

#include <stdio.h>

/**
 * Read a double value from a YAML file
 * Returns 0 on success, non-zero on failure
 */
int yaml_read_double(const char *filename, const char *key, double *value);

/**
 * Read a string value from a YAML file
 * Returns 0 on success, non-zero on failure
 * Caller must free the returned string
 */
int yaml_read_string(const char *filename, const char *key, char **value);

#endif /* YAML_READER_H */
