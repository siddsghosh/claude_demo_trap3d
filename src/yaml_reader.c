#include "yaml_reader.h"
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

#define MAX_LINE_LENGTH 1024

/**
 * Trim leading and trailing whitespace from a string
 */
static char* trim(char *str) {
    char *end;

    // Trim leading space
    while(isspace((unsigned char)*str)) str++;

    if(*str == 0) return str;

    // Trim trailing space
    end = str + strlen(str) - 1;
    while(end > str && isspace((unsigned char)*end)) end--;

    // Write new null terminator
    end[1] = '\0';

    return str;
}

/**
 * Remove quotes from a string
 */
static void remove_quotes(char *str) {
    size_t len = strlen(str);
    if (len >= 2 && ((str[0] == '"' && str[len-1] == '"') ||
                     (str[0] == '\'' && str[len-1] == '\''))) {
        memmove(str, str + 1, len - 2);
        str[len - 2] = '\0';
    }
}

/**
 * Read a double value from a YAML file
 */
int yaml_read_double(const char *filename, const char *key, double *value) {
    FILE *fp = fopen(filename, "r");
    if (!fp) {
        fprintf(stderr, "Error: Cannot open file %s\n", filename);
        return 1;
    }

    char line[MAX_LINE_LENGTH];
    int found = 0;

    while (fgets(line, sizeof(line), fp)) {
        // Skip comments and empty lines
        char *trimmed = trim(line);
        if (trimmed[0] == '#' || trimmed[0] == '\0') continue;

        // Look for key: value pattern
        char *colon = strchr(trimmed, ':');
        if (!colon) continue;

        // Extract key and value
        *colon = '\0';
        char *file_key = trim(trimmed);
        char *file_value = trim(colon + 1);

        // Check if this is the key we're looking for
        if (strcmp(file_key, key) == 0) {
            // Remove quotes if present
            remove_quotes(file_value);

            // Try to parse as double
            char *endptr;
            *value = strtod(file_value, &endptr);
            if (*endptr != '\0' && !isspace((unsigned char)*endptr)) {
                fprintf(stderr, "Error: Cannot parse value '%s' as double for key '%s'\n",
                        file_value, key);
                fclose(fp);
                return 1;
            }
            found = 1;
            break;
        }
    }

    fclose(fp);

    if (!found) {
        fprintf(stderr, "Error: Key '%s' not found in %s\n", key, filename);
        return 1;
    }

    return 0;
}

/**
 * Read a string value from a YAML file
 */
int yaml_read_string(const char *filename, const char *key, char **value) {
    FILE *fp = fopen(filename, "r");
    if (!fp) {
        fprintf(stderr, "Error: Cannot open file %s\n", filename);
        return 1;
    }

    char line[MAX_LINE_LENGTH];
    int found = 0;

    while (fgets(line, sizeof(line), fp)) {
        // Skip comments and empty lines
        char *trimmed = trim(line);
        if (trimmed[0] == '#' || trimmed[0] == '\0') continue;

        // Look for key: value pattern
        char *colon = strchr(trimmed, ':');
        if (!colon) continue;

        // Extract key and value
        *colon = '\0';
        char *file_key = trim(trimmed);
        char *file_value = trim(colon + 1);

        // Check if this is the key we're looking for
        if (strcmp(file_key, key) == 0) {
            // Remove quotes if present
            remove_quotes(file_value);

            // Allocate and copy string
            *value = strdup(file_value);
            if (!*value) {
                fprintf(stderr, "Error: Memory allocation failed\n");
                fclose(fp);
                return 1;
            }
            found = 1;
            break;
        }
    }

    fclose(fp);

    if (!found) {
        fprintf(stderr, "Error: Key '%s' not found in %s\n", key, filename);
        return 1;
    }

    return 0;
}
