#include <assert.h>
#include <ctype.h>
#include <stdarg.h>
#include <stdlib.h>
#include <string.h>

#include "cfg.h"

struct _ConfigEntry {
    char* key;
    char* value;
    struct _ConfigEntry* prev;
    struct _ConfigEntry* next;
};

typedef struct _ConfigEntry* ConfigEntry;

struct _Config {
    int strip_quotes;
    ConfigEntry first;  // first configuration entry, alphabetically
    ConfigEntry recent; // last accessed configuration entry
};


ConfigEntry cfg_new_entry(char* key, char* value, ConfigEntry prev, ConfigEntry next) {
    ConfigEntry entry = (ConfigEntry) malloc(sizeof(struct _ConfigEntry));
    entry->key = key;
    entry->value = value;
    entry->prev = prev;
    if(prev)
        prev->next = entry;
    entry->next = next;
    if(next)
        next->prev = entry;
    return entry;
}

void cfg_destroy_entry(ConfigEntry entry) {
    free(entry->key);
    free(entry->value);
    if(entry->prev)
        entry->prev->next = entry->next;
    if(entry->next)
        entry->next->prev = entry->prev;
    free(entry);
}


ConfigEntry cfg_find_entry(Config cfg, const char* key, ConfigEntry* prev, ConfigEntry* next) {
    int cmp, dir;
    ConfigEntry entry;
    if(cfg->first == NULL) {
        /* Empty database */
        if(prev) *prev = NULL;
        if(next) *next = NULL;
        return NULL;
    }
    else {
        if(cfg->recent == NULL)
            cfg->recent = cfg->first;
        entry = cfg->recent;

        /* Check most recently accessed configuration option. The sign of the
         * comparison result tells us which direction to traverse the database. */
        dir = cmp = strcmp(key, entry->key);

        if(dir == 0) {
            /* Direct cache hit */
            if(prev) *prev = entry->prev;
            if(next) *next = entry->next;
            cfg->recent = entry;
            return entry;
        }
        else if(dir > 0) {
            /* Traverse the list forward */
            while(entry->next != NULL && (cmp = strcmp(key, entry->next->key)) > 0)
                entry = entry->next;
            if(cmp == 0) {
                /* Exact match */
                entry = entry->next;
                if(prev) *prev = entry->prev;
                if(next) *next = entry->next;
                cfg->recent = entry;
                return entry;
            }
            else {
                /* No exact match */
                if(prev) *prev = entry;
                if(next) *next = entry->next;
                cfg->recent = entry;
                return NULL;
            }
        }
        else {
            /* Traverse the list backward */
            while(entry->prev != NULL && (cmp = strcmp(key, entry->prev->key)) < 0)
                entry = entry->prev;
            if(cmp == 0) {
                /* Exact match */
                entry = entry->prev;
                if(prev) *prev = entry->prev;
                if(next) *next = entry->next;
                cfg->recent = entry;
                return entry;
            }
            else {
                /* No exact match */
                if(prev) *prev = entry->prev;
                if(next) *next = entry;
                cfg->recent = entry;
                return NULL;
            }
        }
    }
}

void cfg_insert_entry(Config cfg, char* key, char* value) {
    ConfigEntry entry, prev, next;
    entry = cfg_find_entry(cfg, key, &prev, &next);
    if(entry != NULL) {
        cfg_destroy_entry(entry);
        if(entry == cfg->first)
            cfg->first = next;
        if(entry == cfg->recent)
            cfg->recent = prev;
    }
    entry = cfg_new_entry(key, value, prev, next);
    if(prev == NULL)
        cfg->first = entry;
    if(cfg->recent == NULL)
        cfg->recent = entry;
}

/******************************************************************************/

Config cfg_new() {
    Config cfg = (Config) malloc(sizeof(struct _Config));
    cfg->strip_quotes = 1;
    cfg->first = NULL;
    cfg->recent = NULL;
    return cfg;
}

Config cfg_new_from_file(const char* filename) {
    Config cfg = cfg_new();
    if(cfg_read_file(cfg, filename) != CFG_OKAY) {
        cfg_destroy(cfg);
        cfg = NULL;
    }
    return cfg;
}

Config cfg_new_from_stream(FILE* f) {
    Config cfg = cfg_new();
    if(cfg_read(cfg, f) != CFG_OKAY) {
        cfg_destroy(cfg);
        cfg = NULL;
    }
    return cfg;
}

Config cfg_new_copy(Config orig) {
    ConfigEntry entry;
    Config cfg = cfg_new();
    cfg->strip_quotes = orig->strip_quotes;
    entry = orig->first;
    while(entry != NULL) {
        cfg_insert_entry(cfg, strdup(entry->key), strdup(entry->value));
        entry = entry->next;
    }
    return cfg;
}

Config cfg_new_sub(Config cfg, const char* prefix, int strip) {
    int prefixlen = strlen(prefix);
    ConfigEntry entry;
    Config subcfg = cfg_new();
    subcfg->strip_quotes = cfg->strip_quotes;
    entry = cfg->first;
    while(entry != NULL) {
        if(strncmp(entry->key, prefix, prefixlen) == 0) {
            if(strip)
                cfg_insert_entry(subcfg, strdup(entry->key + prefixlen), strdup(entry->value));
            else
                cfg_insert_entry(subcfg, strdup(entry->key), strdup(entry->value));
        }
        entry = entry->next;
    }
    return subcfg;
}

void cfg_destroy(Config cfg) {
    if(cfg != NULL) {
        ConfigEntry next, entry = cfg->first;
        while(entry != NULL) {
            next = entry->next;
            free(entry->key);
            free(entry->value);
            free(entry);
            entry = next;
        }
        free(cfg);
    }
}

int cfg_read(Config cfg, FILE* f) {
    int err = 0;
    char line[CFG_MAX_LINE_LENGTH];
    if(f == NULL) {
        fprintf(stderr, "Config: file pointer is NULL\n");
        return CFG_STREAM_ERROR;
    }
    while(fgets(line, sizeof(line), f) != NULL)
        err &= cfg_read_line(cfg, line);
    return err;
}

int cfg_read_file(Config cfg, const char* filename) {
    int err;
    FILE* f = fopen(filename, "r");
    if(f == NULL) {
        fprintf(stderr, "Config: could not read file '%s'\n", filename);
        return CFG_FILE_ERROR;
    }
    err = cfg_read(cfg, f);
    fclose(f);
    return err;
}

/* This is a quick replacement for strndup(), which is not available on all systems */
static char *my_strndup(const char *s, size_t n) {
    char *r = strdup(s);
    if(strlen(s) > n)
        r[n] = '\0';
    return r;
}

int cfg_read_line(Config cfg, const char* line) {
    size_t n;
    const char *p, *end;
    const char *eqpos, *keyend;
    const char *valbegin, *valend;
    char *key, *val;

    n = strlen(line);
    p = line;
    end = &line[n];

    /* Strip leading whitespace */
    while(p != end && isspace(*p))
        p++;

    /* Ignore blank lines or lines that start with '#' */
    if(p == end || p[0] == '#')
        return 0;

    /* Split string on '=' character */
    eqpos = strchr(p, '=');
    if(eqpos == NULL || eqpos == p) {
        fprintf(stderr, "Config: '%s' is not a valid configuration line\n", line);
        return CFG_PARSE_ERROR;
    }

    /* Strip whitespace from key name */
    keyend = eqpos;
    while(isspace(*(keyend - 1)))
        keyend--;
    key = my_strndup(p, keyend-p);

    /* Strip whitespace (and optionally quotes) from value */
    valbegin = eqpos + 1;
    while(valbegin != end && isspace(*valbegin))
        valbegin++;
    valend = end;
    while(isspace(*(valend - 1)))
        valend--;
    if(cfg->strip_quotes && valend - valbegin >= 2
            && ((*valbegin == '"' && *(valend-1) == '"')
             || (*valbegin == '\'' && *(valend-1) == '\'')))
    {
        valbegin++;
        valend--;
    }
    val = (valbegin >= valend) ? strdup("") : my_strndup(valbegin, valend-valbegin);

    cfg_insert_entry(cfg, key, val);

    return CFG_OKAY;
}

/* Check whether or not the string starts or ends with whitespace */
static int needs_quotes(const char* s) {
    size_t n = strlen(s);
    return (n > 0 && (isspace(s[0]) || isspace(s[n-1])));
}

int cfg_write(Config cfg, FILE* f) {
    ConfigEntry entry;
    if(f == NULL) {
        fprintf(stderr, "Config: file pointer is NULL\n");
        return CFG_STREAM_ERROR;
    }
    entry = cfg->first;
    while(entry != NULL) {
        if(needs_quotes(entry->value))
            fprintf(f, "%s = \"%s\"\n", entry->key, entry->value);
        else
            fprintf(f, "%s = %s\n", entry->key, entry->value);
        entry = entry->next;
    }
    return CFG_OKAY;
}

int cfg_write_file(Config cfg, const char* filename, const char* mode) {
    FILE* f = fopen(filename, mode);
    if(f == NULL) {
        fprintf(stderr, "Config: could not write to file '%s'\n", filename);
        return CFG_FILE_ERROR;
    }
    cfg_write(cfg, f);
    fclose(f);
    return CFG_OKAY;
}

static char valbuf[CFG_MAX_VALUE_LENGTH];

void cfg_set(Config cfg, const char* key, const char* value) {
    cfg_insert_entry(cfg, strdup(key), strdup(value));
}

#define DEFINE_CFG_SET(name, type, fmt) \
void cfg_set_##name(Config cfg, const char* key, type value) { \
    char valbuf[CFG_MAX_VALUE_LENGTH]; \
    snprintf(valbuf, CFG_MAX_VALUE_LENGTH, fmt, value); \
    cfg_set(cfg, key, valbuf); \
}
DEFINE_CFG_SET(char, char, "%c")
DEFINE_CFG_SET(short, short, "%hi")
DEFINE_CFG_SET(int, int, "%i")
DEFINE_CFG_SET(long, long, "%li")
DEFINE_CFG_SET(uchar, unsigned char, "%c")
DEFINE_CFG_SET(ushort, unsigned short, "%hu")
DEFINE_CFG_SET(uint, unsigned int, "%u")
DEFINE_CFG_SET(ulong, unsigned long, "%lu")
DEFINE_CFG_SET(float, float, "%g")
DEFINE_CFG_SET(double, double, "%g")
#undef DEFINE_CFG_SET

void cfg_set_format(Config cfg, const char* key, const char* fmt, ...) {
    char valbuf[CFG_MAX_VALUE_LENGTH];
    va_list ap;
    va_start(ap, fmt);
    vsnprintf(valbuf, CFG_MAX_VALUE_LENGTH, fmt, ap);
    va_end(ap);
    cfg_insert_entry(cfg, strdup(key), strdup(valbuf));
}

int cfg_has_key(Config cfg, const char* key) {
    return cfg_find_entry(cfg, key, NULL, NULL) != NULL;
}

int cfg_has_keys(Config cfg, const char* keys, const char* sep) {
    int hasall = 1;
    char* s = strdup(keys);
    char* token = strtok(s, sep);
    while(token != NULL) {
        hasall &= cfg_has_key(cfg, token);
        token = strtok(NULL, sep);
    }
    free(s);
    return hasall;
}

const char* cfg_get(Config cfg, const char* key) {
    ConfigEntry e = cfg_find_entry(cfg, key, NULL, NULL);
    return (e == NULL) ? "" : e->value;
}

#define DEFINE_CFG_GET(name, type, fmt) \
type cfg_get_##name(Config cfg, const char* key) { \
    type value = 0; \
    const char* s = cfg_get(cfg, key); \
    int nread = sscanf(s, fmt, &value); \
    assert(nread == 1); \
    return value; \
}
DEFINE_CFG_GET(char, char, "%c")
DEFINE_CFG_GET(short, short, "%hi")
DEFINE_CFG_GET(int, int, "%i")
DEFINE_CFG_GET(long, long, "%li")
DEFINE_CFG_GET(uchar, unsigned char, "%c")
DEFINE_CFG_GET(ushort, unsigned short, "%hu")
DEFINE_CFG_GET(uint, unsigned int, "%u")
DEFINE_CFG_GET(ulong, unsigned long, "%lu")
DEFINE_CFG_GET(float, float, "%f")
DEFINE_CFG_GET(double, double, "%lf")
#undef DEFINE_CFG_GET

#define DEFINE_CFG_GET_ARRAY(name, type, fmt) \
int cfg_get_array_##name(Config cfg, const char* key, int nmax, type* values) { \
    int n = 0, nscanned; \
    char* orig = strdup(cfg_get(cfg, key)); \
    char* s = strtok(orig, " ,"); \
    while(n < nmax && s != NULL) { \
        nscanned = sscanf(s, fmt, &values[n]); \
        assert(nscanned == 1); \
        n++; \
        s = strtok(NULL, " ,"); \
    } \
    free(orig); \
    return n; \
}
DEFINE_CFG_GET_ARRAY(char, char, "%c")
DEFINE_CFG_GET_ARRAY(short, short, "%hi")
DEFINE_CFG_GET_ARRAY(int, int, "%i")
DEFINE_CFG_GET_ARRAY(long, long, "%li")
DEFINE_CFG_GET_ARRAY(uchar, unsigned char, "%c")
DEFINE_CFG_GET_ARRAY(ushort, unsigned short, "%hu")
DEFINE_CFG_GET_ARRAY(uint, unsigned int, "%u")
DEFINE_CFG_GET_ARRAY(ulong, unsigned long, "%lu")
DEFINE_CFG_GET_ARRAY(float, float, "%f")
DEFINE_CFG_GET_ARRAY(double, double, "%lf")
#undef DEFINE_CFG_GET_ARRAY

int cfg_get_format(Config cfg, const char* key, const char* fmt, ...) {
    int nread;
    va_list ap;
    const char* s = cfg_get(cfg, key);
    va_start(ap, fmt);
    nread = vsscanf(s, fmt, ap);
    va_end(ap);
    return nread;
}

#if 0
/* Get the specified configuration option, returning an error code if the
 * configuration key is not present or if the associated value cannot be
 * converted to the desired type. */
int cfg_get_e       (Config cfg, const char* key, const char** value);
int cfg_get_char_e  (Config cfg, const char* key, char* value);
int cfg_get_short_e (Config cfg, const char* key, short* value);
int cfg_get_int_e   (Config cfg, const char* key, int* value);
int cfg_get_long_e  (Config cfg, const char* key, long* value);
int cfg_get_uchar_e (Config cfg, const char* key, unsigned char* value);
int cfg_get_ushort_e(Config cfg, const char* key, unsigned short* value);
int cfg_get_uint_e  (Config cfg, const char* key, unsigned int* value);
int cfg_get_ulong_e (Config cfg, const char* key, unsigned long* value);
#endif
