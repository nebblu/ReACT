#ifndef CFG_H
#define CFG_H

/* TODO:
 * - rename cfg_set_format to cfg_set_printf and add cfg_get_scanf
 * */

#ifdef __cplusplus
extern "C" {
#endif

#include <stdio.h>

/* Defines */
#define CFG_MAX_LINE_LENGTH 8192
#define CFG_MAX_VALUE_LENGTH 2048
//#define CFG_THREAD_SAFE 0

/* Error codes */
#define CFG_OKAY 0
#define CFG_MISSING_KEY_ERROR 1
#define CFG_WRONG_TYPE_ERROR 2
#define CFG_STREAM_ERROR 3
#define CFG_FILE_ERROR 4
#define CFG_PARSE_ERROR 5

/* Opaque data type for configuration databases */
typedef struct _Config*  Config;

/* Create an empty configuration database */
Config cfg_new();

/* Create a new configuration database from the contents of the given file */
Config cfg_new_from_file(const char* filename);
Config cfg_new_from_stream(FILE* f);

/* Create a new configuration database from an existing one */
Config cfg_new_copy(Config orig);

/* Create a new configuration database from an existing one, keeping only those
 * entries that begin with the specified prefix.  If 'strip' is nonzero, the
 * prefix will be stripped from all entries. */
Config cfg_new_sub(Config cfg, const char* prefix, int strip);

/* Destroy a configuration database and free all associated memory */
void cfg_destroy(Config cfg);

/* Read additional configuration options */
int cfg_read     (Config cfg, FILE* f);
int cfg_read_file(Config cfg, const char* filename);
int cfg_read_line(Config cfg, const char* line);

int cfg_write     (Config cfg, FILE* f);
int cfg_write_file(Config cfg, const char* filename, const char* mode);

#if 0
/* Append the contents of one configuration database to another */
void cfg_append(Config dst, Config src);
#endif

/* Set the specified configuration option */
void cfg_set       (Config cfg, const char* key, const char* value);
void cfg_set_char  (Config cfg, const char* key, char value);
void cfg_set_short (Config cfg, const char* key, short value);
void cfg_set_int   (Config cfg, const char* key, int value);
void cfg_set_long  (Config cfg, const char* key, long value);
void cfg_set_uchar (Config cfg, const char* key, unsigned char value);
void cfg_set_ushort(Config cfg, const char* key, unsigned short value);
void cfg_set_uint  (Config cfg, const char* key, unsigned int value);
void cfg_set_ulong (Config cfg, const char* key, unsigned long value);
void cfg_set_float (Config cfg, const char* key, float value);
void cfg_set_double(Config cfg, const char* key, double value);
void cfg_set_format(Config cfg, const char* key, const char* fmt, ...);

/* Return 1 if the specified key exists, 0 otherwise */
int cfg_has_key(Config cfg, const char* key);
int cfg_has_keys(Config cfg, const char* keys, const char* sep);

/* Get the specified configuration option, without error checking */
const char*    cfg_get       (Config cfg, const char* key);
char           cfg_get_char  (Config cfg, const char* key);
short          cfg_get_short (Config cfg, const char* key);
int            cfg_get_int   (Config cfg, const char* key);
long           cfg_get_long  (Config cfg, const char* key);
unsigned char  cfg_get_uchar (Config cfg, const char* key);
unsigned short cfg_get_ushort(Config cfg, const char* key);
unsigned int   cfg_get_uint  (Config cfg, const char* key);
unsigned long  cfg_get_ulong (Config cfg, const char* key);
float          cfg_get_float (Config cfg, const char* key);
double         cfg_get_double(Config cfg, const char* key);

int cfg_get_format(Config cfg, const char* key, const char* fmt, ...);

/* Parse the given configuration option as an array of length at most nmax.
 * Return the number of values actually read. */
int cfg_get_array_int   (Config cfg, const char* key, int nmax, int* values);
int cfg_get_array_float (Config cfg, const char* key, int nmax, float* values);
int cfg_get_array_double(Config cfg, const char* key, int nmax, double* values);

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
int cfg_get_float_e (Config cfg, const char* key, float* value);
int cfg_get_double_e(Config cfg, const char* key, double* value);
#endif

#ifdef __cplusplus
}
#endif

#endif // CFG_H
