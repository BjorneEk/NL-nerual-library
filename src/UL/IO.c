/*==========================================================*
 *
 * @author Gustaf Franzén :: https://github.com/BjorneEk;
 *
 * IO functions
 *
 *==========================================================*/
#include <stdlib.h>
#include <string.h>
#include "IO.h"
#include "assert.h"

i64_t read_file(const char * filename, u8_t * res)
{
  i64_t new_len;
  i64_t buffer_size;
  FILE * fp;

  fp = fopen(filename, "r");
  assert(fp != NULL, "could not read from file");
  /* Go to the end of the file. */
  if (fseek(fp, 0L, SEEK_END) == 0) {
    /* Get the size of the file. */
    buffer_size = ftell(fp);
    assert(buffer_size != -1, "Error reading file");
    /* Allocate our buffer to that size. */
    res = malloc(sizeof(u8_t) * (buffer_size + 1));
    assert(res != NULL, "Out of memory");

    /* Go back to the start of the file. */
    assert(fseek(fp, 0L, SEEK_SET) == 0, "Error reading file");

    /* Read the entire file into memory. */
    new_len = fread(res, sizeof(char), buffer_size, fp);

    assert(ferror( fp ) == 0, "Error reading file");
    res[new_len++] = '\0'; /* Just to be safe. */
  }
  fclose(fp);
  return new_len;
}
i64_t read_sfile(const char * filename, char * res)
{
  return read_file(filename, (u8_t *) res);
}
