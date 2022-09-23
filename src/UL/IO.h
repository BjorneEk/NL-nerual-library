/*==========================================================*
 *
 * @author Gustaf Franzén :: https://github.com/BjorneEk;
 *
 * IO header
 *
 *==========================================================*/

#ifndef _IO_H_
#define _IO_H_

#include "types.h"

/**
 *  allocates space and reads the file contents from the file into the
 *  res array, returns the length of the array
 **/
i64_t read_file(const char * filename, u8_t * res);
/**
 *  allocates space and reads the file contents from the file into the
 *  res string, returns the length of the string
 **/
i64_t read_sfile(const char * filename, char * res);

/**
 *  creates and writes the bytes in data to a file named filename
 **/
u8_t write_file(const char * filename, u8_t * data);
/**
 *  creates and writes the string in data to a file named filename
 **/
u8_t write_sfile(const char * filename, char * data);


#endif /* _IO_H_ */
