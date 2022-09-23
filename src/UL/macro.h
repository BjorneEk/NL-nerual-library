/*==========================================================*
 *
 * @author Gustaf Franzén :: https://github.com/BjorneEk;
 *
 * general macros
 *
 *==========================================================*/

#ifndef _MACRO_H_
#define _MACRO_H_

#define PACKED __attribute__((__packed__))

#define CAT(a, ...) PRIMITIVE_CAT(a, __VA_ARGS__)
#define PRIMITIVE_CAT(a, ...) a ## __VA_ARGS__

#endif /* _MACRO_H_ */
