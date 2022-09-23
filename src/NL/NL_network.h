
/*==========================================================*
 *
 * @author Gustaf Franzén :: https://github.com/BjorneEk;
 *
 * neural network header file for the NL neural network library
 *
 *==========================================================*/

#ifndef _NL_NETWORK_H_
#define _NL_NETWORK_H_

#include "NL_layer.h"
#include <math.h>
#include <stdarg.h>


typedef struct neural_network {

  u32_t layer_count;

  NL_layer * layers;

} NL_network;

f64_t NL_mse(LA_matrix * real, LA_matrix * expected);

LA_matrix * NL_mse_d(LA_matrix * real, LA_matrix * expected);

NL_network NL_new_network(u32_t layer_count, ...);

void NL_free(NL_network self);

#endif /* _NL_NETWORK_H_ */
