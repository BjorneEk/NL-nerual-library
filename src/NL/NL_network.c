/*==========================================================*
 *
 * @author Gustaf Franzén :: https://github.com/BjorneEk;
 *
 * neural network functions for the NL neural network library
 *
 *==========================================================*/


#include <stdlib.h>
#include "../UL/assert.h"
#include "NL_network.h"



f64_t pow2(f64_t f)
{
   return pow(f, 2);
}


f64_t NL_mse(LA_matrix * expected, LA_matrix * real)
{
  f64_t res;
  LA_matrix * tmp;

  tmp = LA_sub(expected, real);
  LA_foreach(tmp, pow2);
  res = LA_mean(tmp);
  LA_free(tmp);

  return res;
}

LA_matrix * NL_mse_d(LA_matrix * expected, LA_matrix * real)
{
  LA_matrix * tmp;

  tmp = LA_sub(real, expected);
  LA_scale_r(tmp, 2.0f/(tmp->cols*tmp->rows));

  return tmp;
}

NL_network NL_new_network(u32_t layer_count, ...)
{
  i32_t i;
  NL_network res;
  va_list layers;

  res.layer_count = layer_count;
  va_start(layers, layer_count);

  res.layers = malloc(layer_count * sizeof(NL_layer));
  assert(res.layers != NULL, "Out of Memory(creating neural network)");

  for(i = 0; i < layer_count; i++) res.layers[i] = va_arg(layers, NL_layer);
  va_end(layers);
  return res;
}

void NL_free(NL_network self)
{
  if(self.layers != NULL) free(self.layers);
}
