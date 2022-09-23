/*==========================================================*
 *
 * @author Gustaf Franzén :: https://github.com/BjorneEk;
 *
 * dense layer functions for the NL neural network library
 *
 *==========================================================*/

#include "NL_layer.h"
#include <stdio.h>

LA_matrix* NL_dense_forward(NL_layer self, LA_matrix* input)
{
  LA_matrix* res;
  LA_copy_r(self.input, input);
  res = LA_mult(self.dense.weights, input);
  LA_add_r(res, self.dense.bias);
  LA_free(input);
  return res;
}

LA_matrix* NL_dense_backward(NL_layer self, LA_matrix* Derr_Dout, f64_t lern_rate)
{
  LA_matrix* res;
  LA_matrix* input_T;
  LA_matrix* weights_T;
  LA_matrix* Derr_Dout_;
  LA_matrix* Derr_Dweights;

  input_T          = LA_transpose(self.input);
  weights_T        = LA_transpose(self.dense.weights);
  Derr_Dweights    = LA_mult(Derr_Dout, input_T);
  Derr_Dout_       = LA_scale(Derr_Dout, lern_rate);

  LA_scale_r(Derr_Dweights, lern_rate);

  LA_sub_r(self.dense.weights, Derr_Dweights);
  LA_sub_r(self.dense.bias, Derr_Dout_);

  LA_free(Derr_Dout_);
  LA_free(input_T);
  LA_free(Derr_Dweights);

  res = LA_mult(weights_T, Derr_Dout);
  LA_free(weights_T);
  LA_free(Derr_Dout);
  return res;
}

NL_layer NL_new_dense(u32_t input_size, u32_t output_size)
{
  NL_layer res;

  res.dense.in_size  = input_size;
  res.dense.out_size = output_size;
  res.forward_func   = NL_dense_forward;
  res.backwards_func = NL_dense_backward;
  res.type           = DENSE;
  res.input          = LA_init_empty(1, input_size);
  res.dense.bias     = LA_random(1, output_size, 0.0f, 1.0f);
  res.dense.weights  = LA_random(input_size, output_size, 0.0f, 1.0f);
  return res;
}
