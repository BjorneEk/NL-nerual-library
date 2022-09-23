/*==========================================================*
 *
 * @author Gustaf Franzén :: https://github.com/BjorneEk;
 *
 * activation layer functions for the NL neural network library
 *
 *==========================================================*/

#include "NL_layer.h"
#include <math.h>
#include <stdio.h>
f64_t tanh_d(f64_t arg)
{
  f64_t th;

  th = tanh(arg);
  return 1.0f - (th * th);
}

void tanh_activate(LA_matrix* m)
{
  LA_foreach(m, tanh);
}

void tanh_activate_d(LA_matrix* m)
{
  LA_foreach(m, tanh_d);
}

LA_matrix* NL_activation_forward(NL_layer self, LA_matrix* input)
{
  LA_matrix* res;
  LA_copy_r(self.input, input);
  res = LA_copy(input);
  self.activation.activation_func(res);
  LA_free(input);
  return res;
}

LA_matrix* NL_activation_backward(NL_layer self, LA_matrix* output_gradient, f64_t lern_rate)
{
  LA_matrix* res;

  self.activation.activation_derivative_func(self.input);

  res = LA_ewise_mult(output_gradient, self.input);
  LA_free(output_gradient);
  return res;
}

NL_layer NL_new_activation(u32_t size, void (*af)(LA_matrix*), void (*af_d)(LA_matrix*))
{
  NL_layer res;

  res.forward_func                          = NL_activation_forward;
  res.backwards_func                        = NL_activation_backward;
  res.type                                  = ACTIVATION;
  res.activation.size                       = size;
  res.activation.activation_func            = af;
  res.activation.activation_derivative_func = af_d;
  res.input                                 = LA_init_empty(1, size);

  return res;
}

NL_layer NL_new_tanh(u32_t size)
{
  return NL_new_activation(size, tanh_activate, tanh_activate_d);
}
