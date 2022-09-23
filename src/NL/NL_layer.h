
/*==========================================================*
 *
 * @author Gustaf Franzén :: https://github.com/BjorneEk;
 *
 * layer header file for the NL neural network library
 *
 *==========================================================*/

#ifndef _NL_LAYER_H_
#define _NL_LAYER_H_

#include "../UL/types.h"
#include "../LA/matrix.h"

struct NL_dense_layer {

  u32_t in_size, out_size;
  LA_matrix * weights;

  LA_matrix * bias;

};

struct NL_activation_layer {

  u32_t size;
  void (*activation_func)(LA_matrix*);
  void (*activation_derivative_func)(LA_matrix*);
};

typedef enum layertype {
  DENSE,
  ACTIVATION
} LN_layertype;


typedef struct neural_network_layer {
  LN_layertype type;

  LA_matrix* input;

  LA_matrix* (*forward_func)(struct neural_network_layer, LA_matrix*);

  LA_matrix* (*backwards_func)(struct neural_network_layer, LA_matrix*, f64_t);

  union {
    struct NL_dense_layer      dense;
    struct NL_activation_layer activation;
  };
} NL_layer;


/**
 * creates a new NL_layer of type DENSE
 **/
NL_layer NL_new_dense(u32_t input_size, u32_t output_size);
/**
 * creates a new NL_layer of type ACTIVATION
 **/
NL_layer NL_new_activation(u32_t size, void (*af)(LA_matrix*), void (*af_d)(LA_matrix*));


/**
 * forward propagation function for the DENSE layer type
 **/
LA_matrix* NL_dense_forward(NL_layer self, LA_matrix* input);

/**
 * forward propagation function for the DENSE layer type
 **/
LA_matrix* NL_dense_backward(NL_layer self, LA_matrix* output_gradient, f64_t lern_rate);

/**
 * forward propagation function for the ACTIVATION layer type
 **/
LA_matrix* NL_activation_forward(NL_layer self, LA_matrix* input);

/**
 * forward propagation function for the NL_dense layer type
 **/
LA_matrix* NL_activation_backward(NL_layer self, LA_matrix* output_gradient, f64_t lern_rate);

NL_layer NL_new_tanh(u32_t size);


#endif /* _NL_LAYER_H_ */
