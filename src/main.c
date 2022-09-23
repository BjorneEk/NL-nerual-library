
////////////////////////////////////////////////////////////////////////////
///        @author Gustaf Franzén :: https://github.com/BjorneEk;        ///
////////////////////////////////////////////////////////////////////////////

#include <stdio.h>
#include <stdlib.h>
#include "NL/NL_network.h"
#include "UL/math.h"



int main(int argc, char const *argv[]) {

  LA_matrix *X, *Y, *output, * grad;
  NL_network nn;
  i32_t epochs, i, j, x, y;
  f64_t lern_rate, err;

  X = LA_new_mat(2, 4, 0.0, 0.0, 0.0, 1.0, 1.0, 0.0, 1.0, 1.0);
  Y = LA_new_mat(1, 4, 0.0, 1.0, 1.0, 0.0);

  epochs = 1000000;
  lern_rate = 0.1;

  nn = NL_new_network(4,
    NL_new_dense(2, 3),
    NL_new_tanh(3),
    NL_new_dense(3, 1),
    NL_new_tanh(1));


  for(i = 0; i < epochs; i++) {
    err = 0;
      for(x=0,y=0; y < Y->rows; y++, x++) {
        output = LA_rowat(X, x);
        SWAP(output->cols, output->rows);
        for(j = 0; j < nn.layer_count; j++) {
          output = nn.layers[j].forward_func(nn.layers[j], output);
        }
        err += NL_mse(&(LA_matrix){1,1,&Y->data[y]}, output);
        j--;
        grad = NL_mse_d(&(LA_matrix){1,1,&Y->data[y]}, output);
        for(;j >= 0; j--) {
          grad = nn.layers[j].backwards_func(nn.layers[j], grad, lern_rate);
        }
        LA_free(output);
        LA_free(grad);
      }
    err /= (f64_t)(X->rows);
  }
  printf("error:%f\n", err);
  printf("training finnished, starting tests\n");
  for(i=0;i<4;i++) {
    output = LA_rowat(X, i);

    SWAP(output->cols, output->rows);
    for(j = 0; j < nn.layer_count; j++) {
      output = nn.layers[j].forward_func(nn.layers[j], output);
    }
    printf("tested: %i, %i, output: %f\n", (i32_t)X->data[i*2], (i32_t)X->data[(i*2)+1],*output->data );

  }

  return 0;
}
