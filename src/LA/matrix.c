/*==========================================================*
 *
 * @author Gustaf Franzén :: https://github.com/BjorneEk;
 *
 *  basic matrix functions for linear algebra library
 *
 *==========================================================*/


#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>
#include "matrix.h"
#include "../UL/assert.h"
#include "../UL/IO.h"
#include "../UL/math.h"



/*==========================================================*
 *  allocates a new matrix of dimension cols x rows and
 *  setts all elements to zero.
 *==========================================================*/
LA_matrix * LA_init_empty(u32_t cols, u32_t rows)
{
  LA_matrix * res;

  assert(cols > 0 && rows > 0, "New matrix must be at least a 1 by 1");
  res = (LA_matrix *) malloc(sizeof(LA_matrix));

  assert(res != NULL, "Out of memory: (creating matrix).");

  res->cols  = cols;
  res->rows = rows;
  res->data   = (f64_t*) malloc(sizeof(f64_t) * cols * rows);

  assert(res->data != NULL, "Out of memory: (creating matrix data).");
  memset(res->data, 0.0, cols * rows * sizeof(f64_t));

  return res;
}

/*==========================================================*
 *  allocates a new matrix of dimension width x height and
 *  setts the elements to the suplied values
 *==========================================================*/
LA_matrix * LA_new_mat(u32_t cols, u32_t rows, ...)
{
  LA_matrix * m;
  u32_t size, i;
  va_list vals;

  size = cols * rows;
  va_start(vals, rows);
  m = LA_init_empty(cols, rows);

  for (i = 0; i < size; i++) m->data[i] = va_arg(vals, f64_t);

  va_end(vals);
  return m;
}
/*==========================================================*
 *  creates a new copy of the input matrix
 *==========================================================*/
LA_matrix * LA_copy(const LA_matrix * const m)
{
  LA_matrix * res;
  res = (LA_matrix *) malloc(sizeof(LA_matrix));
  assert(res != NULL, "Out of memory: (copying matrix).");
  res->cols  = m->cols;
  res->rows = m->rows;
  res->data   = (f64_t*) malloc(sizeof(f64_t) * m->cols * m->rows);

  assert(res->data != NULL, "Out of memory: (copying matrix data).");
  memcpy(res->data, m->data, m->rows*m->cols*sizeof(f64_t));
  return res;
}

/*==========================================================*
 *  copy m into res;
 *==========================================================*/
void LA_copy_r(LA_matrix * res, const LA_matrix * const m)
{
  assert(res->cols == m->cols && res->rows == m->rows,
    "incorrect dimensions: (copying matrix).");
  res->cols  = m->cols;
  res->rows = m->rows;
  memcpy(res->data, m->data, m->rows*m->cols*sizeof(f64_t));
}

/*==========================================================*
 *  frees the memory allocated to the matrix m
 *==========================================================*/
void LA_free(LA_matrix * m)
{
  if (m != NULL) {
    if (m->data != NULL) { free(m->data); m->data = NULL;}
    free(m); m = NULL;
  }
}

/*==========================================================*
 *  prints m to the console
 *==========================================================*/
void LA_printmat(const LA_matrix * m)
{
  i32_t i, j;
  f64_t * ptr = m->data;
  printf("│");
  for (i = 0; i < m->rows; i++) {
    for (j = 0; j < m->cols; j++)
      printf("%5.2f ", *(ptr++));
    printf("│\n");
    if(i<m->rows-1) printf("│");
  }
}

/*==========================================================*
 *  prints m to the console verbously
 *==========================================================*/
void LA_debugmat(const LA_matrix * m)
{
  i32_t i, j;
  f64_t * ptr = m->data;
  printf("%d %d\n", m->cols, m->rows);
  for (i = 0; i < m->rows; i++) {
    for (j = 0; j < m->cols; j++)
      printf("%f ", *(ptr++));
    printf("\n");
  }
}

/*==========================================================*
 *  reads a matrix from a file according to the following
 *  format: cols space rows newline then the matrix
 *  with cols separated by spaces and rows by newlines
 *==========================================================*/
LA_matrix * LA_readmat(const char * filename)
{
  FILE * fp;
  LA_matrix * res;
  f32_t val;
  i32_t cols, rows, i, elements;
  i32_t scan_test;
  f64_t * ptr;

  asserts((fp = fopen(filename, "r")) != NULL, "LA_readmat Cannot open", filename);

  scan_test = fscanf(fp, "%d", &cols);
  assert(scan_test != EOF, "Failed to read from file.");

  scan_test = fscanf(fp, "%d", &rows);
  assert(scan_test != EOF, "Failed to read from file.");

  res = LA_init_empty(cols, rows);
  elements = cols * rows;

  ptr = res->data;
  for (i = 0; i < elements; i++) {
    scan_test = fscanf(fp, "%f", &val);
    assert(scan_test != EOF, "Failed to read from file. File is incomplete.");
    *(ptr++) = val;
  }
  fclose(fp);
  return res;
}



/*==========================================================*
 *  reads a array of matricies from a file
 *==========================================================*/
LA_matrix * LA_read_matricies(const char * filename, i32_t * len)
{
  LA_matrix *res;
  char      *buffer, *p;
  i32_t        i, j, count;
  f64_t    *ptr;

  count = 0;
  read_sfile(filename, buffer);
  p = buffer;
  inc:
  if (*p == '[') count = atoi(p+1);
  else {p++; goto inc;}

  assert(count != 0, "could not read matrix count from file");

  res = (LA_matrix *)malloc(count * sizeof(LA_matrix));
  assert(res != NULL, "Out of memory");
  for(; *p!='\n'; p++);

  i = 0;
  while(i <= count && *p != '\0') {
    switch (*p) {
      case '(':
        res[i].rows = atoi(p+1);
        left_bracket: p++;
        if(*p == ' ') {
          res[i].cols = atoi(p+1);
          res[i].data = (f64_t *)malloc(res[i].rows * res[i].cols * sizeof(f64_t));
          assert(res[i].data != NULL, "Out of memory");
        }
        else goto left_bracket;
        break;
      case '{':
        j = 0;
        while (j < res[i].rows * res[i].cols && *p != '}') {
          if((*p == ' ' || *p == '\n') && (*(p+1) != ' ' && *(p+1) != '\n')) {
            if (*p == '}') break;
            else {
              res[i].data[j] = atof(p+1);
              j++;
            }
          }
          p++;
        }
        i++; break;
      default: p++; break;
    }
  }
  *len = count;
  return res;
}

/*==========================================================*
 *  writes a matrix to a file
 *==========================================================*/
void LA_writemat(const LA_matrix * m, const char * filename)
{
  FILE * fp;
  i32_t i, j;
  f64_t * ptr = m->data;

  asserts((fp = fopen(filename, "w")) != NULL, "LA_writemat Cannot open", filename);

  fprintf(fp, "%d %d\n", m->cols, m->rows);
  for (i = 0; i < m->rows; i++) {
    for (j = 0; j < m->cols; j++)
      fprintf(fp, " %2.5f", *(ptr++));
    fprintf(fp, "\n");
  }
  fclose(fp);
}

/*==========================================================*
 *  sets the diagonal in m to d
 *==========================================================*/
void LA_set_diagonal(LA_matrix *m, f64_t d)
{
  i32_t sz, i;
  assert(m->cols == m->rows, "Matrix must be square (set diagonal)");
  sz = m->rows*m->cols;
  for(i=0; i<sz; m->data[i]=d, i+=m->cols+1);
}

/*==========================================================*
 *  creates a nxn identity matrix
 *==========================================================*/
LA_matrix * LA_identity_mat(u32_t n)
{
  LA_matrix * res;
  i32_t i, size;
  res = LA_init_empty(n, n);
  size = n*n;
  for(i=0; i<size; res->data[i]=1.0, i+=n+1);
  return res;
}

/*==========================================================*
 *  creates a nxm random matrix
 *==========================================================*/
LA_matrix * LA_random(u32_t cols, u32_t rows, f64_t min, f64_t max)
{
  LA_matrix * res;
  i32_t i, size;
  f64_t *p;
  res = LA_init_empty(cols, rows);
  size = cols*rows;
  for(i=0, p=res->data; i<size; *(p++)=randf(min,max), i++);
  return res;
}

/*==========================================================*
 *  creates a transpose the matrix m
 *==========================================================*/
LA_matrix * LA_transpose(const LA_matrix * m)
{
  LA_matrix * res;
  i32_t i, j;
  f64_t * mp, * rp;

  res = LA_init_empty(m->rows, m->cols);
  rp = res->data;

  for (i=0; i<m->cols;i++) {
    mp = &m->data[i];
    for (j=0; j<m->rows; j++, rp++, mp+=m->cols) *rp = *mp;
  }
  return res;
}

/*==========================================================*
 *  multiplies each element in m with s
 *==========================================================*/
void LA_scale_r(LA_matrix * m, f64_t s)
{
  i32_t sz, i;
  f64_t * p;
  sz = m->rows*m->cols;
  for(i=0, p=m->data; i<sz; *(p++)*=s, i++);
}
/*==========================================================*
 *  creates a new matrix with the values of m scaled by s
 *==========================================================*/
LA_matrix * LA_scale(const LA_matrix * m, f64_t s)
{
  LA_matrix * res;
  res = LA_copy(m);
  LA_scale_r(res, s);
  return res;
}

/*==========================================================*
 *  add the matrix b to the matrix a
 *==========================================================*/
void LA_add_r(LA_matrix * a, const LA_matrix * b)
{
  i32_t sz, i, j;
  f64_t *pa = a->data, *pb = b->data;
  assert((a->rows == b->rows) && (a->cols == b->cols),
    "matrix A and B must be same size to perform addition");
  sz = a->rows*a->cols;
  for(i=0; i<sz; *(pa++)+=*(pb++), i++);
}
/*==========================================================*
 *  creates a new matrix that is equal to a + b
 *==========================================================*/
LA_matrix * LA_add(const LA_matrix * a, const LA_matrix * b)
{
  LA_matrix * res;
  res = LA_copy(a);
  LA_add_r(res, b);
  return res;
}

/*==========================================================*
 *  subtract the matrix b to the matrix a
 *==========================================================*/
void LA_sub_r(LA_matrix * a, const LA_matrix * b)
{
  i32_t sz, i, j;
  f64_t *pa = a->data, *pb = b->data;
  assert((a->rows == b->rows) && (a->cols == b->cols),
    "matrix A and B must be same size to perform subtraction");
  sz = a->rows*a->cols;
  for(i=0; i<sz; *(pa++)-=*(pb++), i++);
}
/*==========================================================*
 *  creates a new matrix that is equal to a - b
 *==========================================================*/
LA_matrix * LA_sub(const LA_matrix * a, const LA_matrix * b)
{
  LA_matrix * res;
  res = LA_copy(a);
  LA_sub_r(res, b);
  return res;
}

/*==========================================================*
 *  element wise multiplication of the matrix b to the matrix a
 *==========================================================*/
void        LA_ewise_mult_r(LA_matrix * a, const LA_matrix * b)
{
  i32_t sz, i, j;
  f64_t *pa = a->data, *pb = b->data;
  assert((a->rows == b->rows) && (a->cols == b->cols),
    "matrix A and B must be same size to perform subtraction");
  sz = a->rows*a->cols;
  for(i=0; i<sz; *(pa++)*=*(pb++), i++);
}
/*==========================================================*
 *  creates a new matrix that is equal to element wise
 *  multiplication of a and b
 *==========================================================*/
LA_matrix * LA_ewise_mult(const LA_matrix * a, const LA_matrix * b)
{
  LA_matrix * res;
  res = LA_copy(a);
  LA_ewise_mult_r(res, b);
  return res;
}



void LA_matmult_add_rec(const f64_t *A, const f64_t *B, f64_t *C, i32_t m,   i32_t n,
                             i32_t p,        i32_t fdA,          i32_t fdB, i32_t fdC)
{
  i32_t i, j, k, m2, n2, p2;
  f64_t sum;
  if (m+n+p <= 48) { /* <= 16x16 matrices "on average" */
    for (i = 0; i < m; ++i)
      for (k = 0; k < p; ++k) {
        sum = 0;
        for (j = 0; j < n; ++j) sum += A[i*fdA +j] * B[j*fdB + k];
        C[i*fdC + k] += sum;
      }
  } else { /* divide and conquer */
    m2 = m/2, n2 = n/2, p2 = p/2;
    LA_matmult_add_rec(A, B, C, m2, n2, p2, fdA, fdB, fdC);
    LA_matmult_add_rec(A+n2, B+n2*fdB, C, m2, n-n2, p2, fdA, fdB, fdC);
    LA_matmult_add_rec(A, B+p2, C+p2, m2, n2, p-p2, fdA, fdB, fdC);
    LA_matmult_add_rec(A+n2, B+p2+n2*fdB, C, m2, n-n2, p-p2, fdA, fdB, fdC);
    LA_matmult_add_rec(A+m2*fdA, B, C+m2*fdC, m-m2, n2, p2, fdA, fdB, fdC);
    LA_matmult_add_rec(A+m2*fdA+n2, B+n2*fdB, C+m2*fdC, m-m2, n-n2, p2, fdA, fdB, fdC);
    LA_matmult_add_rec(A+m2*fdA, B+p2, C+m2*fdC+p2, m-m2, n2, p-p2, fdA, fdB, fdC);
    LA_matmult_add_rec(A+m2*fdA+n2, B+p2+n2*fdB, C+m2*fdC, m-m2, n-n2, p-p2, fdA, fdB, fdC);
  }
}

/*==========================================================*
 *  creates a new matrix that is equal to a * b
 *==========================================================*/
LA_matrix * LA_mult(const LA_matrix * a, const LA_matrix * b)
{
  LA_matrix * res;

  if(b->cols != b->rows || a->cols != a->rows ||
    (a->rows <= MIN_MAT_THRESHOLD || a->cols <= MIN_MAT_THRESHOLD ||
    b->rows <= MIN_MAT_THRESHOLD || b->cols <= MIN_MAT_THRESHOLD)) return LA_mult_naive(a, b);

  assert(a->cols == b->rows, "Matrices have incorrect dimensions. a->cols != b->rows");

  res = LA_init_empty(b->cols, a->rows);
  LA_matmult_add_rec(a->data, b->data, res->data, a->rows, b->rows,
    res->rows, b->rows, res->rows, res->rows);
  return res;
}

LA_matrix * LA_mult_naive(const LA_matrix * a, const LA_matrix * b)
{
    i32_t i, j, k;
    f64_t *pa, *pb, *pr;
    LA_matrix * res;

    assert(a->cols == b->rows, "Matrices have incorrect dimensions. a->cols != b->rows");
    res = LA_init_empty(b->cols, a->rows);

    pr = res->data;
    // Multiplying matrix firstMatrix and secondMatrix and storing in array mult.*/
    for(i=0; i<a->rows; i++) {
      for(j=0; j<b->cols; j++, pr++) {
        pa = &a->data[i * a->cols];
        pb = &b->data[j];
        for(k=0; k<a->cols; k++, pa++, pb+=b->cols)
          *pr += *pa * *pb;
      }
    }
    return res;
}

/*==========================================================*
 *  returns true is a is identical to b
 *==========================================================*/
bool LA_equals(const LA_matrix * a, const LA_matrix * b)
{
  i32_t i, sz;
  f64_t *pa, *pb;

  if (a->cols != b->cols || a->rows != b->rows) return false;

  sz = a->rows*a->cols;
  pa = a->data; pb = b->data;
  for(i = 0; i <sz; i++,pb++,pa++) if (*pa != *pb) return false;
  return true;
}

/*==========================================================*
 *  swap row _p with row _q in the matrix m
 *==========================================================*/
void LA_rswap(LA_matrix * m, i32_t _p, i32_t _q)
{
  i32_t i;
  f64_t tmp;
  f64_t* p, * q;

  assert(m->rows >= 2, "Matrix must have at least two rows to swap.");
  assert(_p < m->rows && _q < m->rows,
         "Values p and q must be less than the rows of the matrix.");

  if (_p == _q) return;

  p = m->data + (_p * m->cols);
  q = m->data + (_q * m->cols);

  for (i = 0; i < m->cols; i++, p++, q++) {
    tmp = *p;
    *p  = *q;
    *q  = tmp;
  }
}

/*==========================================================*
 *  swap column _p with column _q in the matrix m
 *==========================================================*/
void LA_cswap(LA_matrix * m, i32_t _p, i32_t _q)
{
  i32_t i;
  f64_t tmp;
  f64_t * p, * q;

  assert(m->rows >= 2, "Matrix must have at least two rows to swap.");
  assert(_p < m->rows && _q < m->rows,
         "Values p and q must be less than the rows of the matrix.");

  if (p == q) return;

  p = &m->data[_p];
  q = &m->data[_q];

  for (i = 0; i < m->rows; i++, p+=m->cols, q+=m->cols) {
    tmp = *p;
    *p  = *q;
    *q  = tmp;
  }
}

/*==========================================================*
 *  returns the row r of a matrix m
 *==========================================================*/
LA_matrix * LA_rowat(const LA_matrix * m, u32_t r)
{
  LA_matrix * res;
  assert(r < m->rows, "Row idx out of bounds");

  res = LA_init_empty(m->cols, 1);
  memcpy(res->data, m->data + r * m->cols, m->cols*sizeof(f64_t));
  return res;
}

/*==========================================================*
 *  returns the col c of a matrix m
 *==========================================================*/
LA_matrix * LA_colat(const LA_matrix * m, u32_t c)
{
  LA_matrix * res;
  i32_t i;
  f64_t *p, *pm;
  assert(c < m->cols, "Col idx out of bounds");

  res = LA_init_empty(1, m->rows);
  p = res->data;
  pm = &m->data[c];
  for(i=0; i<m->rows; i++, p++, pm += m->cols) *p = *pm;
  return res;
}

/*==========================================================*
 *  wultiply al elements in row r with s
 *==========================================================*/
void LA_mult_row_r(LA_matrix * m, u32_t r, f64_t s)
{
  f64_t * p;
  i32_t i;
  assert(r < m->rows, "Row idx out of bounds");
  p = &m->data[r*m->cols];
  for(i=0; i<m->cols; *(p++)*=s,i++);
}

/*==========================================================*
 *  create a new copy of a where row r has been
 *  multiplied by s
 *==========================================================*/
LA_matrix * LA_mult_row(const LA_matrix * m, u32_t r, f64_t s)
{
  LA_matrix * res;
  res = LA_copy(m);
  LA_mult_row_r(res, r, s);
  return res;
}

/*==========================================================*
 *  wultiply al elements in col c with s
 *==========================================================*/
void LA_mult_col_r(LA_matrix * m, u32_t c, f64_t s)
{
  f64_t * p;
  i32_t i;
  assert(c < m->cols, "Col idx out of bounds");
  p = &m->data[c];
  for(i=0; i<m->rows; *p*=s,i++,p+=m->cols);
}

/*==========================================================*
 *  create a new copy of a where col c has been
 *  multiplied by s
 *==========================================================*/
LA_matrix * LA_mult_col(const LA_matrix * m, u32_t c, f64_t s)
{
  LA_matrix * res;
  res = LA_copy(m);
  LA_mult_col_r(res, c, s);
  return res;
}

/*==========================================================*
 * add row r multiplied with the scalar s to the row d
 *==========================================================*/
void LA_row_add_row_r(LA_matrix * m, u32_t d, u32_t r, f64_t s)
{
  f64_t *pd, *pr;
  i32_t i;
  assert(r < m->rows || d < m->rows, "Row idx out of bounds");
  pd = &m->data[d*m->cols]; pr = &m->data[r*m->cols];
  for(i=0; i<m->cols; i++, *pd++ += (s*(*pr++)));
}

/*==========================================================*
 * create a new copy of m where the row r multiplied with
 * the scalar s have been addet to the row d
 *==========================================================*/
LA_matrix * LA_row_add_row(LA_matrix * m, u32_t d, u32_t r, f64_t s)
{
  LA_matrix * res;
  res = LA_copy(m);
  LA_row_add_row_r(res, s, r, s);
  return res;
}

/*==========================================================*
 * calculates dot product of two vectors
 *==========================================================*/
f64_t LA_dot(const LA_matrix * a, const LA_matrix * b)
{
  i32_t i,sz;
  f64_t res, *ap, *bp;
  res = 0;
  assert((a->cols == 1 || a->rows == 1) && (b->cols == 1 || b->rows == 1),
    "Wrong dimensions dot");
  sz = a->cols*a->rows;
  assert(sz == (b->cols * b->rows), "Wrong dimensions dot");

  ap = a->data; bp = b->data;
  for(i=0;i<sz;res+=*ap++**bp++,i++);
  return res;
}
/*==========================================================*
 * calculates dot product of two columns in two matricies
 *==========================================================*/
f64_t LA_dot_col(const LA_matrix * a, u32_t ac, const LA_matrix * b, u32_t bc)
{
  i32_t i;
  f64_t res, *ap, *bp;

  assert(a->rows == b->rows, "Wrong dimensions (dot col)");
  assert((ac < a->cols) && (bc < b->cols), "Index error (dot col)");
  res = 0.0;

  ap = a->data; bp = b->data;

  for(i=0;i<a->rows;res+=*ap**bp, ap+=a->cols,bp+=b->cols,i++);
  return res;
}

/*==========================================================*
 * calculates absolute value (length) of a vector
 *==========================================================*/
f64_t LA_abs(const LA_matrix * a)
{
  i32_t i,sz;
  f64_t res, *ap;
  res = 0;
  assert((a->cols == 1 || a->rows == 1), "Wrong dimensions dot");
  sz = a->cols*a->rows;
  ap = a->data;
  for(i=0;i<sz;res+=*ap**ap,i++,ap++);
  return sqrt(res);
}
/*==========================================================*
 * calculates absolute value (length) of a column
 *==========================================================*/
f64_t LA_abs_col(const LA_matrix * a, u32_t c)
{
  i32_t i;
  f64_t res, *ap;

  assert(c < a->cols, "Index error (dot col)");
  res = 0.0;

  ap = a->data;

  for(i=0;i<a->rows;res+=*ap**ap, ap+=a->cols,i++);
  return sqrt(res);
}

/*==========================================================*
 * calculates absolute value (length) of each column
 * in a matrix and returns a 1 x 3 matrix with the result
 *==========================================================*/
LA_matrix * LA_abs_mat(const LA_matrix * a)
{
  LA_matrix * res;
  i32_t i, j;
  f64_t *ap, *rp;


  res = LA_init_empty(a->cols, 1);
  ap = a->data;
  rp = res->data;


  for(i=0;i<a->cols;i++, rp++) {
    ap = a->data + i;
    for(j=0;j<a->rows;j++, ap+=a->cols)
      *rp += *ap**ap;
    *rp = sqrt(*rp);
  }
  return res;
}

/*==========================================================*
 * applies f to each element in m
 *==========================================================*/
void LA_foreach(LA_matrix * m, f64_t (*f)(f64_t))
{
  i32_t i;
  for(i=0; i<m->cols*m->rows; i++) m->data[i]=f(m->data[i]);
}

/*==========================================================*
 * returns the mean value of m
 *==========================================================*/
f64_t LA_mean(const LA_matrix * m)
{
  i32_t i;
  f64_t res;

  for(i=0; i<m->cols*m->rows; i++) res+=m->data[i];
  return res / (m->cols*m->rows);
}
