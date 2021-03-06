#ifndef SQMS_IO_H
#define SQMS_IO_H
#if defined(__cplusplus)
extern "C" {
#endif


int SqMS_binarysearch(double arr[], int l, int r, double x);
  
int SqMS_biggest_lower_bound(const double arr[], const int n, const double x);

int powi(int base, unsigned int exp);

double SqMS_determinant(double *matrix, int dimension);

#if defined(__cplusplus)
}
#endif
#endif /* SQMS_MISC_H */