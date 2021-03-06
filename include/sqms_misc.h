#ifndef SQMS_MISC_H
#define SQMS_MISC_H
#if defined(__cplusplus)
extern "C" {
#endif


int SqMS_binarysearch(double arr[], int l, int r, double x);
  
int SqMS_biggest_lower_bound(const double arr[], const int n, const double x);

int powi(int base, unsigned int exp);

double SqMS_determinant(double *matrix, int dimension);

void print_ghostbusters();

int compare_ints(const void* a, const void* b);

void SqMS_sorted_set_difference(const int* difference_of_this, const int* from_this, int* checkboard, int num_el);
#if defined(__cplusplus)
}
#endif
#endif /* SQMS_MISC_H */