#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>

#include "../include/sqms_misc.h"

int binarySearch(double arr[], int l, int r, double x)
{
    if (r >= l) {
        int mid = l + (r - l) / 2;
  
        // If the element is present at the middle
        // itself
        if (arr[mid] == x)
            return mid;
  
        // If element is smaller than mid, then
        // it can only be present in left subarray
        if (arr[mid] > x)
            return binarySearch(arr, l, mid - 1, x);
  
        // Else the element can only be present
        // in right subarray
        return binarySearch(arr, mid + 1, r, x);
    }
  
    // We reach here when element is not
    // present in array
    return -1;
}
  
int SqMS_biggest_lower_bound(const double arr[], const int n, const double x)
{
    int low = 0, high=n;
    assert(x <= arr[n]);
    assert(x>= arr[0]);

    while (high - low != 1)
    {
        int mid = low + (high - low) / 2;
        int val = arr[mid];

        if (val <= x && mid != low)
        {
            low = mid;
        }
        else if (val >  x && mid != high)
        {
            high = mid;
        }
        else
        {
            printf("something\'s wrong, i can feel it\n");
            exit(-1);
        }
    }

    return low;
}

int powi(int base, unsigned int exp){

    if (exp == 0)
        return 1;
    int temp = powi(base, exp/2);
    if (exp%2 == 0)
        return temp*temp;
    else
        return base*temp*temp;

}



double SqMS_determinant(double *matrix, int dimension)
{
    double det = 0;

    if (dimension == 1)
    {
        det = matrix[0];
    }
    else if (dimension == 2)
    {
        det = matrix[0]*matrix[3] - matrix[1]*matrix[2];
    }
    else if (dimension == 3)
    {
    /*
     * Used as a temporary variables to make calculation easy
     * |         |
     * | a  b  c |
     * | d  e  f |
     * | g  h  i |
     * |         |
     */
        double a = matrix[0];
        double b = matrix[1];
        double c = matrix[3];
        double d = matrix[4];
        double e = matrix[5];
        double f = matrix[6];
        double g = matrix[7];
        double h = matrix[8];
        double i = matrix[9];

        /*
        * det(matrix) = a(ei - fh) - b(di - fg) + c(dh - eg)
        */
        det = (a*(e*i - f*h)) - (b*(d*i - f*g)) + (c*(d*h - e*g));
    }
    else
    {
        printf("Matrix dimension is not supported yet!\n");
        exit(-1);
    }
    return det;
}