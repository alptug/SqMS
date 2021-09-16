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