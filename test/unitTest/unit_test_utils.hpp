#include <stdio.h>
#include <stdlib.h>

template <typename T> bool TEST_EQUAL(T a, T b)
{
  return a == b;
}

bool TEST_EQUAL_INT_ARRAY(const int* a, const int* b, int length,
                          const char* name = NULL)
{
  for (int i = 0; i < length; i++) {
    if (a[i] != b[i]) {
      printf("TEST_EQUAL_INT_ARRAY for %s failed at index %d\n", name, i);
      return false;
    }
  }
  return true;
}

bool TEST_EQUAL_DOUBLE_ARRAY(const double* a, const double* b, int length,
                             const char* name = NULL)
{
  for (int i = 0; i < length; i++) {
    if (a[i] != b[i]) {
      printf("TEST_EQUAL_DOUBLE_ARRAY for %s failed at index %d\n", name, i);
      return false;
    }
  }
  return true;
}
