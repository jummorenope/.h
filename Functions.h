#include <iostream>
#include <cmath>

double simpson(double (*f)(double), double b, double a, double n)
{
  double sum1=0, sum2=0, sum3=0, c;
  c=b-a;

  for (int ii=1; ii<(n/2-1); ii++)
    {
      sum1=sum1+2*f(a+c*(2*ii)/n);
    }

  for (int ii=1; ii<(n/2); ii++)
    {
      sum2=sum2+4*f(a+c*(2*ii-1)/n);
    }

  sum3=(sum1+sum2+f(a)+f(b))*(c/(3.0*n));

  return sum3;
}
    
