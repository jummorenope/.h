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
    
double fibonacci(double x)
{
  if(n==1||n==2)
    return 1;
  
  else
    return fibonacci(n-1)+fibonacci(n-2);
}

long factorial(int x)
{
  if(n==0)
    return 1;
  
  else
    return n*factorial(n-1);
}

double forward(double (*f)(double), double x, double h)
{
  return (f(x+h)-f(x))/h;
}

double backward(double (*f)(double), double x, double h)

{
  return (f(x)-f(x-h))/h;
}

double halfp(double (*f)(double), double x, double h)

{
  return (f(x+h/2)-f(x-h/2))/h;
}

double richardson(double x, double h)
{
  return (4*halfp(x, h/2)-halfp(x, h))/3;
}



