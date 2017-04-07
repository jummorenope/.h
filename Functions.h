#include <iostream>
#include <cmath>

\ integrales \

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

\otras\

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

\derivadas\

double forward(double (*f)(double), double x, double h)
{
  return (f(x+h)-f(x))/h;
}

double backward(double (*f)(double), double x, double h)

{
  return (f(x)-f(x-h))/h;
}

double central(double (*f)(double), double x, double h)

{
  return (f(x+h/2)-f(x-h/2))/h;
}

double richardson(double x, double h)
{
  return (4*central(x, h/2)-central(x, h))/3;
}

\raices\

double fixed_point(double (*f)(double), double (*g)(double), double x, int MAX, double eps)
{
  double c, d;
  for (int ii=0; ii<=MAX; ii++)
    {
      c=g(x);
      if (std::fabs(f(c))<eps)
        {
          break;
        }
      d=1-c/x;
      std::cout<<ii<<"\t"<<c<<"\t"<<f(c)<<"\t"<<d<<std::endl;
      x=c;
    }
  return c;
}

double secante(double (*f)(double), double po, double pi, int MAX, double eps)
{
  double p;
  for (int a=0; a<=MAX; a++)
    {
      p=pi-(f(pi)*(po-pi)/(f(po)-f(pi)));
      if (std::fabs(f(p))<eps)
	{
	  break;
	}
      po=pi;
      pi=p;
    }
  return p;
}

double newton(double (*f)(double), double (*fp)(double), double po,int MAX,double eps)
{
  double p;
  for(int a=0;a<=MAX; a++)
    {
      p=po-(f(po))/(fp(po));
      
      if (std::fabs(f(p))<eps)
        {
          break;
        }
      po=p;
    }
  return p;
}
