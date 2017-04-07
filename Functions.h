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

long fibonacci(int x)
{
  if(x==1||x==2)
    return 1;
  
  else
    return fibonacci(x-1)+fibonacci(x-2);
}

long factorial(int x)
{
  if(x==0)
    return 1;
  
  else
    return x*factorial(x-1);
}

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

double richardson(double (*f)(double), double x, double h)
{
  return (4*central(f, x, h/2)-central(f, x, h))/3;
}

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

int primo(long int x)
{
  int c=1;
  double d=sqrt(x);
  for(int ii=2; ii<=d; ii++)
    {
      if(x%ii==0)
	c=0;
    }
  return c;
}

void print(double ma[], int m, int n)
{
  for(int ii=0; ii<m; ii++)
    {
      for( int jj=0; jj<n; jj++)
	{
	  std::cout<<ma[ii*n+jj]<<"\t";
	}
      std::cout<<std::endl;
    }
}

double multiply(double m1[], double m2[], int ii, int jj, int N, int M)
{
  double ij=0;
  for (int i=0; i<N; i++)
    {
      ij=ij+m1[ii*N+i]*m2[i*M+jj];
    }
  return ij;
}

int ortogonal(double m1[], int m, int n, double eps)
{
  int c=1;
  for (int ii=0; ii<m; ii++)
    {
      for (int jj=0; jj<n; jj++)
        {
          if(ii==jj)
            {
              if(std::fabs(m1[ii*n+jj]-1)>eps)
                {
                  c=0;
                  break;
        	}
            }

          if(ii!=jj)
            {
              if(std::fabs(m1[ii*n+jj])>eps)
		{
                  c=0;
                  break;
		}
            }
	}
    }
  return c;
}

double media(double a[], int n)
{
  double sum=0;
  for(int ii=0; ii<n; ii++)
    {
      sum=sum+a[ii];
    }
  sum=sum/n;
  return sum;
}

double desviacion(double a[], int n, double media)
{
  double desviacion=0;
  for(int ii=0; ii<n; ii++)
    {
      desviacion=desviacion+(media-a[ii])*(media-a[ii]);
    }
  desviacion=desviacion/(n-1);
  desviacion=std::sqrt(desviacion);
  return desviacion;
}

double maximo(double a[], int n)
{
  double maximo=a[0];
  for(int ii=1; ii<n; ii++)
    {
      if(a[ii]>maximo)
        {
          maximo=a[ii];
        }
    }
  return maximo;
}

double minimo(double a[], int n)
{
  double minimo=a[0];
  for(int ii=1; ii<n; ii++)
    {
      if(a[ii]<minimo)
        minimo=a[ii];
    }
  return minimo;
}

