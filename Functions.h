#include <iostream>
#include <cmath>
#include <vector>
//Integrales

double simpson(const double a, const double b, int n, double (*f)(double))
{
  if (n%2==1)
    {
      n=n+1;
    }
  
  double sum1=0, sum2=0, sum3=0, h;
  h=(b-a)/n;

  for (int ii=1; ii<=(n/2-1); ii++)
    {
      sum1=sum1+f(a+h*(2*ii));
    }

  for (int ii=1; ii<=(n/2); ii++)
    {
      sum2=sum2+f(a+h*(2*ii-1));
    }

  sum3=(2*sum1+4*sum2+f(a)+f(b))*(h/(3.0));

  return sum3;
}

double trapezoid_irregular(const std::vector<double>&x,const std::vector<double>&fx)
{
  double sum=0.0;
  const int n=x.size();
  for(int ii=0; ii<n-1; ii++)
    {
      sum += (x[ii+1]-x[ii])*(fx[ii]+fx[ii+1]);
    }
  sum*=0.5;
  return sum;
}

double trapezoid_regular(const double a, const double b, const int n, double (*f)(double))
{
  double sum=0.5*(f(a)+f(b));
  const double h=(b-a)/n;
  for(int ii=1; ii<n; ii++)
    {
      sum +=f(a+ii*h);
    }
  sum *= h;
  return sum;
}

double richardson_trapezoid(const double a, const double b, const int n, double (*f)(double))
{
  return (4*trapezoid_regular(a,b,2*n,f)-trapezoid_regular(a,b,n,f))/3;
}

double richardson_simpson(const double a, const double b, const int n, double (*f)(double))
{
  return (16*simpson(a,b,2*n,f)-simpson(a,b,n,f))/15;
}

//Funciones varias

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

//Derivadas

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

//Raices

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

//Matrices y arreglos

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
  std::cout<<std::endl;
}

void ortogonal(double m1[], int m, int n, double eps)
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
  if(c==1)
    std::cout<<"Es ortogonal"<<std::endl;
  if(c==0)
    std::cout<<"No es ortogonal"<<std::endl;
  
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

void transpond(double m1[], double m2[], int m , int n)
{
  for(int ii=0; ii<m; ii++)
    {
      for( int jj=0; jj<n; jj++)
	{
          m2[ii*n+jj]=m1[jj*m+ii];

	}
    }

}

void multiply(double m1[], double m2[], double m3[], int m, int n	\
)
{
  double ij=0;
  for (int ii=0; ii<m; ii++)
    {
      for (int jj=0; jj<m; jj++)
        {
          {
            ij=0;
            for (int i=0; i<n; i++)
              {
                ij=ij+m1[ii*n+i]*m2[i*m+jj];
              }
          }

          m3[ii*m+jj]=ij;
        }
      }
}


double determinant (double m1[], int m)
{
  double determinante=0;
  double m2[(m-1)*(m-1)];
  
  if(m==2)
    {
      determinante=m1[0]*m1[3]-m1[1]*m1[2];
    }
  
  else
    {
      for(int jj=0; jj<m; jj++)
	{
	  for (int i=0; i<m-1; i++)
	    {
	      for (int j=0; j<m-1; j++)
		{
		  if(j<jj)
		    m2[i*(m-1)+j]=m1[(i+1)*m+j];
		  
		  else
		    m2[i*(m-1)+j]=m1[(i+1)*m+(j+1)];		    
		}
	    }
	  determinante=determinante+std::pow(-1,jj)*m1[jj]*determinant(m2,m-1);
	}
    }
  
  return determinante;
}


void inverse (double m1[], double m2[], int m)
{  
  double determinante=determinant(m1, m);
  double k=1/determinante;
  if (determinante!=0)
    {
      for(int ii=0; ii<m; ii++)
        {
          for(int jj=0; jj<m; jj++)
            {
              m2[ii*m+jj]=m2[ii*m+jj]*k;
            }
        }
    }

  else
    std::cout<<"La matriz no tiene inversa"<<std::endl;

}


