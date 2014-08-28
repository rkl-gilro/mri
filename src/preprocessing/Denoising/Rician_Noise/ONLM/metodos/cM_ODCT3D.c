/***************************************************************************                  
/* Jose V. Manjon - jmanjon@fis.upv.es                                     */
/* Universidad Politecnica de Valencia, Spain                              */
/* Pierrick Coupe - pierrick.coupe@gmail.com                               */
/* Brain Imaging Center, Montreal Neurological Institute.                  */
/* Mc Gill University                                                      */
/*                                                                         */
/* Copyright (C) 2010 Jose V. Manjon and Pierrick Coupe                    */

/***************************************************************************
* New Methods for MRI Denoising based on Sparseness and Self-Similarity    *
* Jos� V. Manj�n, Pierrick Coup�, Antonio Buades,D. Louis Collins 
* and Montserrat Robles
***************************************************************************/

#include "math.h"
#include "mex.h"
#include <stdlib.h>
#include "matrix.h"
#ifndef PI
#define PI 3.14159265358979323846
#endif
#define B 4
#define BB 16

/* Multithreading stuff*/
#ifdef _WIN32
#include <windows.h>
#include <process.h>
#else
#include <pthread.h>
#endif

typedef struct{
    int rows;
    int cols;
    int slices;
    double * in_image;
    double * out_image;
    double * out_image2;
    double * acu_image;
    int ini;
    int fin;    
    double sigma;
    double * mibias;
}myargument;

double W1,W2; 
double *tabla1,*tabla2; /* tablas de cosenos*/
int rician;

/*Returns the modified Bessel function I0(x) for any real x.*/
double Bessel0(double x)
{
  double ax,ans,a;
  double y; 
  if ((ax=fabs(x)) < 3.75) 
  { 
    y=x/3.75;
    y*=y;
    ans=1.0+y*(3.5156229+y*(3.0899424+y*(1.2067492+y*(0.2659732+y*(0.360768e-1+y*0.45813e-2)))));
  } 
  else 
  {
    y=3.75/ax;
    ans=(exp(ax)/sqrt(ax));
    a=y*(0.916281e-2+y*(-0.2057706e-1+y*(0.2635537e-1+y*(-0.1647633e-1+y*0.392377e-2))));
    ans=ans*(0.39894228 + y*(0.1328592e-1 +y*(0.225319e-2+y*(-0.157565e-2+a))));    
  }
  return ans;
}

/*Returns the modified Bessel function I1(x) for any real x.*/
double Bessel1(double x)
{
  double ax,ans;
  double y; 
  if ((ax=fabs(x)) < 3.75)
  { 
    y=x/3.75;
    y*=y;
    ans=ax*(0.5+y*(0.87890594+y*(0.51498869+y*(0.15084934+y*(0.2658733e-1+y*(0.301532e-2+y*0.32411e-3))))));
  } 
  else 
  {
    y=3.75/ax;
    ans=0.2282967e-1+y*(-0.2895312e-1+y*(0.1787654e-1-y*0.420059e-2));
    ans=0.39894228+y*(-0.3988024e-1+y*(-0.362018e-2+y*(0.163801e-2+y*(-0.1031555e-1+y*ans))));
    ans *= (exp(ax)/sqrt(ax));
  }
  return x < 0.0 ? -ans : ans;
}

void RicianBias(double level,int max,double * mibias,int N)
{
    int cz,j;
    double c,z,i;          
    
    for(i=0;i<max*N;i=i+1) 
    {
        mibias[(int)i]=0;         
    }
    
    for(i=0;i<max;i=i+0.1)
    {      
      c=(i*i)/(level*level);
      z=(PI/2)*exp(-(c/2))*((1+c/2)*Bessel0(c/4)+(c/2)*Bessel1(c/4))*((1+c/2)*Bessel0(c/4)+(c/2)*Bessel1(c/4));
      z=(sqrt(z*(level*level)));
      cz=(int)z;
      if(cz<=0) cz=(int)i;
      if(cz>=0 && cz<max)
      { 
          if(mibias[cz]==0) 
          {  
              for(j=0;j<N;j++) 
              {
                  mibias[j*max+cz]+=i;
              }              
          }
      }
    }
}

static void dct(double *data,double * vol)
{
    int r,k;
    double spec;
	extern double *tabla1;

    /* do the stuff          */   
    for(r=0;r<B;r++) 
    {
	  spec= data[0]*tabla1[r*B] + data[1]*tabla1[r*B+1] + data[2]*tabla1[r*B+2] + data[3]*tabla1[r*B+3];
      vol[r] = W2*spec;
      if(r==0) vol[r] = W1*spec;                          
    }
}

static void dct3(double *data,double * vol)
{
    int r1,r2,r3;
    double v2[4];
    double v3[4];    
    double s2[4];
    double s3[4];
    double tmp[64];    
      
    /* do the stuff     */
    for(r3=0;r3<B;r3++) 
    for(r2=0;r2<B;r2++) 
    {      
      dct(&data[r3*BB+r2*B],&tmp[r3*BB+r2*B]);           
    }
    
    for(r1=0;r1<B;r1++) 
    for(r3=0;r3<B;r3++) 
    {      
      v2[0]=tmp[r3*BB+r1]; 
      v2[1]=tmp[r3*BB+B+r1]; 
      v2[2]=tmp[r3*BB+2*B+r1]; 
      v2[3]=tmp[r3*BB+3*B+r1]; 
      dct(v2,s2);
      tmp[r3*BB+r1]=s2[0]; 
      tmp[r3*BB+B+r1]=s2[1]; 
      tmp[r3*BB+2*B+r1]=s2[2]; 
      tmp[r3*BB+3*B+r1]=s2[3];      
    }
    
    for(r1=0;r1<B;r1++) 
    for(r2=0;r2<B;r2++) 
    {      
      v3[0]=tmp[r2*B+r1];                      
      v3[1]=tmp[BB+r2*B+r1];                     
      v3[2]=tmp[2*BB+r2*B+r1];                     
      v3[3]=tmp[3*BB+r2*B+r1];                     
      dct(v3,s3);
      vol[r2*B+r1]=s3[0];                      
      vol[BB+r2*B+r1]=s3[1];                     
      vol[2*BB+r2*B+r1]=s3[2];                     
      vol[3*BB+r2*B+r1]=s3[3];                                
    }
}

static void idct(double *data,double * vol)
{
    int r,k;
    double spec;
	extern double *tabla2;

    /* do the stuff         */
    for(r=0;r<B;r++) 
    {
      vol[r] = W1*data[0]*tabla2[r*B] + W2*data[1]*tabla2[r*B+1] + W2*data[2]*tabla2[r*B+2] + W2*data[3]*tabla2[r*B+3];                   
    }
}

static void idct3(double *data,double * vol)
{
    int r1,r2,r3;    
    double v2[4];
    double v3[4];    
    double s2[4];
    double s3[4];
    double tmp[64];    
      
    /* do the stuff     */
    for(r3=0;r3<B;r3++) 
    for(r2=0;r2<B;r2++) 
    {      
      idct(&data[r3*BB+r2*B],&tmp[r3*BB+r2*B]);     
    }
    
    for(r1=0;r1<B;r1++) 
    for(r3=0;r3<B;r3++) 
    {      
      v2[0]=tmp[r3*BB+r1]; 
      v2[1]=tmp[r3*BB+B+r1]; 
      v2[2]=tmp[r3*BB+2*B+r1]; 
      v2[3]=tmp[r3*BB+3*B+r1]; 
      idct(v2,s2);
      tmp[r3*BB+r1]=s2[0]; 
      tmp[r3*BB+B+r1]=s2[1]; 
      tmp[r3*BB+2*B+r1]=s2[2]; 
      tmp[r3*BB+3*B+r1]=s2[3];      
    }
    
    for(r1=0;r1<B;r1++) 
    for(r2=0;r2<B;r2++) 
    {      
      v3[0]=tmp[r2*B+r1];                      
      v3[1]=tmp[BB+r2*B+r1];                     
      v3[2]=tmp[2*BB+r2*B+r1];                     
      v3[3]=tmp[3*BB+r2*B+r1];     
      idct(v3,s3);
      vol[r2*B+r1]=s3[0];                      
      vol[BB+r2*B+r1]=s3[1];                     
      vol[2*BB+r2*B+r1]=s3[2];                     
      vol[3*BB+r2*B+r1]=s3[3];             
    }      
}

void* ThreadFunc( void* pArguments )
{
    double *ima,*fima,*acu,sigma,T,z,fac,max;
    int ind,ii,jj,kk,i,j,k,ini,fin,rows,cols,slices,p,p1;   
    double b[64]; 
    double cb[64]; 
 
    myargument arg;
    arg=*(myargument *) pArguments;

    rows=arg.rows;    
    cols=arg.cols;
    slices=arg.slices;
    ini=arg.ini;    
    fin=arg.fin;
    ima=arg.in_image;
    fima=arg.out_image;    
    acu=arg.acu_image;    
    sigma=arg.sigma; 
            
    T=2.7*sigma;

    for(k=ini;k<=fin;k=k+1)
    {  
	for(j=0;j<=rows-B;j=j+1)
	{
		for(i=0;i<=cols-B;i=i+1)
		{	            
            max=0;
            for(kk=0;kk<B;kk++)
                for(jj=0;jj<B;jj++)	
                    for(ii=0;ii<B;ii++)
                    {	
                      p=kk*B*B+jj*B+ii;             
                      b[p] = ima[(k+kk)*rows*cols+(j+jj)*cols+(i+ii)];                     
                      if(b[p]>max) max=b[p];
                    }         
            if(max==0) continue;
                                              
            dct3(b,cb);
                       
            ind=0;
            for(kk=0;kk<B;kk++)
                for(jj=0;jj<B;jj++)	
                    for(ii=0;ii<B;ii++)
                    {
                        p=kk*B*B+jj*B+ii;                
                        if(abs(cb[p])<T) cb[p]=0; /* hardthreshold(cb,T);       */
                        else ind++;                
                    }
            z=1.0/(ind+1); 

            idct3(cb,b);
                                
            for(kk=0;kk<B;kk++)
                for(jj=0;jj<B;jj++)	
                    for(ii=0;ii<B;ii++)
                    {	
                      p1=(k+kk)*rows*cols+(j+jj)*cols+(i+ii);                

                      /* agregation                  */
                      fima[p1]+=z*b[kk*B*B+jj*B+ii];  
                      acu[p1]+=z;     

                    }  
            
		}
	} 
    }   
  
    #ifdef _WIN32
       _endthreadex(0);    
    #else
       pthread_exit(0);    
    #endif

	return 0;
} 

void* ThreadFunc2( void* pArguments )
{
    double *ima,*fima,*fima2,*acu,*mibias,sigma,T,z,fac,val,v1,max;
    int ind,ii,jj,kk,i,j,k,ini,fin,rows,cols,slices,p,p1,iii;
    double b[64]; 
    double cb[64]; 
    double b2[64]; 
    double cb2[64]; 
    extern int rician;    
 
    myargument arg;
    arg=*(myargument *) pArguments;

    rows=arg.rows;    
    cols=arg.cols;
    slices=arg.slices;
    ini=arg.ini;    
    fin=arg.fin;
    ima=arg.in_image;
    fima=arg.out_image;
    fima2=arg.out_image2;
    acu=arg.acu_image;    
    sigma=arg.sigma; 
    mibias=arg.mibias;
    
    /*filter*/
    
    T=sigma;

    for(k=ini;k<=fin;k=k+1)
    {  
	for(j=0;j<=rows-B;j=j+1)
	{
		for(i=0;i<=cols-B;i=i+1)
		{	                                   
            max=0; 
            for(kk=0;kk<B;kk++)
	        for(jj=0;jj<B;jj++)	
		    for(ii=0;ii<B;ii++)
		    {	
              p=kk*B*B+jj*B+ii;
              p1=(k+kk)*rows*cols+(j+jj)*cols+(i+ii);
              b[p] = ima[p1];                     
              b2[p] = fima[p1];     
              if(b[p]>max) max=b[p];
            } 
            if(max==0) continue;
            
            dct3(b,cb);                                                                               
            dct3(b2,cb2);
            
            /*/////////////////////////////////////////////////*/
                       
            ind=0;           
            for(kk=0;kk<B;kk++)
	        for(jj=0;jj<B;jj++)	
		    for(ii=0;ii<B;ii++)
		    {
                p=kk*B*B+jj*B+ii;
               
                if(abs(cb2[p])<T) cb[p]=0; /* hardthreshold(cb,T);                */
                else ind++;                
            }
            z=1.0/(ind+1); 
            
            idct3(cb,b);
                                  
            for(kk=0;kk<B;kk++)
	        for(jj=0;jj<B;jj++)	
		    for(ii=0;ii<B;ii++)
		    {	
              p1=(k+kk)*rows*cols+(j+jj)*cols+(i+ii);  
              p=kk*B*B+jj*B+ii;
              
              if(rician>0)
              {                 
                iii=(int)(b[p]);                 
                val=mibias[iii];
                if(val>2*b[p]) val=b[p];                
                fima2[p1]+=z*val;
              }
              
              else fima2[p1]+=z*b[p];  /* agregation     */
              acu[p1]+=z;                                               
            }
		}
	} 
    }
 
    #ifdef _WIN32
       _endthreadex(0);    
    #else
       pthread_exit(0);    
    #endif

	return 0;
} 


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
/*Declarations*/
mxArray *xData,*pv;
double *ima,*fima,*fima2,*acu,*mibias;
double sigma,val,v1,max;
int i,ndim,Nthreads,ini,fin,N,r,k,ind,imax;
const int *dims;
mwSize *bdims;

myargument *ThreadArgs; 
 
#ifdef _WIN32
HANDLE *ThreadList; /* Handles to the worker threads*/
#else
pthread_t * ThreadList;
int *ThreadIndex;
#endif

/*Copy input pointer x*/
xData = prhs[0];
ima = mxGetPr(prhs[0]);
ndim = mxGetNumberOfDimensions(prhs[0]);
dims= mxGetDimensions(prhs[0]);

/*Copy input parameters*/
pv = prhs[1];
sigma = (double)(mxGetScalar(pv));

pv = prhs[2];
rician = (int)(mxGetScalar(pv));

/*Allocate memory and assign output pointer*/
plhs[0] = mxCreateNumericArray(ndim,dims,mxDOUBLE_CLASS, mxREAL);
fima = mxGetPr(plhs[0]);

xData = mxCreateNumericArray(ndim,dims,mxDOUBLE_CLASS, mxREAL);
acu = mxGetPr(xData);
xData = mxCreateNumericArray(ndim,dims,mxDOUBLE_CLASS, mxREAL);
fima2 = mxGetPr(xData);

max=0;
for(i=0;i<dims[0]*dims[1]*dims[2];i++)
{        
    fima[i]=0;    
    acu[i]=0;
    if(ima[i]>max) max=ima[i];
}					
imax=(int)ceil(2*max);

Nthreads = 8;/*floor((double)dims[2]/(double)B);*/

N=B;
tabla1=mxMalloc(N*N*sizeof(double));
tabla2=mxMalloc(N*N*sizeof(double));
for(r=0;r<N;r++) for(k=0;k<N;k++) tabla1[r*N+k]=cos((PI*(2*k+1)*r)/(2*N));               
for(r=0;r<N;r++) for(k=0;k<N;k++) tabla2[r*N+k]=cos((PI*(2*r+1)*k)/(2*N));      
W1=sqrt(1/(double)N);
W2=sqrt(2/(double)N);	


/*mibias=mxMalloc(max*sizeof(double));*/
bdims = (mwSize *) mxMalloc (2 * sizeof(mwSize));
bdims[0] = 1;
bdims[1] = imax*Nthreads;
pv = mxCreateNumericArray (2, bdims, mxDOUBLE_CLASS, mxREAL);
mibias = mxGetPr(pv);

if(rician>0) RicianBias(sigma,imax,mibias,Nthreads);


#ifdef _WIN32

/* Reserve room for handles of threads in ThreadList*/
ThreadList = (HANDLE*) mxMalloc(Nthreads*sizeof( HANDLE ));
ThreadArgs = (myargument*) mxMalloc(Nthreads*sizeof(myargument));

for (i=0; i<Nthreads; i++)
{       
	/* Make Thread Structure*/
    ini=ceil((i*(dims[2]-B))/Nthreads);
    if(i==Nthreads-1) fin=floor(((i+1)*(dims[2]-B))/Nthreads);
    else fin=floor(((i+1)*(dims[2]-B))/Nthreads)-1;   
       
    ThreadArgs[i].rows=dims[1];
	ThreadArgs[i].cols=dims[0];
   	ThreadArgs[i].slices=dims[2];
    ThreadArgs[i].in_image=ima;
	ThreadArgs[i].out_image=fima;
    ThreadArgs[i].acu_image=acu;
    ThreadArgs[i].ini=ini;
    ThreadArgs[i].fin=fin;    
    ThreadArgs[i].sigma=sigma;        
	
    ThreadList[i] = (HANDLE)_beginthreadex( NULL, 0, &ThreadFunc, &ThreadArgs[i] , 0, NULL );    
}
    
for (i=0; i<Nthreads; i++) { WaitForSingleObject(ThreadList[i], INFINITE); }
for (i=0; i<Nthreads; i++) { CloseHandle( ThreadList[i] ); }

for(i=0;i<dims[0]*dims[1]*dims[2];i++)
{
   if(acu[i]==0) acu[i]=1; 
   fima[i]/=acu[i];
}

/*//////////////////////////////////////////////*/
/* second pass*/
/*//////////////////////////////////////////////*/

for(i=0;i<dims[0]*dims[1]*dims[2];i++)
{        
    fima2[i]=0;    
    acu[i]=0;
}					
for (i=0; i<Nthreads; i++)
{       
	/* Make Thread Structure*/
    ini=ceil((i*(dims[2]-B))/Nthreads);
    if(i==Nthreads-1) fin=floor(((i+1)*(dims[2]-B))/Nthreads);
    else fin=floor(((i+1)*(dims[2]-B))/Nthreads)-1;    
       
    ThreadArgs[i].rows=dims[1];
	ThreadArgs[i].cols=dims[0];
   	ThreadArgs[i].slices=dims[2];
    ThreadArgs[i].in_image=ima;
	ThreadArgs[i].out_image=fima;
   	ThreadArgs[i].out_image2=fima2;
    ThreadArgs[i].acu_image=acu;
    ThreadArgs[i].ini=ini;
    ThreadArgs[i].fin=fin;   
    ThreadArgs[i].sigma=sigma;      
    ThreadArgs[i].mibias=&mibias[i*imax];      
	
    ThreadList[i] = (HANDLE)_beginthreadex( NULL, 0, &ThreadFunc2, &ThreadArgs[i] , 0, NULL );    
}
    
for (i=0; i<Nthreads; i++) { WaitForSingleObject(ThreadList[i], INFINITE); }
for (i=0; i<Nthreads; i++) { CloseHandle( ThreadList[i] ); }


#else

/* Reserve room for handles of threads in ThreadList*/
ThreadList = (pthread_t *) calloc(Nthreads,sizeof(pthread_t));
ThreadArgs = (myargument*) calloc( Nthreads,sizeof(myargument));
ThreadIndex = (int *)malloc(Nthreads* sizeof( int ));
    
for (i=0; i<Nthreads; i++)
{       
	/* Make Thread Structure*/
    ini=ceil((i*(dims[2]-B))/Nthreads);
    if(i==Nthreads-1) fin=floor(((i+1)*(dims[2]-B))/Nthreads);
    else fin=floor(((i+1)*(dims[2]-B))/Nthreads)-1;   
       
    ThreadArgs[i].rows=dims[1];
	ThreadArgs[i].cols=dims[0];
   	ThreadArgs[i].slices=dims[2];
    ThreadArgs[i].in_image=ima;
	ThreadArgs[i].out_image=fima;
    ThreadArgs[i].acu_image=acu;
    ThreadArgs[i].ini=ini;
    ThreadArgs[i].fin=fin;    
    ThreadArgs[i].sigma=sigma;        
}

for (i=0; i<Nthreads; i++)
{
    if(ThreadIndex[i] = pthread_create(&ThreadList[i], NULL, ThreadFunc, (void*) &ThreadArgs[i]))
        /*pthread_create( &ThreadList[i], NULL, ThreadFunc, (void*) &ThreadArgs[i]);*/
    {
       printf("Threads cannot be created\n");
       exit(1);
    }        
}
   
for (i=0; i<Nthreads; i++)
{
  pthread_join(ThreadList[i],NULL);
}

for(i=0;i<dims[0]*dims[1]*dims[2];i++)
{
   if(acu[i]==0) acu[i]=1; 
   fima[i]/=acu[i];
}

/*//////////////////////////////////////////////*/
/* second pass*/
/*//////////////////////////////////////////////*/

for(i=0;i<dims[0]*dims[1]*dims[2];i++)
{        
    fima2[i]=0;    
    acu[i]=0;
}	
				
for (i=0; i<Nthreads; i++)
{       
	/* Make Thread Structure*/
    ini=ceil((i*(dims[2]-B))/Nthreads);
    if(i==Nthreads-1) fin=floor(((i+1)*(dims[2]-B))/Nthreads);
    else fin=floor(((i+1)*(dims[2]-B))/Nthreads)-1;    
       
    ThreadArgs[i].rows=dims[1];
	ThreadArgs[i].cols=dims[0];
   	ThreadArgs[i].slices=dims[2];
    ThreadArgs[i].in_image=ima;
	ThreadArgs[i].out_image=fima;
   	ThreadArgs[i].out_image2=fima2;
    ThreadArgs[i].acu_image=acu;
    ThreadArgs[i].ini=ini;
    ThreadArgs[i].fin=fin;   
    ThreadArgs[i].sigma=sigma;      
    ThreadArgs[i].mibias=&mibias[i*imax];     
}

	
for (i=0; i<Nthreads; i++)
{
    if(ThreadIndex[i] = pthread_create(&ThreadList[i], NULL, ThreadFunc2,&ThreadArgs[i]))
    {
       printf("Threads cannot be created\n");
       exit(1);
    }        
}
   
for (i=0; i<Nthreads; i++)
{
  pthread_join(ThreadList[i],NULL);
}


#endif

 
for(i=0;i<dims[0]*dims[1]*dims[2];i++) 
{
   if(acu[i]==0) acu[i]=1;
   fima[i]=fima2[i]/acu[i];   
   if(rician) if(fima[i]<0) fima[i]=0; /* background*/
}


mxFree(ThreadArgs); 
mxFree(ThreadList);
mxFree(tabla1);
mxFree(tabla2);   

return;

}

