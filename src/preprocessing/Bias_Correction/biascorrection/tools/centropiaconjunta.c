#include <stdlib.h>
#include <stdio.h>
#include "math.h"
#include <float.h>
#include "mex.h"

void gradiente2d(double *ima,double* fima, int *dims )
{

int i,j,x1,x2,y1,y2;

double dx,dy;
        
    for(i=0;i<dims[1];i++){
        
        for(j=0;j<dims[0];j++){
            
            if(ima[(i*dims[0])+j]>0)
            {
            x1=((j+1)<dims[0])?(j+1):j;
            x2=(j-1)<0?0:(j-1);
            dx=(ima[(i*dims[0])+x1]-ima[(i*dims[0])+x2])/(x2-x1);
   
            y1=((i+1)<dims[1])?(i+1):i;
            y2=(i-1)<0?0:(i-1);
            dy=(ima[(y1*dims[0])+j]-ima[(y2*dims[0])+j])/(y2-y1);
   
            fima[(i*dims[0])+j]=sqrt(dx*dx+dy*dy);
            }
            else fima[(i*dims[0])+j]=0;
            }
        }
    return fima;
}

void gradiente3d(double *ima,double* fima, int *dims )
{

int i,j,k,x1,x2,y1,y2,z1,z2;

double dx,dy,dz;
    
for(k=0;k<dims[2];k++){
    
    for(i=0;i<dims[1];i++){
        
        for(j=0;j<dims[0];j++){
            
            if(ima[k*(dims[0]*dims[1])+(i*dims[0])+j]>0)
            {
            x1=((j+1)<dims[0])?(j+1):j;
            x2=(j-1)<0?0:(j-1);
            dx=(ima[k*(dims[0]*dims[1])+(i*dims[0])+x1]-ima[k*(dims[0]*dims[1])+(i*dims[0])+x2])/(x2-x1);
   
            y1=((i+1)<dims[1])?(i+1):i;
            y2=(i-1)<0?0:(i-1);
            dy=(ima[k*(dims[0]*dims[1])+(y1*dims[0])+j]-ima[k*(dims[0]*dims[1])+(y2*dims[0])+j])/(y2-y1);
   
            z1=((k+1)<dims[2])?(k+1):k;
            z2=(k-1)<0?0:(k-1);
            dz=(ima[z1*(dims[0]*dims[1])+(i*dims[0])+j]-ima[z2*(dims[0]*dims[1])+(i*dims[0])+j])/(z2-z1);
   
            fima[k*(dims[0]*dims[1])+(i*dims[0])+j]=sqrt(dx*dx+dy*dy+dz*dz);
            }
            else fima[k*(dims[0]*dims[1])+(i*dims[0])+j]=0;
            }
        }
    }
}

double maxi(double* array,int leng)
{
    int i;
    double maxi;
    
    maxi = array[0];       /* first element is the max */
 
    for(i = 1; i<leng; i++)
    {
        if(array[i] > maxi)
            maxi = array[i];
    }
    return maxi;               

}

void conv2(double *h, double *w, int m1, int m2)
{

int i,j,ii,jj,ni,nj,v;
double indice, media;
double *temp;
temp = (double*) calloc(m1*m2,sizeof(double)); 
        
v = 2;

    for(i=0;i<m1;i++)
    {
        for(j=0;j<m2;j++)
        {
            media=0;
            indice=0;
            for(ii=-v;ii<=v;ii++)
            {
                for(jj=-v;jj<=v;jj++)
                {
                    ni=i+ii;
                    nj=j+jj;		   		  		   		  
                    if(ni>=0 && nj>=0 && ni<m1 && nj<m2)
                    {
                        media = media + w[(ii+2) *5 + (jj+2)]*h[(ni*m2)+nj];
                        indice = indice + w[(ii+2) *5 + (jj+2)];
                    }
                }
            }
        
        media=media/indice;
        temp[(i*m2)+j] = media;
        }
    }
    for (i=0;i<m1*m2; i++) h[i] = temp[i];

free(temp);
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{

mxArray *xData;
double *ima, *h,*temp, *B, binsize, *e, suma;
const int  *dims;
int i,ndim;
int f, m1,m2,g, g1, g2;
double f1,f2, sum;

double w[] = {0.5, 0.75, 1, 0.75, 0.5, 0.75  , 1.5, 2 , 1.5, 0.75, 1  , 2  , 3 , 2  , 1,0.75  , 1.5, 2 , 1.5, 0.75,0.5, 0.75  , 1 , 0.75  , 0.5};
binsize = 0.01;
   
/*read imput parameters*/

/*Copy input pointer x */
xData = prhs[0];
/* Get matrix x */
ima = mxGetPr(xData);

/*number of dimensions of the input matrix*/       
ndim = mxGetNumberOfDimensions(prhs[0]);

/*size of the current matrix*/ 
dims= mxGetDimensions(prhs[0]);  

/*printf("ndim = %d\n", ndim);*/

/*Allocate memory and assign output pointer*/
/*plhs[0] = mxCreateNumericArray(ndim,dims,mxDOUBLE_CLASS, mxREAL);*/
plhs[0] = mxCreateDoubleMatrix(1,1, mxREAL);

/*Get a pointer to the data space in our newly allocated memory*/
e = mxGetPr(plhs[0]);
nlhs = 1;

if(ndim == 2)
{
  temp = (double*) malloc(dims[1]*dims[0]*sizeof(double));
  B = (double*) malloc(dims[1]*dims[0]*sizeof(double)); /*reservar memoria*/
  
  /*make a copy */
  for(i=0;i<dims[1]*dims[0]; i++)
    temp[i]=ima[i];

  /*data=abs(log(data+1))/binsize; */ 
  for(i=0;i<dims[1]*dims[0]; i++) 
    temp[i]=fabs(log(temp[i]+1))/binsize; /* fabs - double abs (double number)  */
    
  gradiente2d(temp,B,dims); /* vao ser atribuidos valores a B */
  
  for(i=0;i<dims[1]*dims[0]; i++) B[i] = B[i] + 1;
  
  m1 = ceil(maxi(temp, dims[1]*dims[0]) + 1); /* warning � normal visto que se perde a parte fraccionaria ao fazer o ceil*/
  m2 = ceil(maxi(B, dims[1]*dims[0]) + 1);
  
  h = (double*) calloc(m1*m2,sizeof(double));  
  
  for(i=0;i<dims[1]*dims[0]; i++) 
    {
        if(floor(temp[i])>0 && floor(B[i])>0)
        {
            f = floor(temp[i]);   
            g = floor(B[i]);                    
            f2 = temp[i]-f;
            f1 = 1-f2;
            g2 = B[i]-g;
            g1 = 1-g2;
      
            /* h(f,g)    = h(f,g)    + f1*g1;   */
            h[f*(m2-1)+g] = h[f*(m2-1)+g] + f1*g1;   
      
            /* h(f+1,g)  = h(f+1,g)  + f2*g1;*/
            h[(f+1)*(m2-1)+g] = h[(f+1)*(m2-1)+g] + f2*g1;
      
            /* h(f,g+1)  = h(f,g+1)  + f1*g2;*/
            h[f*(m2-1)+(g+1)] = h[f*(m2-1)+(g+1)] + f1*g2;
      
            /* h(f+1,g+1)= h(f+1,g+1)+ f2*g2;  */
            h[(f+1)*(m2-1)+g+1] = h[(f+1)*(m2-1)+g+1] + f2*g2;
        }
    } 
  
}

if(ndim == 3)
{
  temp = (double*) malloc(dims[2]*dims[1]*dims[0]*sizeof(double));
  B = (double*) malloc(dims[2]*dims[1]*dims[0]*sizeof(double)); /*reservar memoria*/

  /*make a copy */
  for(i=0;i<dims[2]*dims[1]*dims[0]; i++) temp[i]=ima[i];

  /*data=abs(log(data+1))/binsize; */ 
  for(i=0;i<dims[2]*dims[1]*dims[0]; i++) 
     temp[i]=fabs(log(temp[i]+1))/binsize; /* fabs - double fabs (double number) */
  
  gradiente3d(temp,B,dims); /* vao ser atribuidos valores a B */
  
  for(i=0;i<dims[2]*dims[1]*dims[0]; i++) B[i] = B[i] +1;
  
  m1 = ceil(maxi(temp, dims[2]*dims[1]*dims[0]) + 1); /* warning � normal visto que se perde a parte fraccionaria ao fazer o ceil*/
  m2 = ceil(maxi(B, dims[2]*dims[1]*dims[0]) + 1); 

  h = (double*) calloc(m1*m2,sizeof(double));  
  
   for(i=0;i<dims[2]*dims[1]*dims[0]; i++) 
    {
        if(floor(temp[i])>0 && floor(B[i])>0)   
        {
            f = floor(temp[i]);   
            g = floor(B[i]);                    
            f2 = temp[i]-f; 
            f1 = 1-f2;
            g2 = B[i]-g;
            g1 = 1-g2;
      
            /* h(f,g)    = h(f,g)    + f1*g1;   */
            h[f*(m2-1)+g] = h[f*(m2-1)+g] + f1*g1;
            
            /* h(f+1,g)  = h(f+1,g)  + f2*g1; */
            h[(f+1)*(m2-1)+g] = h[(f+1)*(m2-1)+g] + f2*g1;
            
            /* h(f,g+1)  = h(f,g+1)  + f1*g2; */
            h[f*(m2-1)+(g+1)] = h[f*(m2-1)+(g+1)] + f1*g2;
            
            /* h(f+1,g+1)= h(f+1,g+1)+ f2*g2; */
            h[(f+1)*(m2-1)+g+1] = h[(f+1)*(m2-1)+g+1] + f2*g2;
            
        }
    } 
}

conv2(h, w, m1, m2);

sum = 0;
for (i=0; i<m1*m2; i++) sum = sum + h[i];

for (i = 0 ; i < m1*m2  ; i++) h[i]=h[i]/sum;

e[0] = 0;
for (i = 0 ; i < m1*m2  ; i++)
    if(h[i]>0)
    {
        e[0] = e[0] + h[i]*log(h[i]);
    }     
    e[0] = - e[0];  /* � um valor negativo*/
    
    
  free(temp);
  free(B);
  free(h);

}


