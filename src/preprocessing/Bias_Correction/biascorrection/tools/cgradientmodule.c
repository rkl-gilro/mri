#include "math.h"
#include "mex.h"
#include <stdlib.h>
#include "matrix.h"


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{

//Declarations
mxArray *xData;
double *ima, *fima;
mxArray *pv;
double media,dx,dy,dz;
int x1,x2,y1,y2,z1,z2,i,j,k,ndim;
const int  *dims;

//Copy input pointer x
xData = prhs[0];

//Get matrix x
ima = mxGetPr(xData);

ndim = mxGetNumberOfDimensions(prhs[0]);
dims= mxGetDimensions(prhs[0]);

//Allocate memory and assign output pointer

//printf("ancho=%d\r",dims[0]);
//printf("alto=%d\r",dims[1]);
//printf("profundo=%d\r",dims[2]);

plhs[0] = mxCreateNumericArray(ndim,dims,mxDOUBLE_CLASS, mxREAL);

//Get a pointer to the data space in our newly allocated memory
fima = mxGetPr(plhs[0]);

//precomputo de momentos de orden 1 y 2

for(k=0;k<dims[2];k++)
{
for(i=0;i<dims[1];i++)
{
  for(j=0;j<dims[0];j++)
  {        
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


return;

}

