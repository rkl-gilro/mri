#include "math.h"
#include "mex.h"
#include <stdlib.h>
#include "matrix.h"


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{

/* Declarations*/
mxArray *xData;
double *ima, *fima;
mxArray *pv;
double media;
int i,j,k,ii,jj,kk,ni,nj,nk,v,f,ndim,indice;
const int  *dims;

/* Copy input pointer x*/
xData = prhs[0];

/*Get matrix x */
ima = mxGetPr(xData);

ndim = mxGetNumberOfDimensions(prhs[0]);
dims= mxGetDimensions(prhs[0]);

/* Copy input parameters */
pv = prhs[1];
/* Get the Integer*/
v = (int)(mxGetScalar(pv));

/*Allocate memory and assign output pointer*/

/*printf("ancho=%d\r",dims[0]);*/
/*printf("alto=%d\r",dims[1]);*/
/*printf("profundo=%d\r",dims[2]);*/

plhs[0] = mxCreateNumericArray(ndim,dims,mxDOUBLE_CLASS, mxREAL);

/*Get a pointer to the data space in our newly allocated memory*/
fima = mxGetPr(plhs[0]);


/*precomputo de momentos de orden 1 */

for(k=0;k<dims[2];k++)
{
for(i=0;i<dims[1];i++)
{
  for(j=0;j<dims[0];j++)
  {
        media=0;
        indice=0;
      	for(ii=-v;ii<=v;ii++)
        {
		  for(jj=-v;jj<=v;jj++)
          {
            for(kk=-v;kk<=v;kk++)
            {
              ni=i+ii;
		      nj=j+jj;		   		  
              nk=k+kk;		   		  
		      if(ni>=0 && nj>=0 && nk>0 && ni<dims[1] && nj<dims[0] && nk<dims[2])
             {
                media = media + ima[nk*(dims[0]*dims[1])+(ni*dims[0])+nj];
                indice=indice+1;
             }
            }
         }
        }
        media=media/indice;
        fima[k*(dims[0]*dims[1])+(i*dims[0])+j]=media;
  }
}
}


return;

}

