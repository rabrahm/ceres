#include <Python.h>
#include <numpy/arrayobject.h>
#include <gsl/gsl_linalg.h>
#include <sys/time.h>
#define ARRAYD(p) ((double *) (((PyArrayObject *)p)->data)) 

/* 
 *                                [INITIALIZATION]
 * ------------------ PROTOTYPES FOR FUNCTIONS AND EXTERNAL VARIABLES -----------------------
 *
 */

int CheckAperture(double* center, int len_rows, int len_cen, double S, int CentralAperture);
double Polyval(double x,double* coeffs, int len_coeffs);
double** MakeArray(int rows, int columns);                           /* Function that makes/allocates an array of pointers    */
double* MakeVector(int nelements);                                   /* Function that makes/allocates a pointer in memory     */
double* RangeDetector(double *v,double C,double S,int len_v,int K,int len_rows,int len_cols,int range); /* Detects the range of the                                                                                                      calculations                */
double* SimpleRangeDetector(double *v,double Length,int len_v,int len_rows,int len_cols,int range);
void FreeArray(double** theArray,int rows);                          /* Function that frees an array of pointers from memory  */
void PixelResampling(double** A,double* pmin,double* pmax,int len_cols);                                    /* Algorithm for pixel resampling */
void BPixelResampling(double** A,double** B,double* pmin,double* pmax,int len_cols);                                    /* Algorithm for pixel resampling */
void SimplePixelResampling(double** A,double* pmin,double* pmax,int len_cols);
void LinearSolver(double** A,double* b,double* x,int length);        /* Linear-algebra system solver (returns the solution
                                                                        in the x vector)                                      */
double** getQ(double *v,double C,double S,int len_v,int K,int len_cols,int range,double* pmin,double* pmax,int mode);                                                 /* Function that gets the Q values (returns a matrix)    */
void getJ(double** J,int N,int len_cols);                            /* Function that calculates the powers of j, the columns */
void getC(double** A,double** Q,double** VarE,double** C_qp,double** J,double* pmin,double* pmax,int K,int N,int len_cols); /* Function that gets the C values (For the L.S.)    */
void getX(double** A,double** Q,double** E,double** VarE,double* X,double** J,double* pmin,double* pmax,int K,int N,int len_cols); /* Function that gets the X values            */
void getP(double* B,double** Q,double** P,double** J,double* pmin,double* pmax,int K,int N,int len_rows,int len_cols,int mode);               /*     Computes the light fraction matrix (matrix)       */
void getRowSum(double** A,double* RS,double *v,double* pmin,double* pmax,int len_cols,int len_v);                               /* Computes the sum of the vertical elements in a row    */
void getE(double** A,double* RS,double** E,double* pmin,double* pmax,double *v,int len_rows,int len_cols);                         /* Computes the light fraction estimates E_{ij}          */
void getVarRowSum(double** A,double** V,double* VarRS,double* pmin,double* pmax,double RON,double GAIN,int len_cols);              /* Same as RowSum but for variances                      */
void getVarE(double** A,double** E,double* RS,double* VarRS,double** VarE,double** V,double* pmin,double* pmax,double RON,double GAIN,int len_rows,int len_cols); /* Same as E for Var(E)                */
double CalculateQ(double *v,int k,int i,int j,double C,double S,int len_v,int range);  /* Calculates the value of Q for given parameters        */
double CalculateC(double** A,double** Q,double** VarE,int m, int l, int n, int k,double** J,double* pmin,double* pmax,int len_cols); /* Same as above for C          */
double CalculateX(double** A,double** Q,double** E,double** VarE,int n, int k,double** J,double* pmin,double* pmax,int len_cols); /* Same as above for X             */
double Trace(int x,double *v,int len_v,int range);                                              /* Computes the trace for a given value of x             */
void getMinMax(double* pmin,double* pmax,double *v,double C, double S,int K,int len_v,int len_cols,int range);
void getSimpleMinMax(double* pmin,double* pmax,double *v,int len_v,int len_cols,int range,double Length);

double** MarshIteration(double** M,double** VarImage,double* pmin,double* pmax,double *v,double C,double S,double RON,double GAIN,int len_v,int len_rows,int len_cols,int N,int K,int mode,int range,int debug,double NSigma,int isvarimage);                                 /* Function that iterates the whole Marsh's algorithm    */
double** Transpose(double** O,int rowsIN, int colsIN);               /* Function that obtains the tranpose of a matrix        */
void Renormalize(double** P,double* pmin,double* pmax,int len_rows,int len_cols);     
void VarRevision(double** A,double** VarImage,double** V,double* RS,double** P,double* pmin,double* pmax,int len_cols,double RON,double GAIN,int isvarimage);       /* Function that revise the variance estimates           */
int OutlierRejection(double** A,double** V,double** P,double** E,double** VarE,double* RS,double* pmin,double* pmax,double NSigma,int len_cols,int debug); /* Out. Rej. algorithm           */
double** getSpectrum(double** A,double** VarImage,double** P,double *v,int len_v,double* pmin,double* pmax,double RON,double GAIN,int len_rows,int len_cols,double CosmicSigma,int debug,int isvarimage);                         /* Function that returns the 1-D spectrum matrix given P.
																												     The first row of the matrix are the pixel values in 
																												     matrix coordinates (lambda), the second row is F_{lambda}
																												     and the third is var(F_{lambda}) */
double** BgetSpectrum(double** A,double** B,double** P,double *v,int len_v,double* pmin,double* pmax,double RON,double GAIN,int len_rows,int len_cols,double CosmicSigma,int debug);                         /* Function that returns the 1-D spectrum matrix given P.
																										The first row of the matrix are the pixel values in 
																										matrix coordinates (lambda), the second row is F_{lambda}
																										and the third is var(F_{lambda}) */
void getImageVariances(double** A,double** V,double** VarImage,double* pmin,double* pmax,double RON,double GAIN,int len_cols,int isvarimage);   /* Function that returns a matrix with the values of the 
																		   initial variances in each pixel                       */
void getRowSumW(double** A,double** V,double** P,double* pmin,double* pmax,double* RSW,int len_cols);        /* Function that returns the denominator in eq. (4) of
														Marsh (1989)                                          */
void getW(double** A,double** P,double** V,double* RSW,double** W,double* pmin,double* pmax,int len_cols);  /* Function that returns the weigths for the weighted 
													       averaged spectra                                     */
void getF(double** W,double** A,double* F,double* pmin,double* pmax,int len_cols);                          /* Vector containing the fluxes calculated from the 
													       weighted average                                      */
void getBF(double** W,double** B,double* BF,double* pmin,double* pmax,int len_cols);                          /* Vector containing the fluxes calculated from the 
														 weighted average                                      */
void getVarF(double** A,double** W,double** V,double* VarF,double* pmin,double* pmax,int len_cols);         /* Vector containing the variances of fluxes calculated 
													       from the weighted average                             */
int CosmicRayRejection(double** A,double** V,double** P,double* F,double* VarF,double CosmicSigma,double* pmin,double* pmax,int len_cols,int debug);  /* Cosmic Rays rejection algorithm         */

int BCosmicRayRejection(double** A,double** B,double** V,double** P,double* F,double* VarF,double CosmicSigma,double* pmin,double* pmax,int len_cols,int debug);  /* Cosmic Rays rejection algorithm         */

/*                  [INITIALIZATION OF A METHOD]
 *------------------------THE OBTAINP METHOD-----------------------------
 * PyObject initialization: We define a PyObject which defines a Method 
 * for the Marsh module: The ObtainP method which returns the P[i][j] spatial
 * light fractions back to Python. BE CAREFUL WITH THE INDEXES, remember
 * that i=rows, j=columns. According to Marsh's (1989) variables, j=X and
 * i=Y.
 *----------------------------------------------------------------------
 */

static PyObject *Marsh_ObtainP(PyObject *self, PyObject *args){
  struct timeval tim;
  gettimeofday(&tim, NULL);
  double t1=tim.tv_sec+(tim.tv_usec/1000000.0);
  double *m,dvalue=0;
  double *v,*pmin,*pmax; // pmin,pmax: FREED (on the MarshIteration function)
  double C,S,Length,RON,GAIN,NSigma;
  int len_v,K,N,len_cols,len_rows,mode,range=-1,debug,isvarimage=0;
  int i,j,real_cols,real_rows,col_min,col_max;
  PyObject *marray, *varray;
	
  /* 
   *--------------------------------THE DATA---------------------------------------
   * After initialization of the PyObject pointers, we wish to recover the following inputs:
   *
   * marray:   Vector of the flattened-matrix of the Echelle Spectra data.
   *
   * varray:   Vector containing the coefficients of the trace of the spectrum.
   *
   * len_rows: Length of the rows of the flattened-matrix.
   *
   * len_cols: Length of the columns of the flattened-matrix.
   *
   * len_v:    Length of the coefficient-vector.
   *
   * Length:   Given that the trace is centered on the spectrum, Length gives the length,
   *           in pixels from the center, where we'll consider the spectrum.
   *
   * RON:      The Read Out Noise of our measurements in electrons.
   *
   * GAIN:     The Gain of the detector, in electrons/ADU.
   * 
   * mode:     Set this to 0 to apply Marsh's Algorithm (curved trace) to the data.
   *           Set this to 1 to apply Horne's Algorithm (trace paralell to the rows) to the data.
   *
   * ------------------------------------------------------------------------------
   */
  PyArg_ParseTuple(args,"OOiiidddddiiii",&marray,&varray,&len_rows,&len_cols,&len_v,&Length,&RON,&GAIN,&NSigma,&S,&N,&mode,&col_min,&col_max);
  if(col_min!=0){
    col_min=col_min+1;
  }
  if(col_max!=0){
    col_max=col_max-1;
  }
  debug=0;
  m = ARRAYD(marray);                               /* We convert our PyObject struct pointers to C-vector array  */
  v = ARRAYD(varray);
  real_cols=len_cols;
  real_rows=len_rows;
  if(mode!=0)
    S=1.0;                                          /* If Horne's algorithm is used, the spacing between
						       polynomials is 1                                              */
  Length = CheckAperture(v, len_rows, len_v, S, Length);
  K=(int)(2*(int)((Length/S)+0.5)+1);               /* The number of polynomials. The idea is to form n=Length/S 
						       polynomials, rounding the decimal number, from the center.
						       Then we multiply this number by 2 (to get the total number
						       of polynomials up and down the center) and add 1 to get the
						       middle polynomial (the one corresponding to the center).      */
  if(debug!=0){
    printf("Number of polynomial to fit: %d \n",K);
  }
  C=-(S*(1.0+(((double)K-1.0)/2.0)));               /* We get the C constant for the Q coefficients. This number
						       is the polynomial number of the center.                       */

  double* Range = RangeDetector(v,C,S,len_v,K,len_rows,len_cols,range);  /* This vector containts the detected starting 
									    and final column pixels that "swims in" and "out" of the image, respectively
									    (FREED here, later). */
  if(debug!=0){
    printf("Range[0]: %f\n",Range[0]);
    printf("Range[1]: %f\n",Range[1]);
    printf("len_cols=%d, len_rows=%d \n",len_cols,len_rows);
  }
  range=Range[0];

  if(Range[0]!=-1){
    if(col_min==0 && col_max==0){
      len_cols=Range[1]-Range[0]+1;                   /* If the column presents problems, we start by creating a    
							 cutted image                                               */
    }
    else if((col_min!=0 && col_max!=0) || (col_min==0 && col_max!=0)){
      if(Range[0]>=col_min && Range[1]<=col_max){
	len_cols=Range[1]-Range[0]+1;
      }
      else if(Range[0]<col_min && Range[1]<=col_max){
	range=col_min;
	Range[0]=col_min;
	len_cols=Range[1]-Range[0]+1;
      }
      else if(Range[0]<col_min && Range[1]>col_max){
	range=col_min;
	Range[0]=col_min;
	Range[1]=col_max;
	len_cols=Range[1]-Range[0]+1;
      }
      else if(Range[0]>=col_min && Range[1]>col_max){
	Range[1]=col_max;
	len_cols=Range[1]-Range[0]+1;
      }
    }
    else if(col_min!=0 && col_max==0){
      if(Range[0]>=col_min){
	len_cols=Range[1]-Range[0]+1;
      }
      else if(Range[0]<col_min){
	range=col_min;
	Range[0]=col_min;
	len_cols=Range[1]-Range[0]+1;
      }
    }
  }
  else if((col_min!=0 && col_max!=0) || (col_min==0 && col_max!=0)){
    range=col_min;
    Range[0]=col_min;
    Range[1]=col_max;
    len_cols=Range[1]-Range[0]+1;
  }
  else if(col_min!=0 && col_max==0){
    range=col_min;
    Range[0]=col_min;
    Range[1]=len_cols-1;
    len_cols=Range[1]-Range[0]+1;
  }

  double** M = MakeArray(len_rows,len_cols);        /* Create a matrix array for our image (FREED on MarshIt.)    */
  double** DummyM = MakeArray(len_rows,len_cols);   /* This is just a dummy array to be consistent with the variance image */
  if(Range[0]==-1){
    for(i=0;i<len_rows;i++){                        /* If the column doesn't present problems, proceed as usual   */
      for(j=0;j<len_cols;j++){ 
	M[i][j]=m[i*real_cols+j]-dvalue;
	DummyM[i][j]=m[i*real_cols+j]-dvalue;
      }
    }
  }
  else{
    for(i=0;i<len_rows;i++){                        /* If it does, we fill our cutted image                       */
      for(j=range;j<Range[1]+1;j++){
	M[i][j-range]=m[i*real_cols+j]-dvalue;
	DummyM[i][j-range]=m[i*real_cols+j]-dvalue;
      }
    }
  }

  /* Start of the Marsh's Algorithm iteration: */

  pmin=MakeVector(len_cols);
  pmax=MakeVector(len_cols);
  getMinMax(pmin,pmax,v,C,S,K,len_v,len_cols,range);
  double** PF = MarshIteration(Transpose(M,len_rows,len_cols),Transpose(DummyM,len_rows,len_cols),pmin,pmax,v,C,S,RON,GAIN,len_v,len_rows,len_cols,N,K,mode,range,debug,NSigma,isvarimage); // FREED (here, later).
  double** PFinal = MakeArray(real_rows,real_cols);  // FREED (here, later).
  if(Range[0]==-1){
    for(i=0;i<len_rows;i++){                       
      for(j=0;j<len_cols;j++){ 
	PFinal[i][j]=PF[i][j];
      }
    }
  }
  else{
    for(i=0;i<len_rows;i++){                        /* If it does, we fill our cutted image                       */
      for(j=range;j<Range[1]+1;j++){ 
	PFinal[i][j]=PF[i][j-range];
      }
    }
    len_cols=real_cols;
  }
  FreeArray(PF,len_rows);
  gettimeofday(&tim, NULL);
  double t2=tim.tv_sec+(tim.tv_usec/1000000.0);
  if(debug!=0){
    printf("%2.2lf seconds elapsed\n", t2-t1);
  }
  // /* End of the Marsh's Algorithm iteration! */
  // 
  // /* Start of the matrix-to-vector conversion part of the code */    
  double* theArray;  // FREED (here, later).
  theArray = (double*) malloc((len_cols*len_rows)*sizeof(double));
  for(i=0;i<len_rows;i++){
    for(j=0;j<len_cols;j++){
      theArray[i*len_cols+j]=PFinal[i][j];
      if(PFinal[i][j]>1){
	theArray[i*len_cols+j]=0;
      }
    } 
  }
  FreeArray(PFinal,len_rows);
  free(Range);
  /* End of the matrix-to-vector conversion part of the code. Everything's ok up to here...now I have my theArray[i][j] array (matrix) ready to be called...*/


  /* Finally, we create a Python "Object" List that contains the P coefficients and return it back to Python */

  PyObject *lst = PyList_New(len_rows*len_cols);
  PyObject *num;
  if (!lst)
    return NULL;
  for (i = 0; i < len_rows*len_cols; i++) {
    num=PyFloat_FromDouble(theArray[i]);
    if (!num) {
      Py_DECREF(lst);
      return NULL;
    }
    PyList_SET_ITEM(lst, i, num);
  }
  free(theArray);
  PyObject *MyResult = Py_BuildValue("O",lst);
  Py_DECREF(lst);
  return MyResult;
}

/*                        [INITIALIZATION OF A METHOD]
 * ------------------------THE OBTAINSpectrum METHOD-----------------------------
 */

static PyObject *Marsh_ObtainSpectrum(PyObject *self, PyObject *args){
  struct timeval tim;
  gettimeofday(&tim, NULL);
  int i,j,real_cols,real_rows,col_min,col_max;
  double t1=tim.tv_sec+(tim.tv_usec/1000000.0),dvalue=0;
  double *mflat,*pflat;
  double *v,*pmin,*pmax; // (pmin,pmax) freed here, later.
  double Length,RON,GAIN,CosmicSigma,S,C;
  int len_v,len_cols,len_rows,K,range=-1,debug,isvarimage=0;
  PyObject *parray, *marray,*varray;
  /* 
   *--------------------------------THE DATA---------------------------------------
   * After initialization of the PyObject pointers, we wish to recover the following inputs:
   *
   * marray:   Vector of the flattened-matrix of the Echelle Spectra data.
   *
   * varray:   Vector containing the coefficients of the trace of the spectrum.
   *
   * parray:   Vector of the flattened-matrix of the Echelle Spectra P_{ij} coefficents.
   *
   * len_rows: Length of the rows of the flattened-matrix.
   *
   * len_cols: Length of the columns of the flattened-matrix.
   *
   * len_v:    Length of the coefficient-vector.
   * 
   * Length:   Aperture to be used in the obtention of the spectrum.
   *
   * RON:      The Read Out Noise of our measurements.
   * 
   * GAIN:     The Gain of the detector, in electrons/ADU.
   * 
   * S:        Spacing of the polynomials obtained in the ObtainP function.
   * 
   * CosmicSigma: Number of times sigma is multiplied by to reject Cosmic Rays.
   *
   * ------------------------------------------------------------------------------
   */
  PyArg_ParseTuple(args,"OOOiiidddddii",&marray,&varray,&parray,&len_rows,&len_cols,&len_v,&Length,&RON,&GAIN,&S,&CosmicSigma,&col_min,&col_max);
  if(col_min!=0){
    col_min=col_min+1;
  }
  if(col_max!=0){
    col_max=col_max-1;
  }
  debug=0;
  pflat = ARRAYD(parray);                               /* We convert our PyObject struct pointers to C-vector array  */
  mflat = ARRAYD(marray);                               /* We convert our PyObject struct pointers to C-vector array  */
  v = ARRAYD(varray);
  real_cols=len_cols;
  real_rows=len_rows;
  Length = CheckAperture(v, len_rows, len_v, S, Length);
  /* Start of the algorithm for the spectrum obtention (faster with the transpose!): */
	
  K=(int)(2*(int)((Length/S)+0.5)+1);               /* The number of polynomials. The idea is to form n=Length/S 
						       polynomials, rounding the decimal number, from the center.
						       Then we multiply this number by 2 (to get the total number
						       of polynomials up and down the center) and add 1 to get the
						       middle polynomial (the one corresponding to the center).      */
  C=-(S*(1.0+(((double)K-1.0)/2.0)));               /* We get the C constant for the Q coefficients. This number
						       is the polynomial number of the center.                       */
  range=0;
  double* Range = RangeDetector(v,C,S,len_v,K,len_rows,len_cols,range);                  /* This vector containts the
											    detected starting and final column  pixels that "swims in" and "out" of the image, respectively (FREED here) */
  range=Range[0];
  if(debug!=0){
    printf("Range[0]: %f\n",Range[0]);
    printf("Range[1]: %f\n",Range[1]);
    dvalue=1018.0;
  } 

  if(Range[0]!=-1){
    if(col_min==0 && col_max==0){
      len_cols=Range[1]-Range[0]+1;                   /* If the column presents problems, we start by creating a    
							 cutted image                                               */
    }
    else if((col_min!=0 && col_max!=0) || (col_min==0 && col_max!=0)){
      if(Range[0]>=col_min && Range[1]<=col_max){
	len_cols=Range[1]-Range[0]+1;
      }
      else if(Range[0]<col_min && Range[1]<=col_max){
	range=col_min;
	Range[0]=col_min;
	len_cols=Range[1]-Range[0]+1;
      }
      else if(Range[0]<col_min && Range[1]>col_max){
	range=col_min;
	Range[0]=col_min;
	Range[1]=col_max;
	len_cols=Range[1]-Range[0]+1;
      }
      else if(Range[0]>=col_min && Range[1]>col_max){
	Range[1]=col_max;
	len_cols=Range[1]-Range[0]+1;
      }
    }
    else if(col_min!=0 && col_max==0){
      if(Range[0]>=col_min){
	len_cols=Range[1]-Range[0]+1;
      }
      else if(Range[0]<col_min){
	range=col_min;
	Range[0]=col_min;
	len_cols=Range[1]-Range[0]+1;
      }
    }
  }
  else if((col_min!=0 && col_max!=0) || (col_min==0 && col_max!=0)){
    range=col_min;
    Range[0]=col_min;
    Range[1]=col_max;
    len_cols=Range[1]-Range[0]+1;
  }
  else if(col_min!=0 && col_max==0){
    range=col_min;
    Range[0]=col_min;
    Range[1]=len_cols-1;
    len_cols=Range[1]-Range[0]+1;
  }

  double** M = MakeArray(len_rows,len_cols);        /* Create a matrix array for our image (FREED on getSpectrum) */
  double** P = MakeArray(len_rows,len_cols);        /* Create a matrix array for our P's (FREED on getSpectrum)   */
  double** DummyM = MakeArray(len_rows,len_cols);
  if(Range[0]==-1){
    for(i=0;i<len_rows;i++){                        /* If the column doesn't present problems, proceed as usual   */
      for(j=0;j<len_cols;j++){ 
	M[i][j]=mflat[i*real_cols+j]-dvalue;
	DummyM[i][j]=mflat[i*real_cols+j]-dvalue;
	P[i][j]=pflat[i*real_cols+j];
      }
    }
  }
  else{
    for(i=0;i<len_rows;i++){                        /* If it does, we fill our cutted image                       */
      for(j=range;j<Range[1]+1;j++){ 
	M[i][j-range]=mflat[i*real_cols+j]-dvalue;
	DummyM[i][j-range]=mflat[i*real_cols+j]-dvalue;
	P[i][j-range]=pflat[i*real_cols+j];
      }
    }
  }

  pmin=MakeVector(len_cols);
  pmax=MakeVector(len_cols);
  getMinMax(pmin,pmax,v,C,S,K,len_v,len_cols,range);
  double** Spectrum=getSpectrum(Transpose(M,len_rows,len_cols),Transpose(DummyM,len_rows,len_cols),Transpose(P,len_rows,len_cols),v,len_v,pmin,pmax,RON,GAIN,len_rows,len_cols,CosmicSigma,debug,isvarimage); // FREED (here, later)
  gettimeofday(&tim, NULL);
  double t2=tim.tv_sec+(tim.tv_usec/1000000.0);
  if(debug!=0){
    printf("%2.2lf seconds elapsed\n", t2-t1);
  }
  /* End of the algorithm */

  double* theArray; // FREED (here, later)
  theArray = (double*) malloc((real_cols*3)*sizeof(double));
  for(i=0;i<3;i++){
    for(j=0;j<real_cols;j++){
      if(range==-1){
	theArray[i*real_cols+j]=Spectrum[i][j];
      }
      else if(range!=-1 && i==0){
	theArray[i*real_cols+j]=j;
      }
      else if(range!=-1 && i!=0){
	if(j<Range[0] || j>Range[1])
	  theArray[i*real_cols+j]=0;
	else
	  theArray[i*real_cols+j]=Spectrum[i][j-range];
      }
    } 
  }

  FreeArray(Spectrum,3);
  free(Range);

  /* Finally, we create a Python "Object" List that contains the P coefficients and return it back to Python */

  PyObject *lst = PyList_New(real_cols*3);
  if (!lst)
    return NULL;
  for (i = 0; i < 3*real_cols; i++) {
    PyObject *num = PyFloat_FromDouble(theArray[i]);
    if (!num) {
      Py_DECREF(lst);
      return NULL;
    }
    PyList_SET_ITEM(lst, i, num);
  }
  free(theArray);
  free(pmin); 
  free(pmax);
  PyObject *MyResult = Py_BuildValue("Oi",lst,real_cols);
  Py_DECREF(lst);
  return MyResult;
}


/*                  [INITIALIZATION OF A METHOD]
 *------------------------THE SOBTAINP METHOD-----------------------------
 * PyObject initialization: We define a PyObject which defines a Method 
 * for the Marsh module: The ObtainP method which returns the P[i][j] spatial
 * light fractions back to Python. BE CAREFUL WITH THE INDEXES, remember
 * that i=rows, j=columns. According to Marsh's (1989) variables, j=X and
 * i=Y.
 *----------------------------------------------------------------------
 */

static PyObject *Marsh_SObtainP(PyObject *self, PyObject *args){
  struct timeval tim;
  gettimeofday(&tim, NULL);
  double t1=tim.tv_sec+(tim.tv_usec/1000000.0);
  double *m,*varianceimage,*bimage,dvalue=0;
  double *v,*pmin,*pmax; // pmin,pmax: FREED (on the MarshIteration function)
  double C,S,Length,RON,GAIN,NSigma;
  int len_v,K,N,len_cols,len_rows,mode,range=-1,debug,isvarimage=1;
  int i,j,real_cols,real_rows,col_min,col_max;
  PyObject *marray, *varray,*vararray,*barray;
	
  /* 
   *--------------------------------THE DATA---------------------------------------
   * After initialization of the PyObject pointers, we wish to recover the following inputs:
   *
   * marray:   Vector of the flattened-matrix of the Echelle Spectra data.
   *
   * varray:   Vector containing the coefficients of the trace of the spectrum.
   *
   * len_rows: Length of the rows of the flattened-matrix.
   *
   * len_cols: Length of the columns of the flattened-matrix.
   *
   * len_v:    Length of the coefficient-vector.
   *
   * Length:   Given that the trace is centered on the spectrum, Length gives the length,
   *           in pixels from the center, where we'll consider the spectrum.
   *
   * RON:      The Read Out Noise of our measurements in electrons.
   *
   * GAIN:     The Gain of the detector, in electrons/ADU.
   * 
   * mode:     Set this to 0 to apply Marsh's Algorithm (curved trace) to the data.
   *           Set this to 1 to apply Horne's Algorithm (trace paralell to the rows) to the data.
   *
   * ------------------------------------------------------------------------------
   */
  PyArg_ParseTuple(args,"OOOOiiidddddiiii",&marray,&barray,&vararray,&varray,&len_rows,&len_cols,&len_v,&Length,&RON,&GAIN,&NSigma,&S,&N,&mode,&col_min,&col_max);
  if(col_min!=0){
    col_min=col_min+1;
  }
  if(col_max!=0){
    col_max=col_max-1;
  }
  debug=0;
  m = ARRAYD(marray);                               /* We convert our PyObject struct pointers to C-vector array  */
  v = ARRAYD(varray);
  varianceimage = ARRAYD(vararray);
  bimage = ARRAYD(barray);
  real_cols=len_cols;
  real_rows=len_rows; 
  if(mode!=0)
    S=1.0;                                          /* If Horne's algorithm is used, the spacing between
						       polynomials is 1                                              */
  Length = CheckAperture(v, len_rows, len_v, S, Length);
  K=(int)(2*(int)((Length/S)+0.5)+1);               /* The number of polynomials. The idea is to form n=Length/S 
						       polynomials, rounding the decimal number, from the center.
						       Then we multiply this number by 2 (to get the total number
						       of polynomials up and down the center) and add 1 to get the
						       middle polynomial (the one corresponding to the center).      */
  if(debug!=0){
    printf("Number of polynomial to fit: %d \n",K);
    dvalue=1018.0;
  }
  C=-(S*(1.0+(((double)K-1.0)/2.0)));               /* We get the C constant for the Q coefficients. This number
						       is the polynomial number of the center.                       */

  double* Range = RangeDetector(v,C,S,len_v,K,len_rows,len_cols,range);  /* This vector containts the detected starting 
									    and final column pixels that "swims in" and "out" of the image, respectively
									    (FREED here, later). */
  if(debug!=0){
    printf("Range[0]: %f\n",Range[0]);
    printf("Range[1]: %f\n",Range[1]);
    printf("len_cols=%d, len_rows=%d \n",len_cols,len_rows);
  }
  range=Range[0];

  if(Range[0]!=-1){
    if(col_min==0 && col_max==0){
      len_cols=Range[1]-Range[0]+1;                   /* If the column presents problems, we start by creating a    
							 cutted image                                               */
    }
    else if((col_min!=0 && col_max!=0) || (col_min==0 && col_max!=0)){
      if(Range[0]>=col_min && Range[1]<=col_max){
	len_cols=Range[1]-Range[0]+1;
      }
      else if(Range[0]<col_min && Range[1]<=col_max){
	range=col_min;
	Range[0]=col_min;
	len_cols=Range[1]-Range[0]+1;
      }
      else if(Range[0]<col_min && Range[1]>col_max){
	range=col_min;
	Range[0]=col_min;
	Range[1]=col_max;
	len_cols=Range[1]-Range[0]+1;
      }
      else if(Range[0]>=col_min && Range[1]>col_max){
	Range[1]=col_max;
	len_cols=Range[1]-Range[0]+1;
      }
    }
    else if(col_min!=0 && col_max==0){
      if(Range[0]>=col_min){
	len_cols=Range[1]-Range[0]+1;
      }
      else if(Range[0]<col_min){
	range=col_min;
	Range[0]=col_min;
	len_cols=Range[1]-Range[0]+1;
      }
    }
  }
  else if((col_min!=0 && col_max!=0) || (col_min==0 && col_max!=0)){
    range=col_min;
    Range[0]=col_min;
    Range[1]=col_max;
    len_cols=Range[1]-Range[0]+1;
  }
  else if(col_min!=0 && col_max==0){
    range=col_min;
    Range[0]=col_min;
    Range[1]=len_cols-1;
    len_cols=Range[1]-Range[0]+1;
  }

  double** M = MakeArray(len_rows,len_cols);        /* Create a matrix array for our image (FREED on MarshIt.)    */
  double** Var = MakeArray(len_rows,len_cols);
	
  if(Range[0]==-1){
    for(i=0;i<len_rows;i++){                        /* If the column doesn't present problems, proceed as usual   */
      for(j=0;j<len_cols;j++){
	if(bimage[i*real_cols+j]==1){
	  M[i][j]=m[i*real_cols+j]-dvalue;
	}
	else{
	  M[i][j]=-9999;
	}
	Var[i][j]=varianceimage[i*real_cols+j];
      }
    }
  }
  else{
    for(i=0;i<len_rows;i++){                        /* If it does, we fill our cutted image                       */
      for(j=range;j<Range[1]+1;j++){
	if(bimage[i*real_cols+j]==1){
	  M[i][j-range]=m[i*real_cols+j]-dvalue;
	}
	else{
	  M[i][j-range]=-9999;
	}
	Var[i][j-range]=varianceimage[i*real_cols+j];
      }
    }
  }

  /* Start of the Marsh's Algorithm iteration: */

  pmin=MakeVector(len_cols);
  pmax=MakeVector(len_cols);
  getMinMax(pmin,pmax,v,C,S,K,len_v,len_cols,range);
  double** PF = MarshIteration(Transpose(M,len_rows,len_cols),Transpose(Var,len_rows,len_cols),pmin,pmax,v,C,S,RON,GAIN,len_v,len_rows,len_cols,N,K,mode,range,debug,NSigma,isvarimage); // FREED (here, later).
  double** PFinal = MakeArray(real_rows,real_cols);  // FREED (here, later).
  if(Range[0]==-1){
    for(i=0;i<len_rows;i++){                       
      for(j=0;j<len_cols;j++){ 
	PFinal[i][j]=PF[i][j];
      }
    }
  }
  else{
    for(i=0;i<len_rows;i++){                        /* If it does, we fill our cutted image                       */
      for(j=range;j<Range[1]+1;j++){ 
	PFinal[i][j]=PF[i][j-range];
      }
    }
    len_cols=real_cols;
  }
  FreeArray(PF,len_rows);
  gettimeofday(&tim, NULL);
  double t2=tim.tv_sec+(tim.tv_usec/1000000.0);
  if(debug!=0){
    printf("%2.2lf seconds elapsed\n", t2-t1);
  }
  // /* End of the Marsh's Algorithm iteration! */
  // 
  // /* Start of the matrix-to-vector conversion part of the code */    
  double* theArray;  // FREED (here, later).
  theArray = (double*) malloc((len_cols*len_rows)*sizeof(double));
  for(i=0;i<len_rows;i++){
    for(j=0;j<len_cols;j++){
      theArray[i*len_cols+j]=PFinal[i][j];
      if(PFinal[i][j]>1){
	theArray[i*len_cols+j]=0;
      }
    } 
  }
  FreeArray(PFinal,len_rows);
  free(Range);
  /* End of the matrix-to-vector conversion part of the code. Everything's ok up to here...now I have my theArray[i][j] array (matrix) ready to be called...*/


  /* Finally, we create a Python "Object" List that contains the P coefficients and return it back to Python */

  PyObject *lst = PyList_New(len_rows*len_cols);
  PyObject *num;
  if (!lst)
    return NULL;
  for (i = 0; i < len_rows*len_cols; i++) {
    num=PyFloat_FromDouble(theArray[i]);
    if (!num) {
      Py_DECREF(lst);
      return NULL;
    }
    PyList_SET_ITEM(lst, i, num);
  }
  free(theArray);
  PyObject *MyResult = Py_BuildValue("O",lst);
  Py_DECREF(lst);
  return MyResult;
}

/*                        [INITIALIZATION OF A METHOD]
 * ------------------------THE SOBTAINSpectrum METHOD-----------------------------
 */

static PyObject *Marsh_SObtainSpectrum(PyObject *self, PyObject *args){
  struct timeval tim;
  gettimeofday(&tim, NULL);
  int i,j,real_cols,real_rows,col_min,col_max;
  double t1=tim.tv_sec+(tim.tv_usec/1000000.0),dvalue=0;
  double *mflat,*varflat,*pflat,*bflat;
  double *v,*pmin,*pmax; // (pmin,pmax) freed here, later.
  double Length,RON,GAIN,CosmicSigma,S,C;
  int len_v,len_cols,len_rows,K,range=-1,debug,isvarimage=1;
  PyObject *parray, *marray,*varray,*vararray,*barray;
  /* 
   *--------------------------------THE DATA---------------------------------------
   * After initialization of the PyObject pointers, we wish to recover the following inputs:
   *
   * marray:   Vector of the flattened-matrix of the Echelle Spectra data.
   *
   * varray:   Vector containing the coefficients of the trace of the spectrum.
   *
   * parray:   Vector of the flattened-matrix of the Echelle Spectra P_{ij} coefficents.
   *
   * len_rows: Length of the rows of the flattened-matrix.
   *
   * len_cols: Length of the columns of the flattened-matrix.
   *
   * len_v:    Length of the coefficient-vector.
   * 
   * Length:   Aperture to be used in the obtention of the spectrum.
   *
   * RON:      The Read Out Noise of our measurements.
   * 
   * GAIN:     The Gain of the detector, in electrons/ADU.
   * 
   * S:        Spacing of the polynomials obtained in the ObtainP function.
   * 
   * CosmicSigma: Number of times sigma is multiplied by to reject Cosmic Rays.
   *
   * ------------------------------------------------------------------------------
   */
  PyArg_ParseTuple(args,"OOOOOiiidddddii",&marray,&barray,&vararray,&varray,&parray,&len_rows,&len_cols,&len_v,&Length,&RON,&GAIN,&S,&CosmicSigma,&col_min,&col_max);
  if(col_min!=0){
    col_min=col_min+1;
  }
  if(col_max!=0){
    col_max=col_max-1;
  }
  debug=0;
  pflat = ARRAYD(parray);                               /* We convert our PyObject struct pointers to C-vector array  */
  mflat = ARRAYD(marray);                               /* We convert our PyObject struct pointers to C-vector array  */
  varflat = ARRAYD(vararray);
  bflat = ARRAYD(barray);
  v = ARRAYD(varray);
  real_cols=len_cols;
  real_rows=len_rows;
  Length = CheckAperture(v, len_rows, len_v, S, Length);
  /* Start of the algorithm for the spectrum obtention (faster with the transpose!): */
	
  K=(int)(2*(int)((Length/S)+0.5)+1);               /* The number of polynomials. The idea is to form n=Length/S 
						       polynomials, rounding the decimal number, from the center.
						       Then we multiply this number by 2 (to get the total number
						       of polynomials up and down the center) and add 1 to get the
						       middle polynomial (the one corresponding to the center).      */
  C=-(S*(1.0+(((double)K-1.0)/2.0)));               /* We get the C constant for the Q coefficients. This number
						       is the polynomial number of the center.                       */
  range=0;
  double* Range = RangeDetector(v,C,S,len_v,K,len_rows,len_cols,range);                  /* This vector containts the
											    detected starting and final column  pixels that "swims in" and "out" of the image, respectively (FREED here) */
  range=Range[0];
  if(debug!=0){
    printf("Range[0]: %f\n",Range[0]);
    printf("Range[1]: %f\n",Range[1]);
    dvalue=1018.0;
  } 

  if(Range[0]!=-1){
    if(col_min==0 && col_max==0){
      len_cols=Range[1]-Range[0]+1;                   /* If the column presents problems, we start by creating a    
							 cutted image                                               */
    }
    else if((col_min!=0 && col_max!=0) || (col_min==0 && col_max!=0)){
      if(Range[0]>=col_min && Range[1]<=col_max){
	len_cols=Range[1]-Range[0]+1;
      }
      else if(Range[0]<col_min && Range[1]<=col_max){
	range=col_min;
	Range[0]=col_min;
	len_cols=Range[1]-Range[0]+1;
      }
      else if(Range[0]<col_min && Range[1]>col_max){
	range=col_min;
	Range[0]=col_min;
	Range[1]=col_max;
	len_cols=Range[1]-Range[0]+1;
      }
      else if(Range[0]>=col_min && Range[1]>col_max){
	Range[1]=col_max;
	len_cols=Range[1]-Range[0]+1;
      }
    }
    else if(col_min!=0 && col_max==0){
      if(Range[0]>=col_min){
	len_cols=Range[1]-Range[0]+1;
      }
      else if(Range[0]<col_min){
	range=col_min;
	Range[0]=col_min;
	len_cols=Range[1]-Range[0]+1;
      }
    }
  }
  else if((col_min!=0 && col_max!=0) || (col_min==0 && col_max!=0)){
    range=col_min;
    Range[0]=col_min;
    Range[1]=col_max;
    len_cols=Range[1]-Range[0]+1;
  }
  else if(col_min!=0 && col_max==0){
    range=col_min;
    Range[0]=col_min;
    Range[1]=len_cols-1;
    len_cols=Range[1]-Range[0]+1;
  }

  double** M = MakeArray(len_rows,len_cols);        /* Create a matrix array for our image (FREED on getSpectrum) */
  double** P = MakeArray(len_rows,len_cols);        /* Create a matrix array for our P's (FREED on getSpectrum)   */
  double** Var = MakeArray(len_rows,len_cols);
  if(Range[0]==-1){
    for(i=0;i<len_rows;i++){                        /* If the column doesn't present problems, proceed as usual   */
      for(j=0;j<len_cols;j++){
	if(bflat[i*real_cols+j]==1){
	  M[i][j]=mflat[i*real_cols+j]-dvalue;
	}
	else{
	  M[i][j]=-9999;
	}
	P[i][j]=pflat[i*real_cols+j];
	Var[i][j]=varflat[i*real_cols+j];
      }
    }
  }
  else{
    for(i=0;i<len_rows;i++){                        /* If it does, we fill our cutted image                       */
      for(j=range;j<Range[1]+1;j++){ 
	if(bflat[i*real_cols+j]==1){
	  M[i][j-range]=mflat[i*real_cols+j]-dvalue;
	}
	else{
	  M[i][j-range]=-9999;
	}
	P[i][j-range]=pflat[i*real_cols+j];
	Var[i][j-range]=varflat[i*real_cols+j];
      }
    }
  }

  pmin=MakeVector(len_cols);
  pmax=MakeVector(len_cols);
  getMinMax(pmin,pmax,v,C,S,K,len_v,len_cols,range);
  double** Spectrum=getSpectrum(Transpose(M,len_rows,len_cols),Transpose(Var,len_rows,len_cols),Transpose(P,len_rows,len_cols),v,len_v,pmin,pmax,RON,GAIN,len_rows,len_cols,CosmicSigma,debug,isvarimage); // FREED (here, later)
  gettimeofday(&tim, NULL);
  double t2=tim.tv_sec+(tim.tv_usec/1000000.0);
  if(debug!=0){
    printf("%2.2lf seconds elapsed\n", t2-t1);
  }
  /* End of the algorithm */

  double* theArray; // FREED (here, later)
  theArray = (double*) malloc((real_cols*3)*sizeof(double));
  for(i=0;i<3;i++){
    for(j=0;j<real_cols;j++){
      if(range==-1){
	theArray[i*real_cols+j]=Spectrum[i][j];
      }
      else if(range!=-1 && i==0){
	theArray[i*real_cols+j]=j;
      }
      else if(range!=-1 && i!=0){
	if(j<Range[0] || j>Range[1])
	  theArray[i*real_cols+j]=0;
	else
	  theArray[i*real_cols+j]=Spectrum[i][j-range];
      }
    } 
  }

  FreeArray(Spectrum,3);
  free(Range);

  /* Finally, we create a Python "Object" List that contains the P coefficients and return it back to Python */

  PyObject *lst = PyList_New(real_cols*3);
  if (!lst)
    return NULL;
  for (i = 0; i < 3*real_cols; i++) {
    PyObject *num = PyFloat_FromDouble(theArray[i]);
    if (!num) {
      Py_DECREF(lst);
      return NULL;
    }
    PyList_SET_ITEM(lst, i, num);
  }
  free(theArray);
  free(pmin); 
  free(pmax);
  PyObject *MyResult = Py_BuildValue("Oi",lst,real_cols);
  Py_DECREF(lst);
  return MyResult;
}

static PyObject *Marsh_BObtainSpectrum(PyObject *self, PyObject *args){
  struct timeval tim;
  gettimeofday(&tim, NULL);
  int i,j,real_cols,real_rows,col_min,col_max;
  double t1=tim.tv_sec+(tim.tv_usec/1000000.0),dvalue=0;
  double *mflat,*pflat,*bflat;
  double *v,*pmin,*pmax; // (pmin,pmax) freed here, later.
  double Length,RON,GAIN,CosmicSigma,S,C;
  int len_v,len_cols,len_rows,K,range=-1,debug;
  PyObject *parray, *marray,*varray,*barray;
  /* 
   *--------------------------------THE DATA---------------------------------------
   * After initialization of the PyObject pointers, we wish to recover the following inputs:
   *
   * marray:   Vector of the flattened-matrix of the Echelle Spectra data.
   *
   * varray:   Vector containing the coefficients of the trace of the spectrum.
   *
   * parray:   Vector of the flattened-matrix of the Echelle Spectra P_{ij} coefficents.
   *
   * len_rows: Length of the rows of the flattened-matrix.
   *
   * len_cols: Length of the columns of the flattened-matrix.
   *
   * len_v:    Length of the coefficient-vector.
   * 
   * Length:   Aperture to be used in the obtention of the spectrum.
   *
   * RON:      The Read Out Noise of our measurements.
   * 
   * GAIN:     The Gain of the detector, in electrons/ADU.
   * 
   * S:        Spacing of the polynomials obtained in the ObtainP function.
   * 
   * CosmicSigma: Number of times sigma is multiplied by to reject Cosmic Rays.
   *
   * ------------------------------------------------------------------------------
   */
  PyArg_ParseTuple(args,"OOOOiiidddddii",&marray,&varray,&parray,&barray,&len_rows,&len_cols,&len_v,&Length,&RON,&GAIN,&S,&CosmicSigma,&col_min,&col_max);
  if(col_min!=0){
    col_min=col_min+1;
  }
  if(col_max!=0){
    col_max=col_max-1;
  }
  debug=0;
  pflat = ARRAYD(parray);                               /* We convert our PyObject struct pointers to C-vector array  */
  mflat = ARRAYD(marray);                               /* We convert our PyObject struct pointers to C-vector array  */
  bflat = ARRAYD(barray);                               /* We convert our PyObject struct pointers to C-vector array  */
  v = ARRAYD(varray);
  real_cols=len_cols;
  real_rows=len_rows;
  Length = CheckAperture(v, len_rows, len_v, S, Length);
  /* Start of the algorithm for the spectrum obtention (faster with the transpose!): */
	
  K=(int)(2*(int)((Length/S)+0.5)+1);               /* The number of polynomials. The idea is to form n=Length/S 
						       polynomials, rounding the decimal number, from the center.
						       Then we multiply this number by 2 (to get the total number
						       of polynomials up and down the center) and add 1 to get the
						       middle polynomial (the one corresponding to the center).      */
  C=-(S*(1.0+(((double)K-1.0)/2.0)));               /* We get the C constant for the Q coefficients. This number
						       is the polynomial number of the center.                       */
  range=0;
  double* Range = RangeDetector(v,C,S,len_v,K,len_rows,len_cols,range);                  /* This vector containts the
											    detected starting and final column  pixels that "swims in" and "out" of the image, respectively (FREED here) */
  range=Range[0];
  if(debug!=0){
    printf("Range[0]: %f\n",Range[0]);
    printf("Range[1]: %f\n",Range[1]);
    dvalue=0.0;
  } 

  if(Range[0]!=-1){
    if(col_min==0 && col_max==0){
      len_cols=Range[1]-Range[0]+1;                   /* If the column presents problems, we start by creating a    
							 cutted image                                               */
    }
    else if((col_min!=0 && col_max!=0) || (col_min==0 && col_max!=0)){
      if(Range[0]>=col_min && Range[1]<=col_max){
	len_cols=Range[1]-Range[0]+1;
      }
      else if(Range[0]<col_min && Range[1]<=col_max){
	range=col_min;
	Range[0]=col_min;
	len_cols=Range[1]-Range[0]+1;
      }
      else if(Range[0]<col_min && Range[1]>col_max){
	range=col_min;
	Range[0]=col_min;
	Range[1]=col_max;
	len_cols=Range[1]-Range[0]+1;
      }
      else if(Range[0]>=col_min && Range[1]>col_max){
	Range[1]=col_max;
	len_cols=Range[1]-Range[0]+1;
      }
    }
    else if(col_min!=0 && col_max==0){
      if(Range[0]>=col_min){
	len_cols=Range[1]-Range[0]+1;
      }
      else if(Range[0]<col_min){
	range=col_min;
	Range[0]=col_min;
	len_cols=Range[1]-Range[0]+1;
      }
    }
  }
  else if((col_min!=0 && col_max!=0) || (col_min==0 && col_max!=0)){
    range=col_min;
    Range[0]=col_min;
    Range[1]=col_max;
    len_cols=Range[1]-Range[0]+1;
  }
  else if(col_min!=0 && col_max==0){
    range=col_min;
    Range[0]=col_min;
    Range[1]=len_cols-1;
    len_cols=Range[1]-Range[0]+1;
  }

  double** M = MakeArray(len_rows,len_cols);        /* Create a matrix array for our image (FREED on getSpectrum) */
  double** B = MakeArray(len_rows,len_cols);        /* Create a matrix array for our image (FREED on getSpectrum) */
  double** P = MakeArray(len_rows,len_cols);        /* Create a matrix array for our P's (FREED on getSpectrum)   */
  if(Range[0]==-1){
    for(i=0;i<len_rows;i++){                        /* If the column doesn't present problems, proceed as usual   */
      for(j=0;j<len_cols;j++){ 
	M[i][j]=mflat[i*real_cols+j]-dvalue;
	B[i][j]=bflat[i*real_cols+j];
	P[i][j]=pflat[i*real_cols+j];
      }
    }
  }
  else{
    for(i=0;i<len_rows;i++){                        /* If it does, we fill our cutted image                       */
      for(j=range;j<Range[1]+1;j++){ 
	M[i][j-range]=mflat[i*real_cols+j]-dvalue;
	B[i][j-range]=bflat[i*real_cols+j];
	P[i][j-range]=pflat[i*real_cols+j];
      }
    }
  }

  pmin=MakeVector(len_cols);
  pmax=MakeVector(len_cols);
  getMinMax(pmin,pmax,v,C,S,K,len_v,len_cols,range);
  double** Spectrum=BgetSpectrum(Transpose(M,len_rows,len_cols),Transpose(B,len_rows,len_cols),Transpose(P,len_rows,len_cols),v,len_v,pmin,pmax,RON,GAIN,len_rows,len_cols,CosmicSigma,debug); // FREED (here, later)
  gettimeofday(&tim, NULL);
  double t2=tim.tv_sec+(tim.tv_usec/1000000.0);
  if(debug!=0){
    printf("%2.2lf seconds elapsed\n", t2-t1);
  }
  /* End of the algorithm */

  double* theArray; // FREED (here, later)
  theArray = (double*) malloc((real_cols*4)*sizeof(double));
  for(i=0;i<4;i++){
    for(j=0;j<real_cols;j++){
      if(range==-1){
	theArray[i*real_cols+j]=Spectrum[i][j];
      }
      else if(range!=-1 && i==0){
	theArray[i*real_cols+j]=j;
      }
      else if(range!=-1 && i!=0){
	if(j<Range[0] || j>Range[1])
	  theArray[i*real_cols+j]=0;
	else
	  theArray[i*real_cols+j]=Spectrum[i][j-range];
      }
    } 
  }

  FreeArray(Spectrum,4);
  free(Range);

  /* Finally, we create a Python "Object" List that contains the P coefficients and return it back to Python */

  PyObject *lst = PyList_New(real_cols*4);
  if (!lst)
    return NULL;
  for (i = 0; i < 4*real_cols; i++) {
    PyObject *num = PyFloat_FromDouble(theArray[i]);
    if (!num) {
      Py_DECREF(lst);
      return NULL;
    }
    PyList_SET_ITEM(lst, i, num);
  }
  free(theArray);
  free(pmin); 
  free(pmax);
  PyObject *MyResult = Py_BuildValue("Oi",lst,real_cols);
  Py_DECREF(lst);
  return MyResult;
}

static PyObject *Marsh_SimpleExtraction(PyObject *self, PyObject *args){
  struct timeval tim;
  gettimeofday(&tim, NULL);
  int i,j,real_cols,debug=0,col_min,col_max;
  double t1=tim.tv_sec+(tim.tv_usec/1000000.0),vMin,vMax,dvalue=0;
  double *mflat;
  double *v;
  double Length,C,S;
  int len_v,len_cols,len_rows,K,range;
  PyObject *marray,*varray;
  /* 
   *--------------------------------THE DATA---------------------------------------
   * After initialization of the PyObject pointers, we wish to recover the following inputs:
   *
   * marray:   Vector of the flattened-matrix of the Echelle Spectra data.
   *
   * varray:   Vector containing the coefficients of the trace of the spectrum.
   *
   * len_rows: Length of the rows of the flattened-matrix.
   *
   * len_cols: Length of the columns of the flattened-matrix.
   *
   * len_v:    Length of the coefficient-vector.
   * 
   * Length:   Aperture (in pixels) to be taken from the center to cover the entire spectrum.
   *
   * ------------------------------------------------------------------------------
   */
  PyArg_ParseTuple(args,"OOiiidii",&marray,&varray,&len_rows,&len_cols,&len_v,&Length,&col_min,&col_max);
  if(col_min!=0){
    col_min=col_min+1;
  }
  if(col_max!=0){
    col_max=col_max-1;
  }
  range=0;
  mflat = ARRAYD(marray);                               /* We convert our PyObject struct pointers to C-vector array  */
  v = ARRAYD(varray);
  real_cols=len_cols;
  S=1.0;
  Length = CheckAperture(v, len_rows, len_v, S, Length);
  K=(int)(2*(int)((Length/S))+1);

  C=-(S*(1.0+(((double)K-1.0)/2.0)));               /* We get the C constant for the Q coefficients. This number
						       is the polynomial number of the center.                       */

  double* Range = SimpleRangeDetector(v,Length,len_v,len_rows,len_cols,range);                  /* This vector
												   containts the detected starting and final column pixels that "swims in" and "out" of the image, respectively (F) */
  range=Range[0];

  if(debug!=0){
    printf("Range[0]: %f\n",Range[0]);
    printf("Range[1]: %f\n",Range[1]);
    dvalue=1018.0;
  } 

  if(Range[0]!=-1){
    if(col_min==0 && col_max==0){
      len_cols=Range[1]-Range[0]+1;                   /* If the column presents problems, we start by creating a    
							 cutted image                                               */
    }
    else if((col_min!=0 && col_max!=0) || (col_min==0 && col_max!=0)){
      if(Range[0]>=col_min && Range[1]<=col_max){
	len_cols=Range[1]-Range[0]+1;
      }
      else if(Range[0]<col_min && Range[1]<=col_max){
	range=col_min;
	Range[0]=col_min;
	len_cols=Range[1]-Range[0]+1;
      }
      else if(Range[0]<col_min && Range[1]>col_max){
	range=col_min;
	Range[0]=col_min;
	Range[1]=col_max;
	len_cols=Range[1]-Range[0]+1;
      }
      else if(Range[0]>=col_min && Range[1]>col_max){
	Range[1]=col_max;
	len_cols=Range[1]-Range[0]+1;
      }
    }
    else if(col_min!=0 && col_max==0){
      if(Range[0]>=col_min){
	len_cols=Range[1]-Range[0]+1;
      }
      else if(Range[0]<col_min){
	range=col_min;
	Range[0]=col_min;
	len_cols=Range[1]-Range[0]+1;
      }
    }
  }
  else if((col_min!=0 && col_max!=0) || (col_min==0 && col_max!=0)){
    range=col_min;
    Range[0]=col_min;
    Range[1]=col_max;
    len_cols=Range[1]-Range[0]+1;
  }
  else if(col_min!=0 && col_max==0){
    range=col_min;
    Range[0]=col_min;
    Range[1]=len_cols-1;
    len_cols=Range[1]-Range[0]+1;
  }

  double** M = MakeArray(len_rows,len_cols);        /* Create a matrix array for our image (F)                    */
  if(Range[0]==-1){
    for(i=0;i<len_rows;i++){                        /* If the column doesn't present problems, proceed as usual   */
      for(j=0;j<len_cols;j++){ 
	M[i][j]=mflat[i*real_cols+j]-dvalue;
      }
    }
  }
  else{
    for(i=0;i<len_rows;i++){                        /* If it does, we fill our cutted image                       */
      for(j=range;j<Range[1]+1;j++){
	M[i][j-range]=mflat[i*real_cols+j]-dvalue;
      }
    }
  }

  double* pmin=MakeVector(len_cols); // (F)
  double* pmax=MakeVector(len_cols); // (F)
  getSimpleMinMax(pmin,pmax,v,len_v,len_cols,range,Length);
  SimplePixelResampling(M,pmin,pmax,len_cols);
  //         C=C-0.5;
  //         getMinMax();


  /* Start of the algorithm for the Simple Extraction Algorithm */

  double* FSpectrum=MakeVector(len_cols);  // (F)
  double Sum;
  for(i=0;i<len_cols;i++){
    Sum=0;
    vMin=pmin[i];
    vMax=pmax[i];
    for(j=vMin;j<vMax;j++){
      Sum=M[j][i]+Sum;	    
    }
    FSpectrum[i]=Sum;
  }

  FreeArray(M,len_rows);
  free(pmin);
  free(pmax);
	
  double* Spectrum=MakeVector(real_cols); // (F)
         
  if(range!=-1){
    for(i=0;i<real_cols;i++){
      if(i>=range && i<Range[1]+1){
	Spectrum[i]=FSpectrum[i-range];
      }
      else{
	Spectrum[i]=0.0;
      }
    }
  }
  else{
    for(i=0;i<real_cols;i++){
      Spectrum[i]=FSpectrum[i];
    }
  }
  free(FSpectrum);
  len_cols=real_cols;

  gettimeofday(&tim, NULL);
  double t2=tim.tv_sec+(tim.tv_usec/1000000.0);
  if(debug!=0){
    printf("%2.2lf seconds elapsed\n", t2-t1);
  }

  /* End of the algorithm */

  /* Finally, we create a Python "Object" List that contains the P coefficients and return it back to Python */

  PyObject *lst = PyList_New(len_cols);
  if (!lst)
    return NULL;
  for (i = 0; i < len_cols; i++) {
    PyObject *num = PyFloat_FromDouble(Spectrum[i]);
    if (!num) {
      Py_DECREF(lst);
      return NULL;
    }
    PyList_SET_ITEM(lst, i, num);
  }

  free(Spectrum);
  free(Range);

  PyObject *MyResult = Py_BuildValue("O",lst);
  Py_DECREF(lst);
  return MyResult;
}

static PyMethodDef MarshMethods[] = {
  {"ObtainP", Marsh_ObtainP, METH_VARARGS, "First step in the Method for optimum spectra obtention from Echelle Spectrographs: The obtention of the light fractions P."},
  {"SObtainP", Marsh_ObtainP, METH_VARARGS, "Special function same as ObtainP, but for the case when we have a variance image."},
  {"BObtainSpectrum", Marsh_BObtainSpectrum, METH_VARARGS, "Final step in the Method for optimum spectra obtention from Echelle Spectrographs: The obtention of the spectrum."},
  {"ObtainSpectrum", Marsh_ObtainSpectrum, METH_VARARGS, "Final step in the Method for optimum spectra obtention from Echelle Spectrographs: The obtention of the spectrum."},
  {"SimpleExtraction", Marsh_SimpleExtraction,METH_VARARGS, "Function that performs a simple extraction for comparison with the Optimal Extraction algorithm."},
  {"SObtainSpectrum", Marsh_ObtainSpectrum, METH_VARARGS, "Special function, same as ObtainSpectrum but when we have a variance image."},
  {"SSimpleExtraction", Marsh_SimpleExtraction,METH_VARARGS, "Special function, same as SObtainSpectrum, but when we have a variance image."},
  {NULL, NULL, 0, NULL}
};

#if PY_MAJOR_VERSION >= 3
static struct PyModuleDef moduledef = {
  PyModuleDef_HEAD_INIT,
  "Marsh",
  NULL,
  0, //no state, so re-initialization is fine
  MarshMethods,
  NULL,
  NULL,
  NULL,
  NULL
};

PyMODINIT_FUNC PyInit_Marsh(void)
{
  PyObject* module = PyModule_Create(&moduledef);
  return module;
}

#else
void initMarsh(void){
  Py_InitModule("Marsh", MarshMethods);
}
#endif


/*********************************************************************
 *          [START OF THE FUNCTIONS OF THE MARSH ALGORITHM]          *
 *********************************************************************
 */

double** Transpose(double **O,int rowsIN,int colsIN){
  int i,j;
  double **OT=MakeArray(colsIN,rowsIN);
  for(i=0;i<rowsIN;i++){
    for(j=0;j<colsIN;j++){
      OT[j][i]=O[i][j];
    }
  }
  FreeArray(O,rowsIN);
  return OT;
}
// Must GENERALIZE this
void PixelResampling(double** A,double* pmin,double* pmax,int len_cols){
  int j;
  double MinFluxFraction,MaxFluxFraction=0;
  for(j=0;j<len_cols;j++){
    if(A[j][(int)pmin[j]]!=-9999){
      MinFluxFraction=1.0-(pmin[j]-(double)((int)pmin[j])); // The fraction is 1-(RealPixel+0.5-(int)pmin[j])   
      MaxFluxFraction=pmax[j]-(double)((int)pmax[j]);       // The fraction is RealPixel+0.5-(int)pmax[j]
      A[j][(int)pmin[j]]=MinFluxFraction*A[j][(int)pmin[j]];
      A[j][(int)pmax[j]]=MaxFluxFraction*A[j][(int)pmax[j]];
    }
  }
  
}

void BPixelResampling(double** A,double** B,double* pmin,double* pmax,int len_cols){
  int j;
  double MinFluxFraction,MaxFluxFraction=0;
  for(j=0;j<len_cols;j++){
    MinFluxFraction=1.0-(pmin[j]-(double)((int)pmin[j])); // The fraction is 1-(RealPixel+0.5-(int)pmin[j])   
    MaxFluxFraction=pmax[j]-(double)((int)pmax[j]);       // The fraction is RealPixel+0.5-(int)pmax[j]
    A[j][(int)pmin[j]]=MinFluxFraction*A[j][(int)pmin[j]];
    A[j][(int)pmax[j]]=MaxFluxFraction*A[j][(int)pmax[j]];
    B[j][(int)pmin[j]]=MinFluxFraction*B[j][(int)pmin[j]];
    B[j][(int)pmax[j]]=MaxFluxFraction*B[j][(int)pmax[j]];
  }
  
}

void SimplePixelResampling(double** A,double* pmin,double* pmax,int len_cols){
  int j;
  double MinFluxFraction,MaxFluxFraction=0;
  for(j=0;j<len_cols;j++){
    MinFluxFraction=1.0-(pmin[j]-(double)((int)pmin[j])); // The fraction is 1-(RealPixel+0.5-(int)pmin[j])   
    MaxFluxFraction=pmax[j]-(double)((int)pmax[j]);       // The fraction is RealPixel+0.5-(int)pmax[j]
    A[(int)pmin[j]][j]=MinFluxFraction*A[(int)pmin[j]][j];
    A[(int)pmax[j]][j]=MaxFluxFraction*A[(int)pmax[j]][j];
  }
  
}

/*                                      [FUNCTION]
 * --------------------------------LINEAR SYSTEM SOLVER--------------------------------
 * This function recieves a GSL vector b, a matrix A and returns the solution to the
 * system Ax=b. Here, length is the number of equations of the linear system.
 * ------------------------------------------------------------------------------------
 */


void LinearSolver(double** A,double* b,double* x,int length){
  int i,j;

  gsl_matrix *AGSL = gsl_matrix_alloc(length,length);   /* We allocate memory for our GSL matrices and vectors    */
  gsl_vector *xGSL = gsl_vector_alloc(length);
  gsl_vector *bGSL = gsl_vector_alloc(length);

  for(i=0;i<length;i++){
    for(j=0;j<length;j++){
      gsl_matrix_set(AGSL,i,j,A[i][j]);               /* Set the obtained values of the C_qp matrix to "A"      */
    }
    gsl_vector_set(bGSL,i,b[i]);                       /* Set the obtained values of the X_q vector to "b"       */
  } 
  //        gsl_linalg_HH_solve(AGSL,bGSL,xGSL);                  /* Solve the system...                                    */
  gsl_vector *tau = gsl_vector_alloc(length);
  gsl_linalg_QR_decomp(AGSL,tau);
  gsl_linalg_QR_solve(AGSL,tau,bGSL,xGSL);
  for(i=0;i<length;i++){
    x[i]=gsl_vector_get(xGSL,i);                       /* Set the solution in each B_q vector term               */
  } 
}


/*                                     [FUNCTION]
 * ------------------------------------TRACE VALUE--------------------------------------
 * Function that returns the value of the position of the trace at a given X (given in
 * Matrix Coordiantes).
 * -------------------------------------------------------------------------------------
 */

double Trace(int x,double *v,int len_v,int range){
  //  int i;
  double Sum=0;
  if(range==-1){
    Sum=v[x];
    //     for(i=0;i<len_v;i++)
    //        Sum=v[(len_v-1)-i]*pow(x,(double)i)+Sum;
  }
  else{
    Sum=v[x+range];
    //     for(i=0;i<len_v;i++)
    //        Sum=v[(len_v-1)-i]*pow(x+range,(double)i)+Sum;
  }
  return Sum;
}

/*                               [FUNCTION]
 * ----------------------------Sum of the rows ------------------------------------
 * Here we compute the sum of the row values for each column. This is the "standard
 * spectrum", i.e., the spectrum one would extract with unitary weigths.
 * --------------------------------------------------------------------------------
 */

void getRowSum(double** A,double* RS,double *v,double* pmin,double* pmax,int len_cols,int len_v){
  double Sum,c;
  int i,j,vMin,vMax;
  for(j=0;j<len_cols;j++){
    Sum=0;
    vMin=pmin[j];
    vMax=pmax[j];
    for(i=vMin;i<=vMax;i++){
      c=A[j][i];
      if(c!=-9999){
        Sum=c+Sum;
      }
    }
    RS[j]=Sum;          /* We save the sum of the row (column) of the untransposed 
			   (transposed) matrix.                                    */
  }
}

/*                       [FUNCTION]
 * --------------------E (P Estimators)---------------------------------------
 * Here we compute the light fractions estimators E and save it on our matrix.
 * ---------------------------------------------------------------------------
 */

void getE(double** A,double* RS,double** E,double* pmin,double* pmax,double *v,int len_rows,int len_cols){
  double c;
  int i,j,vMin,vMax;
  for(j=0;j<len_cols;j++){
    vMin=pmin[j];
    vMax=pmax[j];
    for(i=vMin;i<=vMax;i++){
      c=A[j][i];
      if(c!=-9999){
	E[j][i]=c/RS[j]; /* Our estimations are the pixel value divided by 
			    the row sum (column sum of the untransposed (transposed)
			    matrix: E_{ji}=A_{ji}/RS_{j}                          */
      }	
    }
  }
}

/*                             [FUNCTION]
 * --------------------Sum of the rows (Var(E) Estimators)---------------------------
 * Here we compute the sum of the row values of the variances for each column.
 * ----------------------------------------------------------------------------------
 */
void getVarRowSum(double** A,double** V,double* VarRS,double* pmin,double* pmax,double RON,double GAIN,int len_cols){
  int i,j,vMin,vMax;
  double SumVar;
  for(j=0;j<len_cols;j++){
    SumVar=0;
    vMin=pmin[j];
    vMax=pmax[j];
    for(i=vMin;i<=vMax;i++){
      if(A[j][i]!=-9999)
	SumVar=V[j][i]+SumVar;       /* Each variance is taken as Poisson Noise+RON**2 (see
					the getImageVariances() and VarRevision() functions).*/
    }
    VarRS[j]=SumVar;
  }
}


/*                       [FUNCTION]
 * --------------------Var(E) (P Estimators)---------------------------------------
 * Here we compute the light fraction's variance for the estimators E and return
 * them as a matrix. Commented is the delta method.
 * ---------------------------------------------------------------------------
 */

void getVarE(double** A,double** E,double* RS,double* VarRS,double** VarE,double** V,double* pmin,double* pmax,double RON,double GAIN,int len_rows,int len_cols){
  int i,j,vMin,vMax;
  double a1=0,a2=0,c;
  //  double c;
  for(j=0;j<len_cols;j++){
    vMin=pmin[j];
    vMax=pmax[j];
    for(i=vMin;i<=vMax;i++){
      c=A[j][i];
      if(c!=-9999){
	a1=((1.0/pow(RS[j],2.0))-(2.0*c/pow(RS[j],3.0)))*(V[j][i]);
	a2=(pow(E[j][i]/RS[j],2.0))*VarRS[j];
	VarE[j][i]=a1+a2;
	//             VarE[j][i]=V[j][i]/pow(RS[j],2.0);
      }
    }
  }
} 

/*                       [FUNCTION]
 * --------------------Q Coefficients------------------------------------------
 * Here we compute the Q coefficients for our light fractions using the linear
 * interpolation method proposed by Marsh (1989).
 * ---------------------------------------------------------------------------
 */

double CalculateQ(double *v,int k,int i,int j,double C,double S,int len_v,int range){
  double d,x=Trace(j,v,len_v,range)+(S*(double)k)+C;      // We compute the position of the center of the polynomial.
  d=fabs(x-(double)i);
  if(d>= 0.5+S)
    return 0.0;
  else if(d+S <= 0.5)
    return S;
  else if(d <= 0.5){
    return ((S/2.0)+(0.5-d)-(pow(0.5-d,2.0)/(2.0*S)));
  }
  else if(d > 0.5){
    return ((0.5-d)+(pow(0.5-d,2.0)/(2.0*S))+S/2);
  }
  else
    return -999;
}

double** getQ(double *v,double C,double S,int len_v,int K,int len_cols,int range,double* pmin,double* pmax,int mode){
  int i,j,vMin,vMax,contador,k;
  contador=0;
  for(j=0;j<len_cols;j++){
    vMin=pmin[j];         // If the value extends to a lower pixel, we take it into account.
    vMax=pmax[j];         // If the value extends to a higher pixel, we take it into account.
    for(i=vMin;i<=vMax;i++){
      contador+=1;
    }
  }
  double** Q=MakeArray(K,contador);
  for(k=0;k<=(K-1);k++){
    contador=0;
    for(j=0;j<len_cols;j++){
      vMin=pmin[j];         // If the value extends to a lower pixel, we take it into account.
      vMax=pmax[j];         // If the value extends to a higher pixel, we take it into account.
      for(i=vMin;i<=vMax;i++){
	if(mode==0){
          Q[k][contador]=CalculateQ(v,k+1,i,j,C,S,len_v,range);
	}
	else{
	  Q[k][contador]=1.0;
	}
	contador+=1;
      }
    }
  }
  return Q;
}

/*                       [FUNCTION]
 * --------------------C Coefficients------------------------------------------
 * Here we compute the C coefficients for the obtention of the C_p matrix
 * proposed by Marsh (1989).
 * ---------------------------------------------------------------------------
 */

double CalculateC(double** A,double** Q,double** VarE,int m, int l, int n, int k,double** J,double* pmin,double* pmax,int len_cols){
  int i,j,vMin,vMax,contador;
  double TotalSum,Qk,Ql,Power;
  TotalSum=0,Qk=0,Ql=0; 
  contador=0;
  for(j=0;j<len_cols;j++){
    vMin=pmin[j];
    vMax=pmax[j];                // If the value extends to a higher pixel, we take it into account.
    Power=J[n+m-2][j];
    for(i=vMin;i<=vMax;i++){
      Qk=Q[k-1][contador];
      Ql=Q[l-1][contador];
      /* We save computation time if either of Q is 0: */
      if((Qk!=0.0) && (Ql!=0.0) && (A[j][i]!=-9999)){  
        /* If it isn't, we continue the calculation of our sum. */
	TotalSum=((Qk*Ql*Power)/(VarE[j][i]))+TotalSum;
      }

      contador+=1;
    }
  }
  return TotalSum;
}

void getC(double** A,double** Q,double** VarE,double** C_qp,double** J,double* pmin,double* pmax,int K,int N,int len_cols){
  int p,q,m,l,n,k,dummyq;
  double value=0;
  p=0,q=0,dummyq=1;
  for(k=1;k<=K;k++){
    for(n=1;n<=N;n++){
      for(l=1;l<=K;l++){
	for(m=1;m<=N;m++){
          p=N*(l-1)+m;
          q=N*(k-1)+n;
          if(dummyq<q){              // When we change of "q", we go to p=q, because the matrix is symmetric
	    m=n;                    // and therefore all values can be obtained from the upper ones.
	    l=k;
	    p=N*(l-1)+m;
	    q=N*(k-1)+n;
          }
          value=CalculateC(A,Q,VarE,m,l,n,k,J,pmin,pmax,len_cols);
          C_qp[p-1][q-1]=value;
          C_qp[q-1][p-1]=value;
          dummyq=q;
	}
      }
    }
  }
}


/*                       [FUNCTION]
 * --------------------X Coefficients------------------------------------------
 * Here we compute the X coefficients for the obtention of the X_q vector
 * proposed by Marsh (1989).
 * ---------------------------------------------------------------------------
 */

double CalculateX(double** A,double** Q,double** E,double** VarE,int n, int k,double** J,double* pmin,double* pmax,int len_cols){
  int i,j,vMin,vMax,contador;
  double TotalSum,Qk,Power;
  TotalSum=0,Qk=0; 
  contador=0;
  for(j=0;j<len_cols;j++){
    vMin=pmin[j];         // If the value extends to a lower pixel, we take it into account.
    vMax=pmax[j];         // If the value extends to a higher pixel, we take it into account.
    Power=J[n-1][j];
    for(i=vMin;i<=vMax;i++){
      Qk=Q[k-1][contador];
      if(Qk!=0.0 && (A[j][i]!=-9999)){                      // We save computation time if either of Q is 0.
        TotalSum=((E[j][i]*Qk*Power)/(VarE[j][i]))+TotalSum;
      }
      contador+=1;
    }
  }
  return TotalSum;
}

void getX(double** A,double** Q,double** E,double** VarE,double* X,double** J,double* pmin,double* pmax,int K,int N,int len_cols){
  int q,n,k;
  q=0;
  for(k=1;k<=K;k++){
    for(n=1;n<=N;n++){
      q=N*(k-1)+n;
      X[q-1]=CalculateX(A,Q,E,VarE,n,k,J,pmin,pmax,len_cols);
    }
  }
}

/*                       [FUNCTION]
 * --------------------P Coefficients------------------------------------------
 * Here we compute the light fractions P and return them as a matrix
 * ---------------------------------------------------------------------------
 */

void getP(double* B,double** Q,double** P,double** J,double* pmin,double* pmax,int K,int N,int len_rows,int len_cols,int mode){
  double Sum,Gkj;
  int i,j,vMin,vMax,contador,n,k,q;
  contador=0,Sum=0,Gkj=0;
  if(mode==0){
    for(j=0;j<len_cols;j++){
      vMin=pmin[j];         // If the value extends to a lower pixel, we take it into account.
      vMax=pmax[j];         // If the value extends to a higher pixel, we take it into account.
      for(i=vMin;i<=vMax;i++){
	for(k=1;k<=K;k++){
	  for(n=1;n<=N;n++){
	    q=N*(k-1)+n;
	    Gkj=(B[q-1]*J[n-1][j])+Gkj;     /* We obtain the G_{kj} coefficient for a given k and j.*/
	  }
          Sum=(Q[k-1][contador]*Gkj)+Sum;                        /* Therms for the P_{ij} sum (eq (6), Marsh, 1989).     */
          Gkj=0;
	}
	P[j][i]=Sum;
	contador+=1;
	Sum=0;
      }
    }
  }
  else{
    for(j=0;j<len_cols;j++){
      vMin=pmin[j];         // If the value extends to a lower pixel, we take it into account.
      vMax=pmax[j];         // If the value extends to a higher pixel, we take it into account.
      k=1;
      for(i=vMin;i<=vMax;i++){
	for(n=1;n<=N;n++){
	  q=N*(k-1)+n;
	  Gkj=(B[q-1]*J[n-1][j])+Gkj;     /* We obtain the G_{kj} coefficient for a given k and j.*/
	}
	P[j][i]=Gkj;
	Gkj=0;
	k++;
      }
    }
  }
}

void VarRevision(double** A,double** VarImage,double** V,double* RS,double** P,double* pmin,double* pmax,int len_cols,double RON,double GAIN,int isvarimage){
  int i,j,vMin,vMax;
  double c;
  for(j=0;j<len_cols;j++){
    vMin=pmin[j];         // If the value extends to a lower pixel, we take it into account.
    vMax=pmax[j];         // If the value extends to a higher pixel, we take it into account.
    for(i=vMin;i<=vMax;i++){
      c=A[j][i];
      if(c!=-9999){
	if(isvarimage==0){
	  V[j][i]=pow(RON/GAIN,2.0)+fabs(RS[j]*P[j][i])/GAIN;
	}
	else{
	  V[j][i]=VarImage[j][i];
	}
      }
    }
  }
}

int OutlierRejection(double** A,double** V,double** P,double** E,double** VarE,double* RS,double* pmin,double* pmax,double NSigma,int len_cols,int debug){
  int i,j,vMin,vMax,counter=0;
  double TotalSum=0,SigmaSquared,Ratio;
  SigmaSquared=pow(NSigma,2.0);
  for(j=0;j<len_cols;j++){
    vMin=pmin[j];
    vMax=pmax[j];
    for(i=vMin;i<=vMax;i++){
      if(A[j][i]!=-9999){
	Ratio=pow(A[j][i]-RS[j]*P[j][i],2.0)/V[j][i];
	TotalSum=(pow(P[j][i]-E[j][i],2.0)/VarE[j][i])+TotalSum;
	if(Ratio>=SigmaSquared){
	  A[j][i]=-9999;
	  counter++;
	}
      }
    }
  }
  if(debug!=0){
    printf("Number of outliers found: %d \n",counter);
    printf("Chi-squared of the fit: %f \n",TotalSum);
  }
  if(counter==0){
    return 0;
  }
  return 1;
}
/*                           [FUNCTION]
 * ------------------------ Marsh's algorithm----------------------------------
 * Here we obtain images (analogous to the image of pixel count values, M) for
 * the following:
 * 
 * E_ij:      Matrix image (same ij coordiantes as M) containing the estimated light
 *            fractions (recall this is an estimate of P_ij, the model of this light
 *            fraction.
 *
 * Var(E_ij): Matrix image (same ij coordiantes as M) containing the estimated light
 *            fraction's variances.
 *
 * P_ij:      Matrix image (same ij coordiantes as M) containing the modeled light
 *            fractions.
 * ---------------------------------------------------------------------------
 */

double** MarshIteration(double** M,double** VarImage,double* pmin,double* pmax,double *v,double C,double S,double RON,double GAIN,int len_v,int len_rows,int len_cols,int N,int K,int mode,int range,int debug,double NSigma, int isvarimage){
  int Iteration=1,IterNum=0;
  double** P=MakeArray(len_cols,len_rows); // FREED on Transpose.
  double** V=MakeArray(len_cols,len_rows); // FREED here, later.
  double** E=MakeArray(len_cols,len_rows); // FREED here, later.
  double** VarE=MakeArray(len_cols,len_rows); // FREED here, later.
  double** C_qp=MakeArray(N*K,N*K);                              /* Matrix that contains the coefficients for the linear system  (FREED here, later) */
  double** J=MakeArray(2*(N-1)+1,len_cols); // FREED here, later.

  // All from here to the end, FREED, later.

  double*  B=MakeVector(N*K);                                    /* Vector with the solutions to the linear system.              */
  double*  X=MakeVector(N*K);                                    /* Vector with the coefficients for the linear system. Given 
								    these elements, we wish to solve the C_qp B_q=X_q system     */
  double* RS=MakeVector(len_cols);                               /* Vector that saves the sum of pixel values of the rows
								    (columns) of the transposed (untransposed) M matrix.         */
  double* VarRS=MakeVector(len_cols);                            /* Vector that saves the sum of pixel values variances of the rows
								    (columns) of the transposed (untransposed) M matrix          */
  PixelResampling(M,pmin,pmax,len_cols);
  double** Q=getQ(v,C,S,len_v,K,len_cols,range,pmin,pmax,mode);                                             /* First we obtain the Q matrix (note: independant of fit)      */
  getJ(J,N,len_cols);
  getImageVariances(M,V,VarImage,pmin,pmax,RON,GAIN,len_cols,isvarimage);
  getRowSum(M,RS,v,pmin,pmax,len_cols,len_v);
  getE(M,RS,E,pmin,pmax,v,len_rows,len_cols);
  getVarRowSum(M,V,VarRS,pmin,pmax,RON,GAIN,len_cols);
  getVarE(M,E,RS,VarRS,VarE,V,pmin,pmax,RON,GAIN,len_rows,len_cols); 
  while(Iteration != 0){                                 /* If Iteration=1, we continue iterating!                       */
    IterNum+=1;
    if(debug!=0){
      printf("------------------- \n");
      printf("Iteration number %d \n",IterNum);
      printf("------------------- \n");
      printf("Obtaining the profiles...\n");             /* We obtain Var(E_ij), the variances of E_ij on each pixel     */
    }
    getC(M,Q,VarE,C_qp,J,pmin,pmax,K,N,len_cols);        /* We obtain C_qp, the matrix with the C coefficients for
							    the fit                                       */
    getX(M,Q,E,VarE,X,J,pmin,pmax,K,N,len_cols);                                /* We obtain X_q, the vector with the coefficients for the fit  */
    LinearSolver(C_qp,X,B,N*K);                          /* Solve the linear system C_qp*B_q=X_q, obtaining the
							    B_q vector                                                   */
    getP(B,Q,P,J,pmin,pmax,K,N,len_rows,len_cols,mode);                                       /* We obtain P, the model light fractions on each pixel         */
    Renormalize(P,pmin,pmax,len_rows,len_cols);
    VarRevision(M,VarImage,V,RS,P,pmin,pmax,len_cols,RON,GAIN,isvarimage);
    getVarRowSum(M,V,VarRS,pmin,pmax,RON,GAIN,len_cols);
    getVarE(M,E,RS,VarRS,VarE,V,pmin,pmax,RON,GAIN,len_rows,len_cols); 
    if(IterNum!=1)
      Iteration=OutlierRejection(M,V,P,E,VarE,RS,pmin,pmax,NSigma,len_cols,debug);
  }
  FreeArray(E,len_cols);                                        /* Free our vectors and matrices                                 */
  FreeArray(VarE,len_cols);
  FreeArray(C_qp,N*K);
  FreeArray(Q,K);
  FreeArray(J,2*(N-1)+1);
  FreeArray(M,len_cols);
  FreeArray(VarImage,len_cols);
  FreeArray(V,len_cols);
  free(B);
  free(X);
  free(RS);
  free(VarRS);
  free(pmin);
  free(pmax);
  if(debug!=0){
    printf("Going out of MarshIteration...\n");
  }
  return Transpose(P,len_cols,len_rows);
}

void Renormalize(double** P,double* pmin,double* pmax,int len_rows,int len_cols){
  int vMin,vMax,j,i;
  double PartialSum;
  for(j=0;j<len_cols;j++){
    PartialSum=0;
    vMin=pmin[j];
    vMax=pmax[j];
    for(i=vMin;i<=vMax;i++){
      if(P[j][i]<0)
	P[j][i]=0;
      PartialSum=P[j][i]+PartialSum;
    }
    for(i=vMin;i<=vMax;i++){
      P[j][i]=P[j][i]/PartialSum;
    }
    
  }
} 

double** getSpectrum(double** A,double** VarImage,double** P,double *v,int len_v,double* pmin,double* pmax,double RON,double GAIN,int len_rows,int len_cols,double CosmicSigma,int debug,int isvarimage){
  int j;
  int Iteration=1,IterNum=0;
  double** W=MakeArray(len_cols,len_rows); // FREED, later
  double** V=MakeArray(len_cols,len_rows); // FREED, later
  double* F=MakeVector(len_cols); // FREED, later
  double* VarF=MakeVector(len_cols); // FREED, later
  double* RSW=MakeVector(len_cols); // FREED, later
  PixelResampling(A,pmin,pmax,len_cols);
  getRowSum(A,F,v,pmin,pmax,len_cols,len_v);
  while(Iteration == 1){   
    VarRevision(A,VarImage,V,F,P,pmin,pmax,len_cols,RON,GAIN,isvarimage);
    IterNum+=1;
    if(debug!=0){
      printf("------------------- \n");
      printf("Iteration number %d \n",IterNum);
      printf("------------------- \n");
      printf("Obtaining the spectrum...\n");
    }
    getRowSumW(A,V,P,pmin,pmax,RSW,len_cols);    /* We obtain the denominator for the W_i (see Marsh (1989), eq. (4)) */
    getW(A,P,V,RSW,W,pmin,pmax,len_cols);        /* We obtain the weights to estimate the fluxes.                     */
    getF(W,A,F,pmin,pmax,len_cols);
    getVarF(A,W,V,VarF,pmin,pmax,len_cols);
    if(CosmicSigma==0)
      Iteration=0;
    else
      Iteration=CosmicRayRejection(A,V,P,F,VarF,CosmicSigma,pmin,pmax,len_cols,debug);
  }
  double** Spectrum=MakeArray(3,len_cols); // FREED, on the main code
  for(j=0;j<len_cols;j++){
    Spectrum[0][j]=(double)j;
    Spectrum[1][j]=F[j];
    Spectrum[2][j]=(double)1/(VarF[j]);
  }
  FreeArray(W,len_cols);
  FreeArray(V,len_cols);
  FreeArray(P,len_cols);
  FreeArray(A,len_cols);
  FreeArray(VarImage,len_cols);
  free(F);
  free(VarF);
  // free(OutlierDetector);
  free(RSW);
  return Spectrum;
}

double** BgetSpectrum(double** A,double** B,double** P,double *v,int len_v,double* pmin,double* pmax,double RON,double GAIN,int len_rows,int len_cols,double CosmicSigma,int debug){
  int j;
  int Iteration=1,IterNum=0,isvarimage=0;
  double** W=MakeArray(len_cols,len_rows); // FREED, later
  double** V=MakeArray(len_cols,len_rows); // FREED, later
  double* F=MakeVector(len_cols); // FREED, later
  double* BF=MakeVector(len_cols); // REED
  double* VarF=MakeVector(len_cols); // FREED, later
  double* RSW=MakeVector(len_cols); // FREED, later
  BPixelResampling(A,B,pmin,pmax,len_cols);
  getRowSum(A,F,v,pmin,pmax,len_cols,len_v);
  while(Iteration == 1){   
    VarRevision(A,A,V,F,P,pmin,pmax,len_cols,RON,GAIN,isvarimage);
    IterNum+=1;
    if(debug!=0){
      printf("------------------- \n");
      printf("Iteration number %d \n",IterNum);
      printf("------------------- \n");
      printf("Obtaining the spectrum...\n");
    }
    getRowSumW(A,V,P,pmin,pmax,RSW,len_cols);    /* We obtain the denominator for the W_i (see Marsh (1989), eq. (4)) */
    getW(A,P,V,RSW,W,pmin,pmax,len_cols);        /* We obtain the weights to estimate the fluxes.                     */
    getF(W,A,F,pmin,pmax,len_cols);
    getVarF(A,W,V,VarF,pmin,pmax,len_cols);
    if(CosmicSigma==0)
      Iteration=0;
    else
      Iteration=BCosmicRayRejection(A,B,V,P,F,VarF,CosmicSigma,pmin,pmax,len_cols,debug);
  }
  getBF(W,B,BF,pmin,pmax,len_cols);
  double** Spectrum=MakeArray(4,len_cols); // FREED, on the main code
  for(j=0;j<len_cols;j++){
    Spectrum[0][j]=(double)j;
    Spectrum[1][j]=F[j];
    Spectrum[2][j]=(double)1/(VarF[j]);
    Spectrum[3][j]=BF[j];
  }
  FreeArray(W,len_cols);
  FreeArray(V,len_cols);
  FreeArray(P,len_cols);
  FreeArray(A,len_cols);
  FreeArray(B,len_cols);
  free(F);
  free(BF);
  free(VarF);
  // free(OutlierDetector);
  free(RSW);
  return Spectrum;
}

int CosmicRayRejection(double** A,double** V,double** P,double* F,double* VarF,double CosmicSigma,double* pmin,double* pmax,int len_cols,int debug){
  int i,j,vMin,vMax,ii=0,jj=0,counter=0,detection=0;
  double Ratio,DummyRatio=0,count,fluxcount,countsigma,fluxcountsigma;
  for(j=0;j<len_cols;j++){
    vMin=pmin[j];
    vMax=pmax[j];
    for(i=vMin;i<=vMax;i++){
      if(A[j][i]!=-9999){
	count=A[j][i];
	countsigma=sqrt(V[j][i]);
	fluxcount=F[j]*P[j][i];
	fluxcountsigma=sqrt(VarF[j])*P[j][i];
	if(count>=fluxcount){
	  Ratio=(count-CosmicSigma*countsigma)-(fluxcount+CosmicSigma*fluxcountsigma);
	  if(Ratio>0){
	    detection=1;
	  }
	}
	else{
	  Ratio=(fluxcount-CosmicSigma*fluxcountsigma)-(count+CosmicSigma*countsigma);
	  if(Ratio>0){
	    detection=1;
	  }
	}
	if(detection==1){
	  if(DummyRatio<Ratio){
	    ii=i;
	    jj=j;
	    DummyRatio=Ratio;
	    counter=1;
	    detection=0;
	  }
	}
      }
    }
  }

  if(counter==0){
    return 0;
  }
  else{
    if(debug!=0){
      printf("Cosmic ray found at column %d row %d \n",jj,ii);
      printf("Counts=%f, Light fraction (P)=%f , Total flux in column=%f (var=%f) , Variance=%f \n",A[jj][ii],P[jj][ii],F[jj],VarF[jj],V[jj][ii]);
    }
    A[jj][ii]=-9999;
    return 1;
  }

}

int BCosmicRayRejection(double** A,double** B,double** V,double** P,double* F,double* VarF,double CosmicSigma,double* pmin,double* pmax,int len_cols,int debug){
  int i,j,vMin,vMax,ii=0,jj=0,counter=0,detection=0;
  double Ratio,DummyRatio=0,count,fluxcount,countsigma,fluxcountsigma;
  for(j=0;j<len_cols;j++){
    vMin=pmin[j];
    vMax=pmax[j];
    for(i=vMin;i<=vMax;i++){
      if(A[j][i]!=-9999){
	count=A[j][i];
	countsigma=sqrt(V[j][i]);
	fluxcount=F[j]*P[j][i];
	fluxcountsigma=sqrt(VarF[j])*P[j][i];
	if(count>=fluxcount){
	  Ratio=(count-CosmicSigma*countsigma)-(fluxcount+CosmicSigma*fluxcountsigma);
	  if(Ratio>0){
	    detection=1;
	  }
	}
	else{
	  Ratio=(fluxcount-CosmicSigma*fluxcountsigma)-(count+CosmicSigma*countsigma);
	  if(Ratio>0){
	    detection=1;
	  }
	}
	if(detection==1){
	  if(DummyRatio<Ratio){
	    ii=i;
	    jj=j;
	    DummyRatio=Ratio;
	    counter=1;
	    detection=0;
	  }
	}
      }
    }
  }

  if(counter==0){
    return 0;
  }
  else{
    if(debug!=0){
      printf("Cosmic ray found at column %d row %d \n",jj,ii);
      printf("Counts=%f, Light fraction (P)=%f , Total flux in column=%f (var=%f) , Variance=%f \n",A[jj][ii],P[jj][ii],F[jj],VarF[jj],V[jj][ii]);
    }
    A[jj][ii]=-9999;
    B[jj][ii]=-9999;
    return 1;
  }

}

void getImageVariances(double** A,double** V,double** VarImage,double* pmin,double* pmax,double RON,double GAIN,int len_cols,int isvarimage){
  int vMin,vMax,j,i;
  for(j=0;j<len_cols;j++){
    vMin=pmin[j];         // If the value extends to a lower pixel, we take it into account.
    vMax=pmax[j];         // If the value extends to a higher pixel, we take it into account.
    for(i=vMin;i<=vMax;i++){
      //           V[j][i]=fabs(A[j][i]/GAIN)+pow(RON/GAIN,2.0);
      if(isvarimage==0){
	V[j][i]=1.0;
      }
      else{
	if(A[j][i]!=-9999){
	  V[j][i]=VarImage[j][i];
	}
      }
    }
    
  }
} 

void getRowSumW(double** A,double** V,double** P,double* pmin,double* pmax,double* RSW,int len_cols){
  int i,j,vMin,vMax;
  double Sum;
  for(j=0;j<len_cols;j++){
    Sum=0;
    vMin=pmin[j];         // If the value extends to a lower pixel, we take it into account.
    vMax=pmax[j];         // If the value extends to a higher pixel, we take it into account.
    for(i=vMin;i<=vMax;i++){
      if(A[j][i]!=-9999)
	Sum=(pow(P[j][i],2.0)/V[j][i])+Sum;
    }
    RSW[j]=Sum;
  }
}

void getW(double** A,double** P,double** V,double* RSW,double** W,double* pmin,double* pmax,int len_cols){
  int vMin,vMax,j,i;
  double PartialSum;
  for(j=0;j<len_cols;j++){
    PartialSum=0;
    vMin=pmin[j];         // If the value extends to a lower pixel, we take it into account.
    vMax=pmax[j];         // If the value extends to a higher pixel, we take it into account.
    for(i=vMin;i<=vMax;i++){
      if(A[j][i]!=-9999)
	W[j][i]=(P[j][i]/V[j][i])/RSW[j];
    }
    
  }
} 

void getF(double** W,double** A,double* F,double* pmin,double* pmax,int len_cols){
  int vMin,vMax,j,i;
  double PartialSum;
  for(j=0;j<len_cols;j++){
    PartialSum=0;
    vMin=pmin[j];         // If the value extends to a lower pixel, we take it into account.
    vMax=pmax[j];         // If the value extends to a higher pixel, we take it into account.
    for(i=vMin;i<=vMax;i++){
      if(A[j][i]!=-9999)
	PartialSum=(W[j][i]*A[j][i])+PartialSum;
    }
    F[j]=PartialSum;
  }
} 

void getBF(double** W,double** B,double* BF,double* pmin,double* pmax,int len_cols){
  int vMin,vMax,j,i;
  double PartialSum;
  for(j=0;j<len_cols;j++){
    PartialSum=0;
    vMin=pmin[j];         // If the value extends to a lower pixel, we take it into account.
    vMax=pmax[j];         // If the value extends to a higher pixel, we take it into account.
    for(i=vMin;i<=vMax;i++){
      if(B[j][i]!=-9999)
	PartialSum=(W[j][i]*B[j][i])+PartialSum;
    }
    BF[j]=PartialSum;
  }
} 
 
void getVarF(double** A,double** W,double** V,double* VarF,double* pmin,double* pmax,int len_cols){
  int vMin,vMax,j,i;
  double PartialSum;
  for(j=0;j<len_cols;j++){
    PartialSum=0;
    vMin=pmin[j];         // If the value extends to a lower pixel, we take it into account.
    vMax=pmax[j];         // If the value extends to a higher pixel, we take it into account.
    for(i=vMin;i<=vMax;i++){
      if(A[j][i]!=-9999)
	PartialSum=(pow(W[j][i],2.0)*V[j][i])+PartialSum;
    }
    VarF[j]=PartialSum;
  }
} 

void getMinMax(double* pmin,double* pmax,double *v,double C, double S,int K,int len_v,int len_cols,int range){
  int j;
  for(j=0;j<len_cols;j++){
    pmin[j]=((Trace(j,v,len_v,range)+C)+S)+0.5;
    pmax[j]=((Trace(j,v,len_v,range)+C)+(double)(K)*S)+0.5;
  }
}

void getSimpleMinMax(double* pmin,double* pmax,double *v,int len_v,int len_cols,int range,double Length){
  int j;
  for(j=0;j<len_cols;j++){
    pmin[j]=Trace(j,v,len_v,range)-Length+0.5;
    pmax[j]=Trace(j,v,len_v,range)+Length+0.5;
  }
}

void getJ(double** J,int N,int len_cols){
  int ulimit=2*(N-1)+1,power;
  int j;
  for(power=0;power<ulimit;power++){
    for(j=0;j<len_cols;j++){
      J[power][j]=pow((double)j,power);
    }
  }
}

double** MakeArray(int rows, int columns){
  int i,j;
  double** theArray;
  theArray = (double**) malloc(rows*sizeof(double*));
  for(i=0;i<rows;i++)
    theArray[i] = (double*) malloc(columns*sizeof(double));

  /* Fill the array with zeroes (i.e. we clean it) */

  for(i=0;i<rows;i++){
    for(j=0;j<columns;j++){
      theArray[i][j]=0.0;
    }
  }

  return theArray;
}

int CheckAperture(double* center, int len_rows, int len_cen, double S, int CentralAperture){
  int uprow,lowrow;
  int K = (int)((int)((CentralAperture/S)+0.5)); // This is the number of polynomials at each side of the central aperture
  int i;
  int InitialAperture = CentralAperture;
  for(i=0;i<len_cen;i++){
    uprow = (int)(center[i]+S*(K+1.0)+0.5);
    lowrow = (int)(center[i]-S*(K+1.0)-0.5);
    if(uprow>=len_rows){
      CentralAperture = CentralAperture - 1;
      K = (int)((int)((CentralAperture/S)+0.5));
      lowrow = (int)(center[i]-S*(K+1.0)-0.5);
      i = i-1;
    }   
    if(lowrow<0 && uprow>=len_rows){
      CentralAperture = CentralAperture - 1;
    }
    if(lowrow<0 && uprow<len_rows){
      CentralAperture = CentralAperture - 1;
      K = (int)((int)((CentralAperture/S)+0.5));
      i = i-1;
    }
    // If the aperture doesn't stabilize to a value above 1,
    // then it can't be stabilized and we trust the user input value:
    if(CentralAperture<1){
      CentralAperture = InitialAperture;
      break;
    }
  }
  if(CentralAperture != InitialAperture){
    printf("MARSH CODE: Changed central aperture to %d\n",CentralAperture);
  }
  return CentralAperture;
}

double* MakeVector(int nelements){
  double* Vector;
  int j;
  Vector = (double*) malloc(nelements*sizeof(double));

  for(j=0;j<nelements;j++){
    Vector[j]=0.0;
  }
  return Vector;
}

void FreeArray(double** theArray,int rows){
  int i;
  for(i=0;i<rows;i++){
    free(theArray[i]);
  }
  free(theArray);
}

double* RangeDetector(double *v,double C,double S,int len_v,int K,int len_rows,int len_cols,int range){
  int j,s=0,uplim=0,downlim=0;
  double* ReturningVector = MakeVector(2);
  ReturningVector[1]=len_cols;
  if(Trace(0,v,len_v,range)>Trace(len_cols/2,v,len_v,range) && ((int)(((Trace(0,v,len_v,range)+C)+(double)(K)*S)+0.5)>=len_rows || (int)(((Trace(len_cols-1,v,len_v,range)+C)+(double)(K)*S)+0.5)>=len_rows)){              // The problem is on an upper trace
    for(j=0;j<len_cols;j++){
      uplim=((Trace(j,v,len_v,range)+C)+(double)(K)*S)+0.5;
      if( uplim<len_rows && s==0 ){
	s++;
	ReturningVector[0]=j;
      }
      else if( uplim==len_rows && s>0 ){
	ReturningVector[1]=j-1;
	break;
      }
    }
  }
  else if(Trace(0,v,len_v,range)<Trace(len_cols/2,v,len_v,range) && ((int)(((Trace(0,v,len_v,range)+C)+S)+0.5)<0 || (int)(((Trace(len_cols-1,v,len_v,range)+C)+S)+0.5)<0)){                                        // The problem is on an lower trace
    for(j=0;j<len_cols;j++){
      downlim=((Trace(j,v,len_v,range)+C)+S)+0.5;
      if( downlim==0 && s==0 ){
	s++;
	ReturningVector[0]=j;
      }
      else if( downlim<0 && s>0 ){
	ReturningVector[1]=j-1;
	break;
      }
    }
  }
  else{
    ReturningVector[0]=-1;
  }
  return ReturningVector;
}

double* SimpleRangeDetector(double *v,double Length,int len_v,int len_rows,int len_cols,int range){
  int j,s=0,uplim=0,downlim=0;
  double* ReturningVector = MakeVector(2);
  ReturningVector[1]=len_cols;
  if(Trace(0,v,len_v,range)>Trace(len_cols/2,v,len_v,range) && ((int)((Trace(0,v,len_v,range)+Length)+0.5)>=len_rows || (int)((Trace(len_cols-1,v,len_v,range)+Length)+0.5)>=len_rows)){              // The problem is on an upper trace
    for(j=0;j<len_cols;j++){
      uplim=(Trace(j,v,len_v,range)+Length)+0.5;
      if( uplim<len_rows && s==0 ){
	s++;
	ReturningVector[0]=j;
      }
      else if( uplim==len_rows && s>0 ){
	ReturningVector[1]=j-1;
	break;
      }
    }
  }
  else if(Trace(0,v,len_v,range)<Trace(len_cols/2,v,len_v,range) && ((int)((Trace(0,v,len_v,range)-Length)+0.5)<0 || (int)((Trace(len_cols-1,v,len_v,range)-Length)+0.5)<0)){                                        // The problem is on an lower trace
    for(j=0;j<len_cols;j++){
      downlim=(Trace(j,v,len_v,range)-Length)+0.5;
      if( downlim==0 && s==0 ){
	s++;
	ReturningVector[0]=j;
      }
      else if( downlim<0 && s>0 ){
	ReturningVector[1]=j-1;
	break;
      }
    }
  }
  else{
    ReturningVector[0]=-1;
  }
  return ReturningVector;
}

