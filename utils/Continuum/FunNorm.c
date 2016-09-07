#include <Python.h>
#include <arrayobject.h>
#define ARRAYD(p) ((double *) (((PyArrayObject *)p)->data))


static PyObject *FunNorm_Rell(PyObject *self,PyObject *args){

int i,j,n_datos,u=1,fin;
double pref,prel,postl,postf,m,n;


double *CLam;
PyObject *PyLam;

double *CFlux;
PyObject *PyFlux;


PyArg_ParseTuple(args,"OOi",&PyLam,&PyFlux,&n_datos);
CLam = ARRAYD(PyLam);
CFlux = ARRAYD(PyFlux);

if(CFlux[n_datos-1]==-1){
		j=n_datos-2;
		while(j>=0){
			if(CFlux[j]!=-1){
				CFlux[n_datos-1]=CFlux[j];
				break;
			}
			j=j-1;

		}
	}



if(CFlux[0]==-1){
		j=0;
		while(j>=0){
			if(CFlux[j]!=-1){
				CFlux[0]=CFlux[j];
				break;
			}
			j++;

		}
	}


	
i=0;
double* theArray;
theArray = (double*) malloc((n_datos)*sizeof(double));
while(i<n_datos){
	j=i;
	u=1;
	while(u!=0){
		if(CFlux[i]<0){
			if(CFlux[i-1]>0){
				pref=CFlux[i-1];
				prel=CLam[i-1];
				
			}
			if(CFlux[i+1]>0){
				postf=CFlux[i+1];
				postl=CLam[i+1];
				m=(postf-pref)/(postl-prel);
				n=postf-m*postl;
				u=0;
				fin=i+1;
			}
			i=i+1;
			
		}
		else{
			u=0;		
		}

	}
	i=j;
	if(CFlux[i]<0){
		while(i<fin){
			CFlux[i]=CLam[i]*m+n;
			i++;
		}

	}
	i=j;
        theArray[i] = CFlux[i];
        i++;
}

PyObject *lst = PyList_New(n_datos);
if(!lst)
return NULL;
for(i=0;i<n_datos;i++){
	PyObject *num = PyFloat_FromDouble(theArray[i]);
        if(!num){
           Py_DECREF(lst);
           return NULL;
        }
	PyList_SET_ITEM(lst,i,num); 
}
free(theArray);
PyObject *MyResult = Py_BuildValue("O",lst);
Py_DECREF(lst);
return MyResult;

}





static PyMethodDef FunNormMethods[] = {
{"Rell", FunNorm_Rell, METH_VARARGS, "¿que hace esto?"},
{NULL, NULL, 0, NULL}
};

void initFunNorm(void){
(void) Py_InitModule("FunNorm", FunNormMethods);
}


