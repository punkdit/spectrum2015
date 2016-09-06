
#include "Python.h"
#include <assert.h>

#ifdef _LARGEFILE_SOURCE
#undef _LARGEFILE_SOURCE
#endif


static PyObject *
bracket(PyObject *self, PyObject *args)
{
    const char *a, *b;
    int an, bn, n;
    int count, i;
    const int *A, *B;
    PyObject *result;

    if (!PyArg_ParseTuple(args, "s#s#", &a, &an, &b, &bn))
        return NULL;

    assert(an==bn);
    assert(an%8==0);
    n = an/8;

    A = (const int*)a;
    B = (const int*)b;

//    printf("A = ");
//    for(i=0; i<2*n; i++)
//        printf("%d ", A[i]);
//    printf("\n");
//
//    printf("B = ");
//    for(i=0; i<2*n; i++)
//        printf("%d ", B[i]);
//    printf("\n");

    count = 0;
    for(i=0; i<n; i++)
    {
//        printf("%d %d %d %d\n", A[2*i],B[2*i+1] , A[2*i+1],B[2*i]);
        count += A[2*i]*B[2*i+1] + A[2*i+1]*B[2*i];
    }

    count %= 2;
//    printf("result = %d\n", count);

    result = PyLong_FromLong(count);
    return result;
}


static PyMethodDef CBracketMethods[] = 
{
    {"bracket",  bracket, METH_VARARGS, "builtin"},
    {NULL, NULL, 0, NULL}        /* Sentinel */
};


PyMODINIT_FUNC
initcbracket(void)
{
    (void) Py_InitModule("cbracket", CBracketMethods);
}


