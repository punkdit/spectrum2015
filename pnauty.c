
#include "Python.h"

#ifdef _LARGEFILE_SOURCE
#undef _LARGEFILE_SOURCE
#endif


#include "nauty/nausparse.h"    /* which includes nauty.h */



DYNALLSTAT(int,lab,lab_sz);
DYNALLSTAT(int,ptn,ptn_sz);
DYNALLSTAT(int,orbits,orbits_sz);
static DEFAULTOPTIONS_SPARSEGRAPH(options);
static statsblk stats;
static sparsegraph sg;   /* Declare sparse graph structure */
int verts, edges, degree; // regular graph, fixed degree

static PyObject *
init_graph(PyObject *self, PyObject *args)
{
    int i, j, m;

    if (!PyArg_ParseTuple(args, "iii", &verts, &edges, &degree))
        return NULL;

    options.writeautoms = TRUE;
    options.defaultptn = FALSE;

    /* Initialise sparse graph structure. */

    SG_INIT(sg);

    m = SETWORDSNEEDED(verts);
    nauty_check(WORDSIZE, m, verts, NAUTYVERSIONID);

    DYNALLOC1(int, lab, lab_sz, verts, "malloc");
    DYNALLOC1(int, ptn, ptn_sz, verts, "malloc");
    DYNALLOC1(int, orbits, orbits_sz, verts, "malloc");

    for(i=0; i<verts; i++)
    {
        lab[i] = i;
        ptn[i] = 1;
    }

    /* SG_ALLOC makes sure that the v,d,e fields of a sparse graph
       structure point to arrays that are large enough.  This only
       works if the structure has been initialised. */

    assert(verts*degree == edges);
    SG_ALLOC(sg, verts, 2*edges, "malloc");

    sg.nv = verts;              /* Number of vertices */
    sg.nde = 2*edges;           /* Number of directed edges */

    // sg.d is array of vertex degrees
    // sg.e is array of vertecies
    // sg.v[i] indexes the i'th vertecies neighbours:
    //     starting at sg.e[sg.v[i]] to sg.e[sg.v[i]+sg.d[i]-1].

    for(i=0; i<verts; i++)
    {
        sg.d[i] = degree;
        sg.v[i] = degree*i;
        for(j=0; j<degree; j++)
        {
            assert(degree*i + j < 2*edges);
            sg.e[degree*i + j] = -1;
        }
    }

//    for (i = 0; i < n; ++i)
//    {
//        sg.v[i] = 2*i;
//        sg.d[i] = 2;
//        sg.e[2*i] = (i+n-1)%n;      /* edge i->i-1 */
//        sg.e[2*i+1] = (i+n+1)%n;    /* edge i->i+1 */
//    }
//
    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject *
set_partition(PyObject *self, PyObject *args)
{
    int m;

    int i, vert, partition;

    if (!PyArg_ParseTuple(args, "iii", &i, &vert, &partition))
        return NULL;

    assert(0<=i && i<verts);
    assert(0<=vert && vert<verts);
    assert(0<=partition && partition<=1);

    lab[i] = vert;
    ptn[i] = partition;

    Py_INCREF(Py_None);
    return Py_None;
}



static PyObject *
add_edge(PyObject *self, PyObject *args)
{
    int m;

    int src, tgt, i;

    if (!PyArg_ParseTuple(args, "iii", &src, &tgt, &i))
        return NULL;

    assert(0<=src && src<verts);
    assert(0<=tgt && tgt<verts);
    assert(0<=i && i<degree);

    assert(src*degree + i < 2*edges);
    assert(sg.e[src*degree + i]==-1);

    sg.e[src*degree + i] = tgt;

    Py_INCREF(Py_None);
    return Py_None;
}


static PyObject *
search(PyObject *self, PyObject *args)
{
    if (!PyArg_ParseTuple(args, ""))
        return NULL;

    int i, j;
    for(i=0; i<verts; i++)
    {
        for(j=0; j<degree; j++)
            assert(sg.e[degree*i + j]>=0); // missing an edge
    }

    //printf("Generators\n");
    options.writeautoms = FALSE;
    sparsenauty(&sg, lab, ptn, orbits, &options, &stats, NULL);

    PyObject *p_orbits;
    p_orbits = PyList_New(0);
    for(i=0; i<verts; i++)
    {
        PyObject *p_i = PyInt_FromLong(orbits[i]);
        PyList_Append(p_orbits, p_i);
        Py_DECREF(p_i);
    }

    printf("pnauty.search: Automorphism group size = ");
    writegroupsize(stdout, stats.grpsize1, stats.grpsize2);
    printf("\n");

//    Py_INCREF(Py_None);
//    return Py_None;
    return p_orbits;
}


static PyMethodDef PNautyMethods[] = 
{
    {"init_graph",  init_graph, METH_VARARGS, "builtin"},
    {"add_edge",  add_edge, METH_VARARGS, "builtin"},
    {"set_partition",  set_partition, METH_VARARGS, "builtin"},
    {"search",  search, METH_VARARGS, "builtin"},
    {NULL, NULL, 0, NULL}        /* Sentinel */
};


PyMODINIT_FUNC
initpnauty(void)
{
    (void) Py_InitModule("pnauty", PNautyMethods);
}


