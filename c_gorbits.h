

static PyObject *
all_orbits(PyObject *self, PyObject *args)
{
    Bitvec src;
    const char *s_arg;
    if (!PyArg_ParseTuple(args, "s", &s_arg))
        return NULL;

    assert(strlen(s_arg)==NBITS);
    BV_FROMSTRING(&src, s_arg);

    fill_gauge(src);

    PyObject *items = PyList_New(0);
    int i;
    for(i=0; i<NGAUGE; i++)
        list_append(items, all_targets[i]);
    return items;
}




