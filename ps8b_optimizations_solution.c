#include <Python.h>
#include <structmember.h>

#include <time.h>
#include <stdlib.h>

//static PyObject *NoChildException;

/*
 *####################
 *#                  #
 *#  ResistantVirus  #
 *#                  #
 *####################
 */
typedef struct {
    PyObject_HEAD
    double max_birth_prob, clear_prob, mut_prob;
    PyObject *resistances;
} ResistantVirus;

static void
ResistantVirus_dealloc(ResistantVirus *self);

static PyObject *
ResistantVirus_new(PyTypeObject *type, PyObject *args, PyObject *kwds);

static int
ResistantVirus_init(ResistantVirus *self, PyObject *args, PyObject *kwds);

static int
ResistantVirus_isResistantTo(ResistantVirus *self, PyObject *drug);

static PyObject *
ResistantVirus_reproduce(ResistantVirus *self,
    double pop_density, PyObject *active_drugs);

static PyTypeObject ResistantVirus_type = {
    PyObject_HEAD_INIT(NULL)
    0,                                        /*ob_size*/
    "ps8b.ResistantVirus",                    /*tp_name*/
    sizeof(ResistantVirus),                   /*tp_basicsize*/
    0,                                        /*tp_itemsize*/
    (destructor)ResistantVirus_dealloc,       /*tp_dealloc*/
    0,                                        /*tp_print*/
    0,                                        /*tp_getattr*/
    0,                                        /*tp_setattr*/
    0,                                        /*tp_compare*/
    0,                                        /*tp_repr*/
    0,                                        /*tp_as_number*/
    0,                                        /*tp_as_sequence*/
    0,                                        /*tp_as_mapping*/
    0,                                        /*tp_hash*/
    0,                                        /*tp_call*/
    0,                                        /*tp_str*/
    0,                                        /*tp_getattro*/
    0,                                        /*tp_setattro*/
    0,                                        /*tp_as_buffer*/
    Py_TPFLAGS_DEFAULT,                       /*tp_flags*/
    "",                                       /*tp_doc*/
    0,                                        /*tp_traverse*/
    0,                                        /*tp_clear*/
    0,                                        /*tp_richcompare*/
    0,                                        /*tp_weaklistoffset*/
    0,                                        /*tp_iter*/
    0,                                        /*tp_iternext*/
    0,                                        /*tp_methods*/
    0,                                        /*tp_members*/
    0,                                        /*tp_getset*/
    0,                                        /*tp_base*/
    0,                                        /*tp_dict*/
    0,                                        /*tp_descr_get*/
    0,                                        /*tp_descr_set*/
    0,                                        /*tp_dictoffset*/
    (initproc)ResistantVirus_init,            /*tp_init*/
    0,                                        /*tp_alloc*/
    ResistantVirus_new,                       /*tp_new*/
};

static void
ResistantVirus_dealloc(ResistantVirus *self)
{
    Py_XDECREF(self->resistances);
	self->ob_type->tp_free((PyObject *)self);
}

static PyObject *
ResistantVirus_new(PyTypeObject *type, PyObject *args, PyObject *kwds)
{
    ResistantVirus *self = NULL;
    self = (ResistantVirus *)type->tp_alloc(type, 0);
    return (PyObject *)self;
}

static int
ResistantVirus_init(ResistantVirus *self, PyObject *args, PyObject *kwds)
{
	static char *kwlist[] = {"maxBirthProb", "clearProb",
	    "resistances", "mutProb", NULL};
	PyObject *r, *tmp=NULL;

	if (! PyArg_ParseTupleAndKeywords(args, kwds, "ddOd", kwlist,
			&self->max_birth_prob, &self->clear_prob, &r, &self->mut_prob))
        return -1;
	Py_INCREF(r);
    tmp = self->resistances;
    self->resistances = r;
    if (tmp) Py_DECREF(tmp);

    return 0;
}

static int
ResistantVirus_isResistantTo(ResistantVirus *self, PyObject *drug)
{
    PyObject *item = PyDict_GetItem(self->resistances, drug);

    if (item != NULL)
        if (PyObject_IsTrue(item))
            return 0;
    return -1;
}

static PyObject *
ResistantVirus_reproduce(ResistantVirus *self,
    double pop_density, PyObject *active_drugs)
{
    PyObject *new_dict, *key, *value, *result, *args;
    Py_ssize_t pos = 0;
    int i, length = PyList_Size(active_drugs);
    double r = (double)rand()/RAND_MAX;

    for (i = 0; i < length; ++i) {
        value = PyList_GetItem(active_drugs, i);
        if (ResistantVirus_isResistantTo(self, value) < 0) {
            //PyErr_SetString(NoChildException, "");
            //return NULL;
            Py_RETURN_NONE;
        }
    }
    if (r < (self->max_birth_prob * (1 - pop_density))) {
        new_dict = PyDict_Copy(self->resistances);
        while (PyDict_Next(new_dict, &pos, &key, &value)) {
            r = (double)rand()/RAND_MAX;
            if (r < self->mut_prob) {
                if (PyObject_IsTrue(value))
                    PyDict_SetItem(new_dict, key, Py_False);
                else
                    PyDict_SetItem(new_dict, key, Py_True);
            }
        }
        args = Py_BuildValue("ddOd", self->max_birth_prob,
            self->clear_prob, new_dict, self->mut_prob);
        Py_DECREF(new_dict);
        result = PyObject_Call((PyObject *) &ResistantVirus_type, args, NULL);
        Py_DECREF(args);

        return result;
    }

    //PyErr_SetString(NoChildException, "");
    //return NULL;
    Py_RETURN_NONE;
}

/*
 *####################
 *#                  #
 *#  TreatedPatient  #
 *#                  #
 *####################
 */
typedef struct {
    PyObject_HEAD
    int max_pop;
    PyObject *viruses, *drugs;
} TreatedPatient;

static void
TreatedPatient_dealloc(TreatedPatient *self);

static PyObject *
TreatedPatient_new(PyTypeObject *type, PyObject *args, PyObject *kwds);

static int
TreatedPatient_init(TreatedPatient *self, PyObject *args, PyObject *kwds);

static PyObject *
TreatedPatient_getTotalPop(TreatedPatient *self);

static PyObject *
TreatedPatient_addPrescription(TreatedPatient *self, PyObject *args);

static PyObject *
TreatedPatient_update(TreatedPatient *self);

static PyMethodDef TreatedPatient_methods[] = {
    {"getTotalPop", (PyCFunction)TreatedPatient_getTotalPop,
        METH_NOARGS, ""},
    {"addPrescription", (PyCFunction)TreatedPatient_addPrescription,
        METH_VARARGS, ""},
    {"update", (PyCFunction)TreatedPatient_update, METH_NOARGS, ""},
    {NULL}  /* Sentinel */
};

static PyTypeObject TreatedPatient_type = {
    PyObject_HEAD_INIT(NULL)
    0,                                        /*ob_size*/
    "ps8b.TreatedPatient",                    /*tp_name*/
    sizeof(TreatedPatient),                   /*tp_basicsize*/
    0,                                        /*tp_itemsize*/
    (destructor)TreatedPatient_dealloc,       /*tp_dealloc*/
    0,                                        /*tp_print*/
    0,                                        /*tp_getattr*/
    0,                                        /*tp_setattr*/
    0,                                        /*tp_compare*/
    0,                                        /*tp_repr*/
    0,                                        /*tp_as_number*/
    0,                                        /*tp_as_sequence*/
    0,                                        /*tp_as_mapping*/
    0,                                        /*tp_hash*/
    0,                                        /*tp_call*/
    0,                                        /*tp_str*/
    0,                                        /*tp_getattro*/
    0,                                        /*tp_setattro*/
    0,                                        /*tp_as_buffer*/
    Py_TPFLAGS_DEFAULT,                       /*tp_flags*/
    "",                                       /*tp_doc*/
    0,                                        /*tp_traverse*/
    0,                                        /*tp_clear*/
    0,                                        /*tp_richcompare*/
    0,                                        /*tp_weaklistoffset*/
    0,                                        /*tp_iter*/
    0,                                        /*tp_iternext*/
    TreatedPatient_methods,                   /*tp_methods*/
    0,                                        /*tp_members*/
    0,                                        /*tp_getset*/
    0,                                        /*tp_base*/
    0,                                        /*tp_dict*/
    0,                                        /*tp_descr_get*/
    0,                                        /*tp_descr_set*/
    0,                                        /*tp_dictoffset*/
    (initproc)TreatedPatient_init,            /*tp_init*/
    0,                                        /*tp_alloc*/
    TreatedPatient_new,                       /*tp_new*/
};

static void
TreatedPatient_dealloc(TreatedPatient *self)
{
    Py_XDECREF(self->viruses);
    Py_XDECREF(self->drugs);
    self->ob_type->tp_free((PyObject *)self);
}

static PyObject *
TreatedPatient_new(PyTypeObject *type, PyObject *args, PyObject *kwds)
{
    TreatedPatient *self = NULL;
    self = (TreatedPatient *)type->tp_alloc(type, 0);
    return (PyObject *)self;
}

static int
TreatedPatient_init(TreatedPatient *self, PyObject *args, PyObject *kwds)
{
    static char *kwlist[] = {"viruses", "maxPop", NULL};
    PyObject *v, *tmp=NULL;

    if (! PyArg_ParseTupleAndKeywords(args, kwds, "Oi", kwlist,
        &v, &self->max_pop))
        return -1;
    Py_INCREF(v);
    tmp = self->viruses;
    self->viruses = v;
    if (tmp) Py_DECREF(tmp);

    self->drugs = PyList_New(0);

    return 0;
}

static PyObject *
TreatedPatient_getTotalPop(TreatedPatient *self)
{
    return PyInt_FromSsize_t(PyList_Size(self->viruses));
}

static PyObject *
TreatedPatient_addPrescription(TreatedPatient *self, PyObject *args)
{
    PyObject *drug;

    if (! PyArg_ParseTuple(args, "O", &drug))
        return NULL;
    if (! PySequence_Contains(self->drugs, drug))
        PyList_Append(self->drugs, drug);

    Py_RETURN_NONE;
}

static PyObject *
TreatedPatient_update(TreatedPatient *self)
{
    PyObject *v, *obj, *last_obj;
    int i, length = (int)PyList_Size(self->viruses);
    double density, r;

    for (i = 0; i < length; ++i) {
        r = (double)rand()/RAND_MAX;
        obj = PyList_GetItem(self->viruses, i);
        if (r < ((ResistantVirus *)obj)->clear_prob) {
            --length;
            last_obj = PyList_GetItem(self->viruses, length);
            Py_INCREF(last_obj);
            PyList_SetItem(self->viruses, i, last_obj);
            PySequence_DelItem(self->viruses, length);
            --i;
        }
    }

    density = (double)PyList_Size(self->viruses) / self->max_pop;

    for (i = 0; i < length; ++i) {
        obj = PyList_GetItem(self->viruses, i);
        v = ResistantVirus_reproduce((ResistantVirus *)obj,
            density, self->drugs);
        if (v != Py_None)
            PyList_Append(self->viruses, v);
        Py_DECREF(v);
    }

    return PyInt_FromSsize_t(PyList_Size(self->viruses));
}

/*
 *################
 *#              #
 *#    Module    #
 *#              #
 *################
 */
static PyMethodDef ps8b_methods[] = {
    {NULL}  /* Sentinel */
};

#ifndef PyMODINIT_FUNC	/* declarations for DLL import/export */
#define PyMODINIT_FUNC void
#endif
PyMODINIT_FUNC
initps8b(void)
{
    PyObject* m;

    if (PyType_Ready(&ResistantVirus_type) < 0)
        return;
    if (PyType_Ready(&TreatedPatient_type) < 0)
        return;
    m = Py_InitModule3("ps8b", ps8b_methods, "Problem Set 8 classes");
    if (m == NULL)
        return;

    srand(time(NULL));

    //NoChildException = PyErr_NewException("ps8b.NoChildException", NULL, NULL);
    //Py_INCREF(NoChildException);
    //PyModule_AddObject(m, "NoChildException", NoChildException);

    Py_INCREF(&ResistantVirus_type);
    PyModule_AddObject(m, "ResistantVirus", (PyObject *)&ResistantVirus_type);

    Py_INCREF(&TreatedPatient_type);
    PyModule_AddObject(m, "TreatedPatient", (PyObject *)&TreatedPatient_type);
}
