#include <iostream>
#include <cstdint>
#include <cassert>
#include <vector>

#define PY_SSIZE_T_CLEAN
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION

#include <python3.6m/Python.h>

#include "cpi_ndarray.hpp"
#include "cpi_tuple.hpp"
#include "cpi_ndarray.hpp"

#include "cdf_tree_main.h"

static std::vector<CDFTree<float>>* all_data = nullptr;


extern "C" {
static PyObject * insert_sample(PyObject *self, PyObject *args) {
    (void)self;
    int index;
    PyObject *data;

    if (!PyArg_ParseTuple(args, "iO", &index, &data)) {
        std::cerr << "Cannot read inptut" << std::endl;
        return NULL;
    }
    if (all_data == nullptr) {
        PyErr_SetString(PyExc_RuntimeError, "CDFtree wants to write into unallocated memory");
        return NULL;
    }
    try {
        if (index < 0) {
            throw std::runtime_error("CDFtree index can be only positive number");
        }
        if (all_data->size() <= static_cast<unsigned>(index)) {
            all_data->resize(index+1);
        }

        NpyArray<float, 1> input_array(data, false);
        for (unsigned i = 0; i < input_array.dim_sizes[0]; ++i) {
            float sample = input_array.unsafe_get(i);
            (*all_data)[index].insert_sample(sample);
        }

        Py_INCREF(Py_None);
        return Py_None;
    } catch(std::runtime_error& e) {
        PyErr_SetString(PyExc_RuntimeError, e.what());
        return NULL;
    }
}


static PyObject * sample_to_cdf(PyObject *self, PyObject *args) {
    (void)self;
    int index;
    PyObject *data;
    bool insitu = true;

    if (!PyArg_ParseTuple(args, "iO|b", &index, &data, &insitu))
        return NULL;

    if (all_data == nullptr) {
        PyErr_SetString(PyExc_RuntimeError, "CDFtree wants to write into unallocated memory");
        return NULL;
    }
    if (index < 0) {
        PyErr_SetString(PyExc_RuntimeError, "CDFtree index can be only positive number");
        return NULL;
    }
    if (all_data->size() <= static_cast<unsigned>(index)) {
        PyErr_SetString(PyExc_RuntimeError, "CDFtree wants to use unallocated tree [possibly bad index?]");
        return NULL;
    }

    try { 
        NpyArray<float, 1> input_array(data, false);
        auto& tree = (*all_data)[index];

        if (insitu) {
            for (unsigned i = 0; i < input_array.dim_sizes[0]; ++i) {
                float& sample = input_array.unsafe_get(i);
                sample = tree.search_CDF(sample);
            }
            return input_array.pass_to_python();
        }
        else {
            NpyArray<float, 1> output_data(INIT::EMPTY, input_array.dim_sizes[0]);

            for (unsigned i = 0; i < input_array.dim_sizes[0]; ++i) {
                float sample = input_array.unsafe_get(i);
                output_data.unsafe_get(i) = tree.search_CDF(sample);
            }
            return output_data.pass_to_python();
        }
    } catch (std::runtime_error& e) {
        PyErr_SetString(PyExc_RuntimeError, e.what());
        return NULL;
    }
}

static PyObject * search_element_by_cdf(PyObject *self, PyObject *args) {
    (void)self;
    PyObject *data; bool insitu = true; int index;

    if (!PyArg_ParseTuple(args, "iO|b", &index, &data, &insitu))
        return NULL;

    if (all_data == nullptr) {
        PyErr_SetString(PyExc_RuntimeError, "CDFtree wants to use unallocated memory");
        return NULL;
    }
    if (index < 0) {
        PyErr_SetString(PyExc_RuntimeError, "CDFtree index can be only positive number");
        return NULL;
    }
    if (all_data->size() <= static_cast<unsigned>(index)) {
        PyErr_SetString(PyExc_RuntimeError, "CDFtree wants to use unallocated tree [possibly bad index?]");
        return NULL;
    }

    try {
        NpyArray<float, 1> input_array(data, false);
        auto& tree = (*all_data)[index];

        if (insitu) {
            for (unsigned i = 0; i < input_array.dim_sizes[0]; ++i) {
                float& sample = input_array.safe_get(i);
                sample = tree.inverse_search_CDF(sample);
            }
            return input_array.pass_to_python();
        } else {
            NpyArray<float, 1> output_data(INIT::EMPTY, input_array.dim_sizes[0]);

            for (unsigned i = 0; i < input_array.dim_sizes[0]; ++i) {
                float sample = input_array.unsafe_get(i);
                output_data.unsafe_get(i) = tree.inverse_search_CDF(sample);
            }
            return output_data.pass_to_python();
        }
    } catch (std::runtime_error& e) {
        PyErr_SetString(PyExc_RuntimeError, e.what());
        return NULL;
    }
}



static PyObject * init_memory(PyObject *self, PyObject *args) {
    (void)args; (void)self;
    all_data = new std::vector<CDFTree<float>>();
    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject * free_memory(PyObject *self, PyObject *args) {
    (void)args; (void)self;
    delete all_data;
    all_data = nullptr;
    Py_INCREF(Py_None);
    return Py_None;
}


static const char* doc = NULL;

static PyMethodDef Methods[] = {
    {"insert_sample", insert_sample, METH_VARARGS, "doc"},
    {"sample_to_cdf", sample_to_cdf, METH_VARARGS, "doc"},
    {"search_element_by_cdf", search_element_by_cdf, METH_VARARGS, "doc"},
    {"free_memory", free_memory, METH_VARARGS, "doc"},
    {"init_memory", init_memory, METH_VARARGS, "doc"},
    {NULL, NULL, 0, NULL}        /* Sentinel */
};

static struct PyModuleDef module = {
    PyModuleDef_HEAD_INIT,
    "cdftree",   /* name of module */
    doc, /* module documentation, may be NULL */
    -1,       /* size of per-interpreter state of the module,
                 or -1 if the module keeps state in global variables. */
    Methods
};


PyMODINIT_FUNC PyInit_libcdftree(void)
{
    PyObject *m;

    m = PyModule_Create(&module);
    if (m == NULL)
        return NULL;

    import_array();
    if (PyErr_Occurred()) return NULL;
    /* you can create some objects here <errors...> */
    return m;
}

} // extern "C"




