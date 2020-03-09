%{
#include "array.h"
%}

/* Allow Python sequences to be passed as arrays */
%typemap(in) const array& {
    PyArrayObject* pyarray = (PyArrayObject*) PyArray_ContiguousFromObject($input, PyArray_DOUBLE, 1, 1);
    if(pyarray == NULL) {
        PyErr_SetString(PyExc_TypeError, "Expected a sequence of floats.");
        return NULL;
    }

    int n = pyarray->dimensions[0];
    $1 = new array(n);
    for(int i = 0; i < n; i++)
        (*$1)[i] = *(double *)(pyarray->data + i*pyarray->strides[0]);
    Py_DECREF(pyarray);
}
%typemap(freearg) const array& {
    delete $1;
}
%typemap(typecheck) const array& {
    if(!PySequence_Check($input)) {     // Make sure we have a sequence
        printf("Not a Python sequence.\n");
        $1 = 0;
    }
    else {
        int n = PySequence_Size($input);
        $1 = 1;
        /* Make sure each element of the sequence can be cast to float */
        for(int i = 0; i < n; i++) {
            bool okay = false;
            PyObject* obj = PySequence_GetItem($input, i);
            PyObject* floatobj = PyNumber_Float(obj);
            Py_DECREF(obj);
            if(floatobj != NULL)
                Py_DECREF(floatobj);
            else {
                $1 = 0;
                break;
            }
        }
    }
}

/* Return numpy arrays */
%typemap(out) array {
    int n = $1.size();
    npy_intp dims[1];
    dims[0] = (npy_intp) n;
    $result = PyArray_SimpleNew(1, dims, NPY_DOUBLE);
    memcpy(((PyArrayObject*) $result)->data, $1.data(), n*sizeof(real));
}
//%typemap(out) const array& {
//    int n = $1->size();
//    npy_intp dims[1];
//    dims[0] = (npy_intp) n;
//    $result = PyArray_SimpleNew(1, dims, NPY_DOUBLE);
//    memcpy(((PyArrayObject*) $result)->data, $1->data(), n*sizeof(real));
//}

%typemap(in,numinputs=0) array& OUTPUT {
    $1 = new array();
}
%typemap(argout) array& OUTPUT {
    int n = $1->size();
    npy_intp dims[1] = { (npy_intp) n };
    PyObject* arrayobj = PyArray_SimpleNew(1, dims, NPY_DOUBLE);
    memcpy(((PyArrayObject*) arrayobj)->data, $1->data(), n*sizeof(real));
    $result = SWIG_Python_AppendOutput($result, arrayobj);
}
%typemap(freearg) array& OUTPUT {
    delete $1;
}
