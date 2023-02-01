#ifdef OUTPUTGRAPHMODE
#include "Outputfile.hpp"

int *Plot::initNumpy(void) {
    import_array();
    return NULL;
}

void Plot::initPython() {
    Py_Initialize();
    initNumpy();
    PyRun_SimpleString("import sys");
    std::string tmp1 = "\"sys.path.append('";
    std::string tmp2 =  WORKDIRECTORY;
    std::string tmp3 = "')\"";
    std::string all = tmp1.append(tmp2);
    all = all.append(tmp3);
    std::cout << all << std::endl;
    PyRun_SimpleString(all.c_str());
}

void Plot::plot_zeta(int n, vvmArray & myArray) {
    // initPython();

    Py_Initialize();

    if (!Py_IsInitialized()) {
        return;
    }

    initNumpy();
    PyRun_SimpleString("import sys");
    std::string tmp1 = "sys.path.append('";
    std::string tmp2 =  WORKDIRECTORY;
    std::string tmp3 = "')";
    std::string all = tmp1.append(tmp2);
    all = all.append(tmp3);
    PyRun_SimpleString(all.c_str());            // append the python script directory

    PyObject *pDict = nullptr;
    PyObject *pModule = nullptr;
    pModule = PyImport_ImportModule("plotCallByCPlusPlus");
    pDict = PyModule_GetDict(pModule);
    npy_intp Dims[2] = {nx, nz};
    PyObject *PyArray  = PyArray_SimpleNewFromData(2, Dims, NPY_DOUBLE, myArray.zeta);
    PyObject *ArgArray = PyTuple_New(1);
    PyTuple_SetItem(ArgArray, 0, PyArray);
    PyObject *pFuncFive = PyDict_GetItemString(pDict, "zetaPlot");
    PyObject_CallObject(pFuncFive, ArgArray);
    return;
}
#endif