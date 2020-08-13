/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2019 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2019 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2019 Total, S.A
 * Copyright (c) 2019-     GEOSX Contributors
 * All right reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include <Python.h>
#include <numpy/arrayobject.h>
#include <iostream>
#include <sstream>
#include "common/Logger.hpp"

void
init_numpy()
{
  import_array();
}

int
main( int argc, char ** argv )
{
  GEOSX_ERROR_IF( argc != 5, "The arguments are: path module x y" );
  GEOSX_LOG( "Testing c++ python bindings..." );

  // Initialize Python, setup the user-defined search path
  std::cout << "Initializing python..." << std::endl;
  Py_Initialize();
  PyObject * sysPath = PySys_GetObject( (char *)"path" );
  PyObject * modPath = PyString_FromString( (char *)argv[1] );
  PyList_Insert( sysPath, 0, modPath );
  init_numpy();

  // Load the targeted module and function
  PyObject * pModule = PyImport_ImportModule( (char *)argv[2] );
  if( pModule == NULL )
  {
    PyErr_Print();
    GEOSX_ERROR( "Error in loading python module" );
  }

  // Call the function designed to print and return modified string
  PyObject * pFunc =
    PyObject_GetAttrString( pModule, (char *)"printStringFromPython" );
  PyObject * pValue = Py_BuildValue( "(z)", (char *)"blah-blah-blah" );
  PyObject * pResult = PyObject_CallObject( pFunc, pValue );
  GEOSX_LOG( "Returned string: " << PyString_AsString( pResult ) );

  // Cast an array of doubles into a numpy array, modify the values in python
  const int numpy_size = 10;
  double numpy_array[numpy_size];
  npy_intp numpy_dims = numpy_size;
  for( int ii = 0; ii < numpy_size; ++ii )
  {
    numpy_array[ii] = (double)ii;
  }

  PyObject * pArray =
    PyArray_SimpleNewFromData( 1, &numpy_dims, NPY_DOUBLE, (void *)numpy_array );
  PyObject * pFuncArgsA = PyTuple_New( 1 );
  PyTuple_SetItem( pFuncArgsA, 0, pArray );
  pFunc = PyObject_GetAttrString( pModule, (char *)"modifyNumpyArray" );
  pResult = PyObject_CallObject( pFunc, pFuncArgsA );

  // Print modified values
  std::cout << "Modified Array:\n[";
  for( int ii; ii < numpy_size; ++ii )
  {
    std::cout << numpy_array[ii] << ", ";
  }
  std::cout << "]" << std::endl;

  // On the fly, generate a function that will add two floats
  std::stringstream buf;
  buf << "def add( n1 , n2 ) :" << std::endl
      << "    return n1+n2" << std::endl;
  PyObject * pCompiledFn = Py_CompileString( buf.str().c_str(), "", Py_file_input );
  assert( pCompiledFn != NULL );

  // Convert into a module
  PyObject * pModuleB = PyImport_ExecCodeModule( (char *)"test_mod_b", pCompiledFn );
  PyObject * pAddFn = PyObject_GetAttrString( pModuleB, "add" );
  PyObject * pPosArgsB = PyTuple_New( 2 );
  PyObject * pVal1 = PyFloat_FromDouble( atof( argv[3] ) );
  PyObject * pVal2 = PyFloat_FromDouble( atof( argv[4] ) );
  PyTuple_SetItem( pPosArgsB, 0, pVal1 );
  PyTuple_SetItem( pPosArgsB, 1, pVal2 );

  // Create a dict for folding keyword arguments
  PyObject * pKywdArgs = PyDict_New();
  PyObject * pResultB = PyObject_Call( pAddFn, pPosArgsB, pKywdArgs );
  PyObject * pResultRepr = PyObject_Repr( pResultB );
  std::cout << "Python add result:" << std::endl;
  std::cout << argv[3] << "+" << argv[4] << "="
            << PyString_AsString( pResultRepr ) << std::endl;

  // Clean up
  // Note: Tuples have ownership of their children, so no need to doulble decref
  // std::cout << "Cleaning up..." << std::endl;
  // Py_DecRef( pResultRepr );
  // Py_DecRef( pResultB );
  // Py_DecRef( pKywdArgs );
  // Py_DecRef( pPosArgsB );
  // Py_DecRef( pAddFn );
  // Py_DecRef( pModule );
  // Py_DecRef( pCompiledFn );
  // Py_DecRef( pFuncArgsA );
  // Py_DecRef( pResult );
  // Py_DecRef( pValue );
  // Py_DecRef( pFunc );
  // Py_DecRef( pModule );
  // Py_DecRef( modPath );
  // Py_DecRef( sysPath );

  Py_Finalize();
  std::cout << "Done!" << std::endl;
}
