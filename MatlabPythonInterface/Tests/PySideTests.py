from __future__ import print_function
import numpy as np

# Various constants used to test the conversion of python types back to Matlab types
int_value = 1
float_value = 1.0
string_value = 'Hello world!'
bool_value = True
numpy_scalar0_array_value = np.array(1.0)
numpy_scalar0_bool_array_value = np.array(True)
numpy_scalar1_array_value = np.array([1.0])
numpy_scalar1_bool_array_value = np.array([True])
numpy_array_value = np.array([[1.0, 2.0],[3.0,4.0]])
list_value = [1, 1.0, 'Hello world!', True, numpy_array_value]
tuple_value = tuple(list_value)
dict_value = {'int_value':1, 'float_value':1.0, 'string_value':'Hello world!', 'bool_value':True, 'array_value':numpy_array_value}
dict_scalar_arrays_value = {'int_value':np.array([1]), 'float_value':np.array([1.0]), 'string_value':'Hello world!', 'bool_value':np.array([True]), 'array_value':numpy_array_value}
dict_deep_value = {'dict_value': dict_value, 'list_value': list_value, 'int_value': 10, 'float_value': 10.0, 'string_value': 'Goodbye world!', 'array_value': numpy_array_value*10.0}
dict_deep_scalar_arrays_value = {'dict_value': dict_scalar_arrays_value, 'list_value': list_value, 'int_value': np.array([10]), 'float_value': np.array([10.0]), 'string_value': 'Goodbye world!', 'array_value': numpy_array_value*10.0}

def dim_order_slices(arr):
    """
    This function provides two slices of a 2D array, arr: arr[:,0] and arr[0,:]. arr must be a 2D numpy array.
    """
    if not isinstance(arr, np.ndarray) or arr.ndim != 2:
        if not isinstance(arr, np.ndarray):
            in_ndim = 'N/A'
        else:
            in_ndim = arr.ndim
        
        raise TypeError('arr must be a 2D numpy array (input type was {}, input ndim was {}'.format(
                type(arr), in_ndim
            ))

    return arr[:,0], arr[0,:]

def are_collections_equal(c1, c2):
    """
    Returns a boolean indicating whether two dictionaries, possibly including numpy arrays, are equal
    (The usual method of c1 == c2 fails b/c numpy arrays do not return scalar booleans from the == 
    operator.) c1 and c2 must be dictionaries.
    """
    if not isinstance(c1, (list,tuple,dict)) or not isinstance(c2, (list,tuple,dict)):
        raise TypeError('c1 and c2 must both be collections, i.e. tuples, lists, or dicts (input type c1 = {}, input type c2 = {})'.format(
                type(c1), type(c2)
            ))

    # Comparing two collections that contain numpy arrays is tricky because the == operator
    # expects applying == to each item to return True or False, while numpy arrays return a boolean
    # array. Turns out that numpy's assert_equal test can do collections, so we use it by trapping
    # the assertion error and returning False if that happens. 
    #
    # Weird note: it seems like == works if all elements of the numpy boolean array are True, but not
    # if all are false.
    
    try:
        np.testing.assert_equal(c1,c2)
    except AssertionError:
        return False
    else:
        return True

def dict_elements_equal(d1,d2):
    for k in d1.keys():
        if k not in d2.keys():
            print('Key {} not in d2!'.format(k))
        elif not np.array_equal(d1[k], d2[k]):
            print('Values for key {} are different!'.format(k))

    for k in d2.keys():
        if k not in d1.keys():
            print('Key {} not in d1!'.format(k))
            
