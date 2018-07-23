classdef IPyMatTests < matlab.unittest.TestCase
    %IPyMatTests Unit tests for Matlab-python interface
    %   Run the unit tests with "run(IPyMatTests)". This will automatically
    %   execute the unit tests for the interface using Matlab's unit
    %   testing framework.
    
    % Currently a slight shortcoming of these tests is that since the
    % default int type in Python 2 is a 32 bit, if we want to test
    % converting Matlab ints to Python ints, we have to make them 32 bit on
    % this side, but converting Python ints to Matlab ints are
    % automatically converted to 64 bit ints, so the Py --> Mat tests have
    % to specify their test results to have 64 bit ints, whereas the
    % instance properties typically use 32 bit ints.
    
    properties
        mat_scalar;
        py_scalar;
        py_0d_array;
        py_1d_array;
        mat_row_vec;
        mat_col_vec;
        mat_2d_array;
        py_2d_array;
        mat_3d_array;
        mat_3d_array_s1;
        mat_3d_array_s2;
        mat_flat_struct;
        py_flat_dict;
        mat_flat_cell;
        py_flat_list;
        mat_deep_struct;
        py_deep_dict;
    end
    
    methods(TestClassSetup)
        function addTestDirToPythonPath(~) % instance methods need at least one input, but it's not used here.
            % Add the directory containing this file to the Python search
            % path if it isn't there.
            my_dir = fileparts(mfilename('fullpath'));
            if count(py.sys.path, my_dir) == 0
                insert(py.sys.path, int32(0), my_dir);
            end
        end
        
        function makeTestValues(testCase)
            % Python types should be defined in PySideTests. Storing
            % PySideTests values in instance variables here is a holdover
            % from the first draft, it's better to just reference
            % py.PySideTests in the individual test functions.
            
            % Scalar test variables
            testCase.mat_scalar = 1;
            testCase.py_1d_array = py.numpy.array(testCase.mat_scalar);
            testCase.py_1d_array = py.numpy.array(py.list({testCase.mat_scalar}));
            
            % Vector test variables
            testCase.mat_row_vec = zeros(1,10);
            testCase.mat_col_vec = zeros(10,1);
            
            % Array test variables
            array_2d = [1,2; 3 4];
            testCase.mat_2d_array = array_2d;
            testCase.py_2d_array = py.PySideTests.numpy_array_value;
            
            testCase.mat_3d_array = zeros(2,3,4);
            testCase.mat_3d_array_s1 = zeros(1,3,4);
            testCase.mat_3d_array_s2 = zeros(2,1,4);
            
            % Structure and list test variables
            testCase.mat_flat_struct = struct('int_value', int32(1), 'float_value', 1, 'string_value', 'Hello world!', 'bool_value', true, 'array_value', array_2d);
            testCase.py_flat_dict = py.PySideTests.dict_value;
            
            testCase.mat_flat_cell = {int32(1), 1, 'Hello world!', true, testCase.mat_2d_array};
            % py.list can't deal with a cell array with numeric arrays in
            % it; must convert the array to a numpy array first 
            py_input_cell = [testCase.mat_flat_cell(1:4), {testCase.py_2d_array}];
            testCase.py_flat_list = py.list(py_input_cell);
            
            testCase.mat_deep_struct = struct('dict_value', testCase.mat_flat_struct,...
                'list_value', {testCase.mat_flat_cell},...
                'int_value', int32(10),...
                'float_value', 10,...
                'string_value', 'Goodbye world!',...
                'array_value', testCase.mat_2d_array * 10);
            testCase.py_deep_dict = py.PySideTests.dict_deep_value;
        end
    end
    
    methods(Test)
        function testConvertScalar(testCase)
            converted_val = matlab2python(testCase.mat_scalar);
            testCase.verifyEqual(converted_val, py.PySideTests.float_value);
        end
        
        function testConvertScalarToArray(testCase)
            converted_val = matlab2python(testCase.mat_scalar, 'array1');
            testCase.verifyPyArrayEq(converted_val, py.PySideTests.numpy_scalar1_array_value);
        end
        
        function testConvertLogicalScalar(testCase)
            converted_val = matlab2python(true);
            testCase.verifyEqual(converted_val, py.PySideTests.bool_value);
        end
        
        function testConvertLogicalScalarToArray(testCase)
            converted_val = matlab2python(true, 'array1');
            testCase.verifyPyArrayEq(converted_val, py.PySideTests.numpy_scalar1_bool_array_value);
        end
        
        function testConvertArray(testCase)
            conv_py_2d_array = matlab2python(testCase.mat_2d_array);
            testCase.verifyPyArrayEq(conv_py_2d_array, testCase.py_2d_array);
        end
        
        function testArrayOrder(testCase)
            % This tests that matarray2numpyarray's 'match' and 'native'
            % ordering behave as expected, particularly, that given P =
            % matarray2numpyarray(A, dimorder), for dimorder = 'match',
            % then A(:,1) == P[:,0] and A(1,:) == P[0,:], while for
            % dimorder = 'native', A(:,1) == P[0,:], etc.
            
            matched_py_array = matarray2numpyarray(testCase.mat_2d_array, 'match');
            native_py_array = matarray2numpyarray(testCase.mat_2d_array, 'native');
            
            % We can't index numpy arrays in Matlab (at least in R2014b).
            % So we slice them in Python and return the two slices. This
            % function, given array P, will return P[:,0] and P[0,:]
            matched_slices = cell(py.PySideTests.dim_order_slices(matched_py_array));
            native_slices = cell(py.PySideTests.dim_order_slices(native_py_array));
            
            % Now we need to convert individual slices and make sure that
            % the right slices match up. Conversion from Matlab vectors to
            % numpy arrays requires that the vectors are row vectors.
            first_dim = py.numpy.array(testCase.mat_2d_array(:,1)');
            second_dim = py.numpy.array(testCase.mat_2d_array(1,:));
            
            testCase.verifyPyArrayEq(matched_slices{1}, first_dim);
            testCase.verifyPyArrayEq(matched_slices{2}, second_dim);
            testCase.verifyPyArrayEq(native_slices{1}, second_dim);
            testCase.verifyPyArrayEq(native_slices{2}, first_dim);
        end
        
        function testArrayShapeMatched(testCase)
            % This will test that numpy arrays come out with the same shape
            % (in matched mode) as the input Matlab array
            py_3d = matarray2numpyarray(testCase.mat_3d_array, 'match');
            py_3d_s1 = matarray2numpyarray(testCase.mat_3d_array_s1, 'match');
            py_3d_s2 = matarray2numpyarray(testCase.mat_3d_array_s2, 'match');
            testCase.verifyShapes(testCase.mat_3d_array, py_3d, false);
            testCase.verifyShapes(testCase.mat_3d_array_s1, py_3d_s1, false);
            testCase.verifyShapes(testCase.mat_3d_array_s2, py_3d_s2, false);
        end
        
        function testArrayShapeNative(testCase)
            % This will test that numpy arrays come out with the reverse
            % shape (in native mode) as the input Matlab array
            py_3d = matarray2numpyarray(testCase.mat_3d_array, 'native');
            py_3d_s1 = matarray2numpyarray(testCase.mat_3d_array_s1, 'native');
            py_3d_s2 = matarray2numpyarray(testCase.mat_3d_array_s2, 'native');
            testCase.verifyShapes(testCase.mat_3d_array, py_3d, true);
            testCase.verifyShapes(testCase.mat_3d_array_s1, py_3d_s1, true);
            testCase.verifyShapes(testCase.mat_3d_array_s2, py_3d_s2, true);
        end
        
        function testVecAsMatNever(testCase)
            % Test converting a row and a column vector using the 'never'
            % option of LIST_RECURSION. Both should end up as 1D numpy
            % arrays.
            py_row_vec = matarray2numpyarray(testCase.mat_row_vec, [], [], 'never');
            py_col_vec = matarray2numpyarray(testCase.mat_col_vec, [], [], 'never');
            testCase.verifyShapes(numel(testCase.mat_row_vec), py_row_vec, -1);
            testCase.verifyShapes(numel(testCase.mat_col_vec), py_col_vec, -1);
        end
        
        function testVecAsMatRow(testCase)
            % Test converting a row and a column vector using the 'never'
            % option of LIST_RECURSION. Both should end up as 1D numpy
            % arrays.
            py_row_vec = matarray2numpyarray(testCase.mat_row_vec, [], [], 'row');
            py_col_vec = matarray2numpyarray(testCase.mat_col_vec, [], [], 'row');
            testCase.verifyShapes(size(testCase.mat_row_vec), py_row_vec, -1);
            testCase.verifyShapes(numel(testCase.mat_col_vec), py_col_vec, -1);
        end
        
        function testVecAsMatColumn(testCase)
            % Test converting a row and a column vector using the 'never'
            % option of LIST_RECURSION. Both should end up as 1D numpy
            % arrays.
            py_row_vec = matarray2numpyarray(testCase.mat_row_vec, [], [], 'column');
            py_col_vec = matarray2numpyarray(testCase.mat_col_vec, [], [], 'column');
            testCase.verifyShapes(numel(testCase.mat_row_vec), py_row_vec, -1);
            testCase.verifyShapes(size(testCase.mat_col_vec), py_col_vec, -1);
        end
        
        function testVecAsMatAlways(testCase)
            % Test converting a row and a column vector using the 'never'
            % option of LIST_RECURSION. Both should end up as 1D numpy
            % arrays.
            py_row_vec = matarray2numpyarray(testCase.mat_row_vec, [], [], 'always');
            py_col_vec = matarray2numpyarray(testCase.mat_col_vec, [], [], 'always');
            testCase.verifyShapes(size(testCase.mat_row_vec), py_row_vec, -1);
            testCase.verifyShapes(size(testCase.mat_col_vec), py_col_vec, -1);
        end
        
        % TODO: Add tests for arrays in structures and cell arrays
        
        function testConvertCell(testCase)
            % Cell arrays convert to python lists.
            flat_list = matlab2python(testCase.mat_flat_cell);
            testCase.verifyPyCollEq(flat_list, py.PySideTests.list_value);
        end
        
        function testConvertStruct(testCase)
            % Structures convert to python dictionaries. Try four cases:
            % one with no recursion (no substructures), one with recursion
            % (structure with substructures), and (todo) both again with
            % non-scalar structures.
            
            flat_dict = matlab2python(testCase.mat_flat_struct);
            testCase.verifyPyCollEq(flat_dict, testCase.py_flat_dict);
            deep_dict = matlab2python(testCase.mat_deep_struct);
            testCase.verifyPyCollEq(deep_dict, testCase.py_deep_dict);
        end
        
        function testConvertStructScalarArrays(testCase)
            flat_dict = matlab2python(testCase.mat_flat_struct, 'array1');
            testCase.verifyPyCollEq(flat_dict, py.PySideTests.dict_scalar_arrays_value);
            deep_dict = matlab2python(testCase.mat_deep_struct, 'array1');
            testCase.verifyPyCollEq(deep_dict, py.PySideTests.dict_deep_scalar_arrays_value);
        end
        
        function testConvertPyInt(testCase)
            % Verify that python2matlab successfully converts a Python
            % integer to a Matlab scalar number. verifyEqual seems (in
            % R2014b) to check that the types are the same as well.
            conv_int = python2matlab(py.PySideTests.int_value);
            testCase.verifyEqual(conv_int, int64(1))
        end
        
        function testConvertPyFloat(testCase)
            % Verify that python2matlab successfully converts a Python
            % float to a Matlab scalar number. verifyEqual seems (in
            % R2014b) to check that the types are the same as well.
            conv_float = python2matlab(py.PySideTests.float_value);
            testCase.verifyEqual(conv_float, 1);
        end
        
        function testConvertPyScalarFloatArray(testCase)
            conv_float = python2matlab(py.PySideTests.numpy_scalar0_array_value);
            testCase.verifyEqual(conv_float, 1);
            conv_float = python2matlab(py.PySideTests.numpy_scalar1_array_value);
            testCase.verifyEqual(conv_float, 1);
        end
        
        function testConvertPyString(testCase)
            % Verify that python2matlab successfully converts a Python
            % string to a Matlab string. verifyEqual seems (in
            % R2014b) to check that the types are the same as well.
            conv_str = python2matlab(py.PySideTests.string_value);
            testCase.verifyEqual(conv_str, 'Hello world!');
        end
        
        function testConvertPyBool(testCase)
            % Verify that python2matlab successfully converts a Python
            % boolean to a Matlab scalar logical. verifyEqual seems (in
            % R2014b) to check that the types are the same as well.
            conv_bool = python2matlab(py.PySideTests.bool_value);
            testCase.verifyEqual(conv_bool, true);
        end
        
        function testConvertPyScalarBoolArray(testCase)
            conv_bool = python2matlab(py.PySideTests.numpy_scalar0_bool_array_value);
            testCase.verifyEqual(conv_bool, true);
            conv_bool = python2matlab(py.PySideTests.numpy_scalar1_bool_array_value);
            testCase.verifyEqual(conv_bool, true);
        end
        
        function testConvertPyArray(testCase)
            conv_array = python2matlab(py.PySideTests.numpy_array_value);
            testCase.verifyEqual(conv_array, testCase.mat_2d_array);
        end
        
        function testConvertPyList(testCase)
            % Verify that python2matlab successfully converts a Python
            % list to a Matlab cell array. verifyEqual seems (in
            % R2014b) to check that the types are the same as well.
            conv_list = python2matlab(py.PySideTests.list_value);
            test_cell = {int64(1), 1, 'Hello world!', true, testCase.mat_2d_array};
            testCase.verifyEqual(conv_list, test_cell);
        end
        
        function testConvertPyTuple(testCase)
            % Verify that python2matlab successfully converts a Python
            % tuple to a Matlab cell array. verifyEqual seems (in
            % R2014b) to check that the types are the same as well.
            conv_list = python2matlab(py.PySideTests.tuple_value);
            test_cell = {int64(1), 1, 'Hello world!', true, testCase.mat_2d_array};
            testCase.verifyEqual(conv_list, test_cell);
        end
        
        function testConvertPyDict(testCase)
            % Verify that python2matlab successfully converts a Python
            % dict to a Matlab scalar struct. verifyEqual seems (in
            % R2014b) to check that the types are the same as well.
            conv_dict = python2matlab(py.PySideTests.dict_value);
            test_struct = struct('int_value', int64(1),...
                                 'float_value', 1,...
                                 'string_value', 'Hello world!',...
                                 'bool_value', true,...
                                 'array_value', testCase.mat_2d_array);
            testCase.verifyEqual(conv_dict, test_struct);
        end
    end
    
    methods
        function verifyPyArrayEq(testCase, array1, array2)
            % Matlab's ISEQUAL() function returns true for any two Python
            % Numpy array. Since that seems to be what verifyEqual uses
            % internally, any time we expect Numpy arrays, we need to check
            % for equality this way.
            testCase.verifyTrue(py.numpy.array_equal(array1, array2));
        end
        
        function verifyPyCollEq(testCase, c1, c2)
            % Collections containing numpy arrays are a little tricky to
            % evaluate equality - it's easier to do it on the Python side
            % (see PySideTests.are_collections_equal for more detail).
            testCase.verifyTrue(py.PySideTests.are_collections_equal(c1, c2));
        end
        
        function verifyShapes(testCase, mat_array, py_array, is_native)
            % This just has the logic to convert a Numpy array's shape
            % tuple to a Matlab numeric array that we can compare to the
            % size of the input mat array
            %
            % is_native controls how the order of dimensions is compared,
            % if 0, the size of mat_array and py_array are directly
            % compared. If 1, the size of py_array is reversed first. If
            % -1, mat_array is treated as a size itself, not the array to
            % take the size of.
            if is_native >= 0
                mat_size = size(mat_array);
            else
                mat_size = mat_array;
            end
            py_size = cellfun(@double, cell(py_array.shape));
            if is_native > 0
                py_size = fliplr(py_size);
            end
            testCase.verifyEqual(mat_size, py_size);
        end
    end
    
end

