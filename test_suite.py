import unittest

from test.argparse_tests import ArgparseTests
from test.io_tests import IOTests
from test.demultiplex_tests import DemultiplexTests

def generate_test_suite() -> unittest.TestSuite:
    test_classes_to_run = [ArgparseTests, IOTests, DemultiplexTests]
    
    loader = unittest.TestLoader()
    
    suites_list = []
    for test_class in test_classes_to_run:
        suite = loader.loadTestsFromTestCase(test_class)
        suites_list.append(suite)
        
    return unittest.TestSuite(suites_list)

if __name__ == '__main__':
    runner = unittest.TextTestRunner()
    runner.run(generate_test_suite())