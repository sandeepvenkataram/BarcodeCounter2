import unittest

from lib import argparse


class ArgparseTests(unittest.TestCase):
    
    def test_basics(self):
        expected_keys = [arg['dest'] for arg in argparse.AVAILABLE_ARGUMENTS]
        res_dict = argparse.generate_argparser('dummy description')
        for k, v in res_dict.items():
            self.assertTrue(k in expected_keys)
            expected_val = [arg['default'] for arg in argparse.AVAILABLE_ARGUMENTS if arg['dest']==k and 'default' in arg]
            self.assertLessEqual(len(expected_val), 1)
            if len(expected_val) == 1:
                self.assertEqual(expected_val[0], v)
                

if __name__ == '__main__':
    unittest.main()