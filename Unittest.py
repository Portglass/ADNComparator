import unittest
import ADN_library


class TestADN_library(unittest.TestCase):
    def test_nbMismatch(self):
        self.assertEqual(ADN_library.nbMismatch(list("ATCAATATCCACCTGCAGAT"),list("TAGTTATAGGTGGACGTCTA")), 0)
        self.assertEqual(ADN_library.nbMismatch(list("ATCAATATCCACCTGCAGAT"),list("TAAATATAGGTGGACGTCTA")), 2)


if __name__ == '__main__':
    unittest.main()
