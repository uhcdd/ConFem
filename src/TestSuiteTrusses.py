# TestSuite0 -- 2014-01-13
# Copyright (C) [2014] [Ulrich Haeussler-Combe]
# This program is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License (GNU GPLv3) as published by the Free Software Foundation; either version 3 of the License, or (at your option) any later version.
# This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
# You should have received a copy of the GNU General Public License along with this program; if not, see <http://www.gnu.org/licenses
#
'''
Created on 23.05.2011

@author: uhc
'''
import unittest
import ConFem
import ConSimFem
LinAlgFlag = True

class Test(unittest.TestCase):

    def setUp(self):
        """
        set up data used in the tests.
        setUp is called before each test function execution.
        """
        self.ConFem_ = ConFem.ConFem()
        self.ConSimFem_ = ConSimFem.ConSimFem()
        self.NameLog = "../_DataBeams/tmp"
       
    def test_01(self):
        Name="../_DataTrusses/staebe_1_1d"                 #                                        7s
        self.assertEqual(self.ConFem_.Run(Name, self.NameLog, False, LinAlgFlag, False, "elemout", None, None, [], False), 'b56e38907144e6a1b4c1715240f69348')

    def test_02(self):
        Name="../_DataTrusses/staebe_4"                 #                                        7s
        self.assertEqual(self.ConFem_.Run(Name, self.NameLog, False, LinAlgFlag, False, "elemout", None, None, [], False), '7b31868867dd72a706ec5cde3f89dcf0')

    def test_03(self):
        Name="../_DataTrusses/staebe3D2"                 #                                        7s
        self.assertEqual(self.ConFem_.Run(Name, self.NameLog, False, LinAlgFlag, False, "elemout", None, None, [], False), '8731b6d5b16822eb0881138aed717684')

if __name__ == "__main__":
#    numpy.seterr(all='raise')
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()
