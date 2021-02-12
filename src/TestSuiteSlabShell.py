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
        self.NameLog = "../_DataShellsSlabs/tmp"
       
    def testshell_01(self):
        Name="../_DataShellsSlabs/Shell"
        self.assertEqual(self.ConFem_.Run(Name,    self.NameLog, False, LinAlgFlag, False, "elemout", None, None, [], False), '82cf7598deffe469ab409e8d5aa29b38')

    def testshell_02(self):
        Name="../_DataShellsSlabs/bridge_el05m"                 #                                # folded plate
        self.assertEqual(self.ConFem_.Run(Name,    self.NameLog, False, LinAlgFlag, False, "elemout", None, None, [], False), '6069f16a9ce2983ef91d34001cc89623')

    def testshell_03(self):
        Name="../_DataShellsSlabs/c_1461(0.08)_2.1e5_0.3_segment load"                 # with SH3
        self.assertEqual(self.ConSimFem_.Run(Name, False, LinAlgFlag, "elemout"), '0dfa01a06e61e169c2d8b49dc4361b9d')

    def testshell_04(self):
        Name, nEig, niterEig ="../_DataBuckling/c_264_210000", 5, 200                 #                                # buckle
        self.assertEqual(self.ConFem_.Run(Name,    self.NameLog, False, LinAlgFlag, False, "protocol",  [ nEig, niterEig ], None, [], False), 'd3f2809672ec5f61eb2a4f81c271dc1c')

#    def testE8_02a(self):
#        Name="../DataExamples/E10/E8-01a"                # ConFem Slab as RC shell            
#        self.assertEqual(self.ConFem_.Run(Name, self.NameLog, False, LinAlgFlag, False,"elemout", None, None, []), 'fc78ec9e2e8c8aa70473ee69e5439acd')

if __name__ == "__main__":
#    numpy.seterr(all='raise')
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()
