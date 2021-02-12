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
       
    def testE3_02_B21(self):
        Name="../_DataBeams/E3-02_B21"
        self.assertEqual(self.ConSimFem_.Run(Name, False, LinAlgFlag,"elemout"), 'eb592ab8387d41e656a6df55643430c1')

    def testE3_02_B21E(self):
        Name="../_DataBeams/E3-02_B21E"                 #
        self.assertEqual(self.ConFem_.Run(Name, self.NameLog, False, LinAlgFlag, False, "elemout", None, None, [], False), '421630929819f81404fdc062c157cb53')

    def testE3_02_B23(self):
        Name="../_DataBeams/E3-02_B23"                 #
        self.assertEqual(self.ConSimFem_.Run(Name, False, LinAlgFlag,"elemout"), 'f7551ce64521dd8abd316253b094484f')

    def testE3_02_B21E_PolyL(self):
        Name="../_DataBeams/E3-02_B23E_PolyL"                 #
        self.assertEqual(self.ConFem_.Run(Name, self.NameLog, False, LinAlgFlag, False, "elemout", None, None, [], False), 'd9a985117e6d9e8d32071a7d056c3a0b')

    def testE3_02_Tens(self):
        Name="../_DataBeams/E3-02_Tens"                 #
        self.assertEqual(self.ConFem_.Run(Name, self.NameLog, False, LinAlgFlag, False, "elemout", None, None, [], False), '3930c470d6f3d0eb47109002ea884fe6')

if __name__ == "__main__":
#    numpy.seterr(all='raise')
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()
