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
        self.NameLog = "../LogFiles"
       
    def testshell_01(self):
        Name="../_DataShellsSlabs/Shell"
        self.assertEqual(self.ConFem_.Run(Name,    self.NameLog, False, LinAlgFlag, False, "elemout", None, [], False), '291afe2ecd94f1e13b2c16dcd7704350')

    def testshell_02(self):
        Name="../_DataShellsSlabs/bridge_el05m"                 #                                # folded plate
        self.assertEqual(self.ConFem_.Run(Name,    self.NameLog, False, LinAlgFlag, False, "elemout", None, [], False), '1a5da2db1089b8fb6d8829dadb41528b')

    def testshell_03(self):
        Name="../_DataShellsSlabs/c_1461(0.08)_2.1e5_0.3_segment load"                 # with SH3
        self.assertEqual(self.ConSimFem_.Run(Name, False, LinAlgFlag, "elemout"), 'bc711661fe266e64c5dde22b001605ad')

    def testshell_04(self):
        Name ="../_DataBuckling/c_264_210000"
        self.assertEqual(self.ConFem_.Run(Name,    self.NameLog, False, LinAlgFlag, False, "protocol",  None, [], False), '6537386be97166d196d0c2d2fa658398')

#    def testE8_02a(self):
#        Name="../DataExamples/E10/E8-01a"                # ConFem Slab as RC shell            
#        self.assertEqual(self.ConFem_.Run(Name, self.NameLog, False, LinAlgFlag, False,"elemout", None, None, []), 'fc78ec9e2e8c8aa70473ee69e5439acd')

if __name__ == "__main__":
#    numpy.seterr(all='raise')
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()
