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
        self.LogData = True

    def testE3_02_B21(self):
        Name="../_DataBeams/E3-02_B21"
        self.assertEqual(self.ConSimFem_.Run(Name, False, LinAlgFlag,"elemout"), '5bf5e17936110e02dbd2594ee03a8c87')

    def testE3_02_B21E(self):
        Name="../_DataBeams/E3-02_B21E"                 #
        self.assertEqual(self.ConFem_.Run(Name, self.LogData,self.NameLog, False, LinAlgFlag, False, "elemout", [], False), '5c143cebb3834d0784c9dea533a5689c')

    def testE3_02_B23(self):
        Name="../_DataBeams/E3-02_B23"                 #
        self.assertEqual(self.ConSimFem_.Run(Name, False, LinAlgFlag,"elemout"), '23ba45c110feadfd0e931339c74fda42')

    def testE3_02_B23E_PolyL(self):
        Name="../_DataBeams/E3-02_B23E_PolyL"                 #
        self.assertEqual(self.ConFem_.Run(Name, self.LogData,self.NameLog, False, LinAlgFlag, False, "elemout", [], False), '4540d85d584c1cef118f11aef688115c')

    def testE3_02_Tens(self):
        Name="../DataExamples/E04/E4-02"                 # introducing example for RC beams whereby regarding concrete tensile strength
        self.assertEqual(self.ConFem_.Run(Name, self.LogData,self.NameLog, False, LinAlgFlag, False, "elemout", [], False), '3bc2c7588cb6dd28bd7cbc6b0e613685')

if __name__ == "__main__":
#    numpy.seterr(all='raise')
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()
