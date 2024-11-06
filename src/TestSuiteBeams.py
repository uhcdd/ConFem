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
       
    def testE3_02_B21(self):
        Name="../_DataBeams/E3-02_B21"
        self.assertEqual(self.ConSimFem_.Run(Name, False, LinAlgFlag,"elemout"), 'd28ea08ad56a44300e2feea6a936d3a4')

    def testE3_02_B21E(self):
        Name="../_DataBeams/E3-02_B21E"                 #
        self.assertEqual(self.ConFem_.Run(Name, self.NameLog, False, LinAlgFlag, False, "elemout", None, None, [], False), '72b9694c4e8e981af105beb5033a6a10')

    def testE3_02_B23(self):
        Name="../_DataBeams/E3-02_B23"                 #
        self.assertEqual(self.ConSimFem_.Run(Name, False, LinAlgFlag,"elemout"), '1ee342ee800fb7331f28bbde4e7b7379')

    def testE3_02_B23E_PolyL(self):
        Name="../_DataBeams/E3-02_B23E_PolyL"                 #
        self.assertEqual(self.ConFem_.Run(Name, self.NameLog, False, LinAlgFlag, False, "elemout", None, None, [], False), 'dc08427e96d81e453d65e7087c618a8d')

    def testE3_02_Tens(self):
        Name="../DataExamples/E04/E4-02"                 # introducing example for RC beams whereby regarding concrete tensile strength
        self.assertEqual(self.ConFem_.Run(Name, self.NameLog, False, LinAlgFlag, False, "elemout", None, None, [], False), '995928fe76bbe4b20043d0b9b33d8586')

if __name__ == "__main__":
#    numpy.seterr(all='raise')
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()
