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
import ConPlaD
import ConSimplex
LinAlgFlag = True

class Test(unittest.TestCase):
    def setUp(self):
        """
        set up data used in the tests.
        setUp is called before each test function execution.
        """
#        self.ConFem_ = ConFem.ConFem()
        self.ConSimFem_ = ConSimFem.ConSimFem()
        self.ConPlaD_ = ConPlaD.ConPlaD()
        self.ConSimplex_ = ConSimplex.ConSimplex()
        self.NameLog = "../LogFiles"
        self.LogData = True
        """
        bond
        """
    def testBond_0(self):
        Name="../_DataBond/bond_T2"                         # CPS4, T2D3
        self.assertEqual(ConFem.ConFem().Run(Name, self.LogData,self.NameLog, False, LinAlgFlag, False,"elemout", [], False), 'f993cf99d69851b0bbf2c2268af5ffc0')
    def testBond_1(self):
        Name="../_DataBond/Pullout2"                        # CPS4, T2D2 Pullout of elastic bar in (isodamage) CP-block - eccentric bar placement along quad element edge
        self.assertEqual(ConFem.ConFem().Run(Name, self.LogData,self.NameLog, False, LinAlgFlag, False,"elemout", [], False), '8fdf324fa646d9b2af0ad249c2b3409c')
    def testBond_2(self):
        Name="../_DataC3D/Cube8"                            # C3D8, T3D3 fibers in 3D block
        self.assertEqual(ConFem.ConFem().Run(Name, self.LogData,self.NameLog, False, LinAlgFlag, False,"elemout", [], False), '4ed9a592900a0554aab50b79ca06c57d')
    def testBond_3(self):
        Name="../DataExamples/E08/E8-03B23E"               # CPS4, B23E single fiber connecting a dissected continuum
        self.assertEqual(ConFem.ConFem().Run(Name, self.LogData,self.NameLog, False, LinAlgFlag, False,"elemout", [], False), '6d0e8f63a4faaf0cdcd1927423f5e28f')
    def testBond_4(self):
        Name="../_DataBond/PulloutAxiSym"                  # CAX4, TAX2
        self.assertEqual(ConFem.ConFem().Run(Name, self.LogData,self.NameLog, False, LinAlgFlag, False,"elemout", [], False), 'b447294da04587b11961753e1a996922')

if __name__ == "__main__":
#    numpy.seterr(all='raise')
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()
