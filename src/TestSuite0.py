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
        self.NameLog = "../DataExamples/tmp"
        """
        tensile bars
        """
    def testE2_02C1(self):
        self.ConFem1 = ConFem.ConFem()
        Name="../DataExamples/E04/E2-02C1"               # ConFem Tensile bar with creep - case 1 imposed immediate stress
        self.assertEqual(self.ConFem1.Run(Name, self.NameLog, False, LinAlgFlag, False,"elemout", None, None, [], False), '380ad4e7c0c84ddd56e0bbecec55ce9b')
    def testE2_02C2(self):
        self.ConFem2 = ConFem.ConFem()
        Name="../DataExamples/E04/E2-02C2"               # ConFem Tensile bar with creep - case 2 imposed immediate strain
        self.assertEqual(self.ConFem2.Run(Name, self.NameLog, False, LinAlgFlag, False,"elemout", None, None, [], False), '5932183a73b76a483ae3dba637e23d78')
    def testE2_02C3(self):
        self.ConFem3 = ConFem.ConFem()
        Name="../DataExamples/E04/E2-02C3"               # ConFem Tensile bar with creep - case 3 gradually imposed strain
        self.assertEqual(self.ConFem3.Run(Name, self.NameLog, False, LinAlgFlag, False,"elemout", None, None, [], False), 'e678d8c60aff5d46409356c4487bdac5')
    def testE2_04(self):
        self.ConFem4 = ConFem.ConFem()
        Name="../DataExamples/E04/E2-04"                 # ConFem Simple reinforced concrete tension bar
        self.assertEqual(self.ConFem4.Run(Name, self.NameLog, False, LinAlgFlag, False,"elemout", None, None, [], False), '233acb7630a361fdab6007e1bf042176')
        """
        beams
        """
    def testE3_02(self):
        self.ConFem5 = ConFem.ConFem()
        Name="../DataExamples/E05/E3-02"                 # ConFem Simple reinforced concrete beam
        self.assertEqual(self.ConFem5.Run(Name, self.NameLog, False, LinAlgFlag, False,"elemout", None, None, [], False), '7312a6a326144433e5d3e6e6dae9c743')
    def testE3_03(self):
        self.ConFem6 = ConFem.ConFem()
        Name="../DataExamples/E05/E3-03"                 # Confem Creep deformations of reinforced concrete beam 
        self.assertEqual(self.ConFem6.Run(Name, self.NameLog, False, LinAlgFlag, False,"elemout", None, None, [], False), '1e1940cca56828340d6b96dd1c974410')
    def testE3_04(self):
        self.ConFem7 = ConFem.ConFem()
        Name="../DataExamples/E05/E3-04"                 # ConFem Temperature actions on reinforced concrete beam
        self.assertEqual(self.ConFem7.Run(Name, self.NameLog, False, LinAlgFlag, False, "elemout", None, None, [], False), '1a9e818d1aa011587d134afd9c5ded71')
    def testE3_06(self):
        self.ConFem8 = ConFem.ConFem()
        Name="../DataExamples/E05/E3-06"                 # ConFem Prestressed reinforced concrete beam
        self.assertEqual(self.ConFem8.Run(Name, self.NameLog, False, LinAlgFlag, False,"elemout", None, None, [], False), 'c658cdca0287b74c2995076abf414928')
    def testE3_08(self):
        self.ConFem8 = ConFem.ConFem()
        Name="../DataExamples/E05/E3-08"                 # ConFem Cantilever column 2nd order
        self.assertEqual(ConFem.ConFem().Run(Name, self.NameLog, False, LinAlgFlag, False,"elemout", None, None, [], False), '4ab6ed350251fd7b3aed1b6b9b9f86fd')
    def testE3_09(self):
        self.ConFem9 = ConFem.ConFem()
        Name="../DataExamples/E05/E3-09"                 # ConFem Reinforced concrete beam under impact load
        self.assertEqual(self.ConFem9.Run(Name, self.NameLog, False, LinAlgFlag, False,"elemout", None, None, [], False), '569ff8cbb5f74d1dd1a8f22fda75910f')
    def testE3_10(self):
        self.ConFem10 = ConFem.ConFem()
        Name="../DataExamples/E05/E3-02_CircMises"                 # ConFem circ cross section with mises plasticity
        self.assertEqual(self.ConFem10.Run(Name, self.NameLog, False, LinAlgFlag, False,"elemout", None, None, [], False), '6c8e30864c73edefcb47f7ac893df0c1')
        """
        plate / strut-and-tie
        """
    def testE4_01plate(self):
        self.ConFem10 = ConFem.ConFem()
        Name="../DataExamples/E06/E4-01plate"            # ConFem Deep Beam as linear elastic plate
        self.assertEqual(self.ConFem10.Run(Name, self.NameLog, False, LinAlgFlag, False,"elemout", None, None, [], False), 'a9960234fbe52304f62a2836eb0144f2')
    def testE4_01(self):
        self.ConFem11 = ConFem.ConFem()
        Name="../DataExamples/E06/E4-01"                 # ConFem Deep Beam as strut and tie model nonlinear arc length
        self.assertEqual(self.ConFem11.Run(Name, self.NameLog, False, LinAlgFlag, False,"elemout", None, None, [], False), 'b3656161e69a6cf6d4738fe086ee427f')
    def testE4_01Simplex(self):
        Name="../DataExamples/E06/E4-01"                 # ConSim Plate strut and tie E4-01 ideal elasto plastic limit design
        self.assertEqual(self.ConSimplex_.Run(Name, LinAlgFlag, False ), 'babd4a7dc0658160e1b4dc691dabd1ff')
    def testE4_02u03(self):
        self.ConFem12 = ConFem.ConFem()
        Name="../DataExamples/E06/E4-02u03"              # ConFem Corbel as strut and tie nonlinear arc length              2s            
        self.assertEqual(self.ConFem12.Run(Name, self.NameLog, False, LinAlgFlag, False,"elemout", None, None, [], False), '01be4c0f4d437661412487820dd9eb1c')
        """
        plates
        """
    def testE6_01(self):
        Name="../DataExamples/E08/E6-02"                 # ConSimFem Plate E4_01 linear elastic
        self.assertEqual(self.ConSimFem_.Run(Name, False, LinAlgFlag, "elemout"), 'a9960234fbe52304f62a2836eb0144f2')
    def testE6_02Design(self):
        Name="../DataExamples/E08/E6-02"                 # ConPlad Plate E4_01 reinforcement design
        KeyElset, KeyTime, KeyMat = 'EL1', '1.0000', 'MAT1'
        self.assertEqual(self.ConPlaD_.Run(Name, KeyElset,KeyTime,KeyMat, False, 'plate'), '0257a2e67ba66ea679ffec3c9d820444')
        """
        slabs
        """
    def testE7_01(self):
        Name="../DataExamples/E09/E7-01"                 # SimFem Slab linear elastic
        self.assertEqual(self.ConSimFem_.Run(Name, False, LinAlgFlag,"elemout"), 'ce327dede7897ef81a5b77f3a1e20e93')
    def testE7_02u03Design(self):
        Name="../DataExamples/E09/E7-01"                 # ConPlad Slab E7_01 linear elastic, reinforcement design
        KeyElset, KeyTime, KeyMat = 'PROP1', '1.0000', 'MAT1'
        self.assertEqual(self.ConPlaD_.Run(Name, KeyElset,KeyTime,KeyMat, False, "slab"), '5512f1502194db5650689fc1e7934c85')
    def testE7_04(self):
        self.ConFem13 = ConFem.ConFem()
        Name="../DataExamples/E09/E7-04"                 # ConFem Slab nonlinear                                       
        self.assertEqual(self.ConFem13.Run(Name, self.NameLog, False, LinAlgFlag, False,"elemout", None, None, [], False), '9dff98e3f77236dbea65d97922e33800')
        """
        shells -> moved to TestSuiteSlabShell -> back
        """
    def testE8_02a(self):
        Name="../DataExamples/E10/E8-01a"                # ConFem Slab as RC shell            
        self.assertEqual(ConFem.ConFem().Run(Name, self.NameLog, False, LinAlgFlag, False,"elemout", None, None, [], False), 'fc78ec9e2e8c8aa70473ee69e5439acd')
        """
        bond
        """
    def testBond_0(self):
        Name="../_DataDiv/bond_T2"                #            
        self.assertEqual(ConFem.ConFem().Run(Name, self.NameLog, False, LinAlgFlag, False,"elemout", None, None, [], False), '08d35cf7391e72910d1f571a03e4e396')
        """
        specimen
        """
    def testSpec_0(self):
        Name="../_DataSpecimen/One3D/WillamsTest"           # microplane with crack band              #            
        self.assertEqual(ConFem.ConFem().Run(Name, self.NameLog, False, LinAlgFlag, False,"elemout", None, None, [], False), '8a13f8de5850d23436a106e65bb41ef1')

if __name__ == "__main__":
#    numpy.seterr(all='raise')
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()
