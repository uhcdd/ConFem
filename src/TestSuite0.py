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
        """
        tensile bars
        """
    def testE3_02C1(self):
        self.ConFem1 = ConFem.ConFem()
        Name="../DataExamples/E03/E3-02C1"               # ConFem Tensile bar with creep - case 1 imposed immediate stress
        self.assertEqual(self.ConFem1.Run(Name, self.NameLog, False, LinAlgFlag, False,"elemout", None, None, [], False), 'c5ed5fc0064ac2e283dd2f04664dc7d9')
    def testE3_02C2(self):
        self.ConFem2 = ConFem.ConFem()
        Name="../DataExamples/E03/E3-02C2"               # ConFem Tensile bar with creep - case 2 imposed immediate strain
        self.assertEqual(self.ConFem2.Run(Name, self.NameLog, False, LinAlgFlag, False,"elemout", None, None, [], False), '693ba44f767166abaa5b16de2cf4de0e')
    def testE3_02C3(self):
        self.ConFem3 = ConFem.ConFem()
        Name="../DataExamples/E03/E3-02C3"               # ConFem Tensile bar with creep - case 3 gradually imposed strain
        self.assertEqual(self.ConFem3.Run(Name, self.NameLog, False, LinAlgFlag, False,"elemout", None, None, [], False), 'ed9402914f00e0a396bc0e2484f1af0a')
    def testE3_04(self):
        self.ConFem4 = ConFem.ConFem()
        Name="../DataExamples/E03/E3-04"                 # ConFem Simple reinforced concrete tension bar
        self.assertEqual(self.ConFem4.Run(Name, self.NameLog, False, LinAlgFlag, False,"elemout", None, None, [], False), '81f1bdd92c7203431800ba2e23aedbfa')
        """
        beams
        """
    def testE4_02(self):
        self.ConFem5 = ConFem.ConFem()
        Name="../DataExamples/E04/E3-02"                 # ConFem Simple reinforced concrete beam - old version
        self.assertEqual(self.ConFem5.Run(Name, self.NameLog, False, LinAlgFlag, False,"elemout", None, None, [], False), '2554e3af22311af866090970a535fa6a')
    def testE4_03(self):
        self.ConFem6 = ConFem.ConFem()
        Name="../DataExamples/E04/E4-03"                 # Confem Creep deformations of reinforced concrete beam 
        self.assertEqual(self.ConFem6.Run(Name, self.NameLog, False, LinAlgFlag, False,"elemout", None, None, [], False), '9e4a27be6fc2ec4ec296e9b3e6f0ccb2')
    def testE4_04(self):
        self.ConFem7 = ConFem.ConFem()
        Name="../DataExamples/E04/E4-04"                 # ConFem Temperature actions on reinforced concrete beam
        self.assertEqual(self.ConFem7.Run(Name, self.NameLog, False, LinAlgFlag, False, "elemout", None, None, [], False), '14919558a817061a1a7967085f817067')
    def testE4_06(self):
        self.ConFem8 = ConFem.ConFem()
        Name="../DataExamples/E04/E4-06"                 # ConFem Prestressed reinforced concrete beam
        self.assertEqual(self.ConFem8.Run(Name, self.NameLog, False, LinAlgFlag, False,"elemout", None, None, [], False), '698134caac0963a0a821880be492d06c')
    def testE4_08(self):
        self.ConFem8 = ConFem.ConFem()
        Name="../DataExamples/E04/E4-08"                 # ConFem Cantilever column 2nd order
        self.assertEqual(ConFem.ConFem().Run(Name, self.NameLog, False, LinAlgFlag, False,"elemout", None, None, [], False), 'e487212a2819cf9718ef5478520d7525')
    def testE4_09(self):
        self.ConFem9 = ConFem.ConFem()
        Name="../DataExamples/E04/E4-09"                 # ConFem Reinforced concrete beam under impact load
        self.assertEqual(self.ConFem9.Run(Name, self.NameLog, False, LinAlgFlag, False,"elemout", None, None, [], False), '58ac4c4c2b5e47da07d7d3172f721c6f')
    def testE4_10(self):
        self.ConFem10 = ConFem.ConFem()
        Name="../DataExamples/E04/E3-02_CircMises"                 # ConFem circ cross section with mises plasticity
        self.assertEqual(self.ConFem10.Run(Name, self.NameLog, False, LinAlgFlag, False,"elemout", None, None, [], False), '286166846a1f62789d2c1d86dd870b38')
        """
        plate / strut-and-tie
        """
    def testE5_01plate(self):
        self.ConFem10 = ConFem.ConFem()
        Name="../DataExamples/E05/E5-01plate"            # ConFem Deep Beam as linear elastic plate
        self.assertEqual(self.ConFem10.Run(Name, self.NameLog, False, LinAlgFlag, False,"elemout", None, None, [], False), 'ff12989226a1a873f0eeda2c38a6d9b7')
    def testE5_01(self):
        Name="../DataExamples/E05/E5-02"                 # ConFem Deep Beam as strut and tie model nonlinear arc length
        self.assertEqual(ConFem.ConFem().Run(Name, self.NameLog, False, LinAlgFlag, False, "elemout", None, None, [], False),'f084961ab296a5696882dfb4bb0e3e44')
    def testE5_01Simplex(self):
        Name="../DataExamples/E05/E5-02"                 # ConSim Plate strut and tie E4-01 ideal elasto plastic limit design
        self.assertEqual(self.ConSimplex_.Run(Name, LinAlgFlag, False ), '666cb390b4db073e0ac044fb28e77ab7')
    def testE5_02u03(self):
        self.ConFem12 = ConFem.ConFem()
        Name="../DataExamples/E05/E5-03"              # ConFem Corbel as strut and tie nonlinear arc length              2s            
        self.assertEqual(self.ConFem12.Run(Name, self.NameLog, False, LinAlgFlag, False,"elemout", None, None, [], False), '7cefcf64d60beb9741d6e7203bc2066c')
        """
        plates
        """
    def testE8_00(self):
        Name="../DataExamples/E08/E8-01"                 # ConSimFem Plate E4_01 linear elastic
        self.assertEqual(self.ConSimFem_.Run(Name, False, LinAlgFlag, "elemout"), 'ff12989226a1a873f0eeda2c38a6d9b7')
    def testE8_02Design(self):
        Name="../DataExamples/E08/E8-01"                 # ConPlad Plate E4_01 reinforcement design
        KeyElset = 'EL1'
        self.assertEqual(self.ConPlaD_.Run(Name, KeyElset, False, 'plate'), '0257a2e67ba66ea679ffec3c9d820444')
        """
        slabs
        """
    def testE9_01(self):
        Name="../DataExamples/E09/E9-01"                 # SimFem Slab linear elastic
        self.assertEqual(self.ConSimFem_.Run(Name, False, LinAlgFlag,"elemout"), '9d4850e2a7a029e0c7ce3458e706296d')
    def testE9_02u03Design(self):
        Name="../DataExamples/E09/E9-01"                 # ConPlad Slab E7_01 linear elastic, reinforcement design
        KeyElset = 'PROP1'
        self.assertEqual(self.ConPlaD_.Run(Name, KeyElset, False, "slab"), '56d08e553b7712a5bcf3a02374736263')
    def testE9_04(self):
        self.ConFem13 = ConFem.ConFem()
        Name="../DataExamples/E09/E9-04"                 # ConFem Slab nonlinear                                       
        self.assertEqual(self.ConFem13.Run(Name, self.NameLog, False, LinAlgFlag, False,"elemout", None, None, [], False), '6fc9671a48335d85c1cc21114ee56132')
        """
        shells
        """
    def testE10_02a(self):
        Name="../DataExamples/E10/E10-02a"                # ConFem Slab as RC shell            
        self.assertEqual(ConFem.ConFem().Run(Name, self.NameLog, False, LinAlgFlag, False,"elemout", None, None, [], False), 'c38dc136dd7e2ce9509a635e339f3d44')
        """
        bond
        """
    def testBond_0(self):
        Name="../_DataBond/bond_T2"                #
        self.assertEqual(ConFem.ConFem().Run(Name, self.NameLog, False, LinAlgFlag, False,"elemout", None, None, [], False), 'f993cf99d69851b0bbf2c2268af5ffc0')
        """
        specimen
        """
    def testSpec_0(self):
        Name="../_DataSpecimen/One3D/WillamsTest"           # microplane with crack band              #            
        self.assertEqual(ConFem.ConFem().Run(Name, self.NameLog, False, LinAlgFlag, False,"elemout", None, None, [], False), 'eb43d6286e371efa17b92c3f84e9623c')

if __name__ == "__main__":
#    numpy.seterr(all='raise')
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()
