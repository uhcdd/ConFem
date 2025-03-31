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
        tensile bars
        """
    def testE3_02C1(self):
        self.ConFem1 = ConFem.ConFem()
        Name="../DataExamples/E03/E3-02C1"               # ConFem Tensile bar with creep - case 1 imposed immediate stress
        self.assertEqual(self.ConFem1.Run(Name, self.LogData,self.NameLog, False, LinAlgFlag, False,"elemout", [], False), '7406076a8460cceff4eed5d221ef23fb')
    def testE3_02C2(self):
        self.ConFem2 = ConFem.ConFem()
        Name="../DataExamples/E03/E3-02C2"               # ConFem Tensile bar with creep - case 2 imposed immediate strain
        self.assertEqual(self.ConFem2.Run(Name, self.LogData,self.NameLog, False, LinAlgFlag, False,"elemout", [], False), 'fd497c8139f44f9f629baeb5e102fa07')
    def testE3_02C3(self):
        self.ConFem3 = ConFem.ConFem()
        Name="../DataExamples/E03/E3-02C3"               # ConFem Tensile bar with creep - case 3 gradually imposed strain
        self.assertEqual(self.ConFem3.Run(Name, self.LogData,self.NameLog, False, LinAlgFlag, False,"elemout", [], False), '65d0888f9293181aa261d6e59ab902db')
    def testE3_04(self):
        self.ConFem4 = ConFem.ConFem()
        Name="../DataExamples/E03/E3-04"                 # ConFem Simple reinforced concrete tension bar
        self.assertEqual(self.ConFem4.Run(Name, self.LogData,self.NameLog, False, LinAlgFlag, False,"elemout", [], False), '9c9d2bef8e2d84594de8b259a4c60386')
        """
        beams
        """
    def testE4_02(self):
        self.ConFem5 = ConFem.ConFem()
        Name="../DataExamples/E04/E3-02"                 # ConFem Simple reinforced concrete beam - old version
        self.assertEqual(self.ConFem5.Run(Name, self.LogData,self.NameLog, False, LinAlgFlag, False,"elemout", [], False), 'dcb94757c72f339f60880059eb4cb0a7')
    def testE4_03(self):
        self.ConFem6 = ConFem.ConFem()
        Name="../DataExamples/E04/E4-03"                 # Confem Creep deformations of reinforced concrete beam 
        self.assertEqual(self.ConFem6.Run(Name, self.LogData,self.NameLog, False, LinAlgFlag, False,"elemout", [], False), 'ec36505c704ec68f6f8ec7844cbecdec')
    def testE4_04(self):
        self.ConFem7 = ConFem.ConFem()
        Name="../DataExamples/E04/E4-04"                 # ConFem Temperature actions on reinforced concrete beam
        self.assertEqual(self.ConFem7.Run(Name, self.LogData,self.NameLog, False, LinAlgFlag, False, "elemout", [], False), 'a2ec469fba8ebce7d2bfaef37179637c')
    def testE4_06(self):
        self.ConFem8 = ConFem.ConFem()
        Name="../DataExamples/E04/E4-06"                 # ConFem Prestressed reinforced concrete beam
        self.assertEqual(self.ConFem8.Run(Name, self.LogData,self.NameLog, False, LinAlgFlag, False,"elemout", [], False), '8cad1a7a1d59237c9c35121fe26f3e9b')
    def testE4_08(self):
        self.ConFem8 = ConFem.ConFem()
        Name="../DataExamples/E04/E4-08"                 # ConFem Cantilever column 2nd order
        self.assertEqual(ConFem.ConFem().Run(Name, self.LogData,self.NameLog, False, LinAlgFlag, False,"elemout", [], False), '616f4fe40a5a073cc0cb7e52c836a997')
    def testE4_09(self):
        self.ConFem9 = ConFem.ConFem()
        Name="../DataExamples/E04/E4-09"                 # ConFem Reinforced concrete beam under impact load
        self.assertEqual(self.ConFem9.Run(Name, self.LogData,self.NameLog, False, LinAlgFlag, False,"elemout", [], False), '9a4b21b900a7ca5e67a37fb77aab2d28')
    def testE4_10(self):
        self.ConFem10 = ConFem.ConFem()
        Name="../DataExamples/E04/E3-02_CircMises"                 # ConFem circ cross section with mises plasticity
        self.assertEqual(self.ConFem10.Run(Name, self.LogData,self.NameLog, False, LinAlgFlag, False,"elemout", [], False), 'f623691ea4a79b91e4f2cde273f8cbd5')
        """
        plate / strut-and-tie
        """
    def testE5_01plate(self):
        self.ConFem10 = ConFem.ConFem()
        Name="../DataExamples/E05/E5-01plate"            # ConFem Deep Beam as linear elastic plate
        self.assertEqual(self.ConFem10.Run(Name, self.LogData,self.NameLog, False, LinAlgFlag, False,"elemout", [], False), 'e547af43b68e5f580b84f776406dfee8')
    def testE5_01(self):
        Name="../DataExamples/E05/E5-02"                 # ConFem Deep Beam as strut and tie model nonlinear arc length
        self.assertEqual(ConFem.ConFem().Run(Name, self.LogData,self.NameLog, False, LinAlgFlag, False, "elemout", [], False),'1b5db2cb6cce2d13450b4f645621ac70')
    def testE5_01Simplex(self):
        Name="../DataExamples/E05/E5-02"                 # ConSim Plate strut and tie E4-01 ideal elasto plastic limit design
        self.assertEqual(self.ConSimplex_.Run(Name, LinAlgFlag, False ), '666cb390b4db073e0ac044fb28e77ab7')
    def testE5_02u03(self):
        self.ConFem12 = ConFem.ConFem()
        Name="../DataExamples/E05/E5-03"              # ConFem Corbel as strut and tie nonlinear arc length              2s            
        self.assertEqual(self.ConFem12.Run(Name, self.LogData,self.NameLog, False, LinAlgFlag, False,"elemout", [], False), 'eca9e193a5da7b969b7b2674b66433b6')
        """
        plates
        """
    def testE8_00(self):
        Name="../DataExamples/E08/E8-01"                 # ConSimFem Plate E4_01 linear elastic
        self.assertEqual(self.ConSimFem_.Run(Name, False, LinAlgFlag, "elemout"), 'e547af43b68e5f580b84f776406dfee8')
    def testE8_02Design(self):
        Name="../DataExamples/E08/E8-01"                 # ConPlad Plate E4_01 reinforcement design
        KeyElset = 'EL1'
        self.assertEqual(self.ConPlaD_.Run(Name, KeyElset, False, 'plate'), '0257a2e67ba66ea679ffec3c9d820444')
        """
        slabs
        """
    def testE9_01(self):
        Name="../DataExamples/E09/E9-01"                 # SimFem Slab linear elastic
        self.assertEqual(self.ConSimFem_.Run(Name, False, LinAlgFlag,"elemout"), '43b8146b5a7d2a1843d8b3cabece9281')
    def testE9_02u03Design(self):
        Name="../DataExamples/E09/E9-01"                 # ConPlad Slab E7_01 linear elastic, reinforcement design
        KeyElset = 'PROP1'
        self.assertEqual(self.ConPlaD_.Run(Name, KeyElset, False, "slab"), '56d08e553b7712a5bcf3a02374736263')
    def testE9_04(self):
        self.ConFem13 = ConFem.ConFem()
        Name="../DataExamples/E09/E9-04"                 # ConFem Slab nonlinear                                       
        self.assertEqual(self.ConFem13.Run(Name, self.LogData,self.NameLog, False, LinAlgFlag, False,"elemout", [], False), '08b99b30f2ba7561a2e283de3a95f91e')
        """
        shells
        """
    def testE10_02a(self):
        Name="../DataExamples/E10/E10-02a"                # ConFem Slab as RC shell            
        self.assertEqual(ConFem.ConFem().Run(Name, self.LogData,self.NameLog, False, LinAlgFlag, False,"elemout", [], False), '1f0abe63ff2bd3e24c69638102e35ff5')
        """
        bond
        """
    def testBond_0(self):
        Name="../_DataBond/bond_T2"                #
        self.assertEqual(ConFem.ConFem().Run(Name, self.LogData,self.NameLog, False, LinAlgFlag, False,"elemout", [], False), '1a9768df03c1e841e488b0f2abe62324')
        """
        specimen
        """
    def testSpec_0(self):
        Name="../_DataSpecimen/One3D/WillamsTest"           # microplane with crack band              #            
        self.assertEqual(ConFem.ConFem().Run(Name, self.LogData,self.NameLog, False, LinAlgFlag, False,"elemout", [], False), '04e826f770b3c9b90f5b0ed6f4f13644')

if __name__ == "__main__":
#    numpy.seterr(all='raise')
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()
