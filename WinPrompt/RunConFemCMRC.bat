@echo off
echo 
echo 0 .\DataExamplesCMRC\2.1\E3-04    reinforced concrete tension bar
echo 0 .\DataExamplesCMRC\6.1\E4-01, prop1, mat1,  0, -1, -2   Reinforced concrete cross section
echo 0 .\DataExamplesCMRC\6.2\E4-02    Simple reinforced concrete beam
echo 0 .\DataExamplesCMRC\6.3\E4-04    effect of temperature actions on RCbeam
echo 0 .\DataExamplesCMRC\6.4\E4-06    prestressed RC beam
echo 0 .\DataExamplesCMRC\6.5\E4-09    beam under impact load
echo 0 .\DataExamplesCMRC\7.1\E3-01    tension bar with localization
echo 0 .\DataExamplesCMRC\7.3\UniaxTenBar_6 tensile bar with strain softening 
echo 0 .\DataExamplesCMRC\7.4\E7-01a   notched plate 
echo 0 .\DataExamplesCMRC\7.5\UniaxTenBar_6 tensile bar with strain softening and regularization 
echo 0 .\DataExamplesCMRC\7.6\E7-01a   notched plate with regularization
echo 0 .\DataExamplesCMRC\7.7\UniaxTenBar_12R Numerical reinforced tensile bar
echo 0 .\DataExamplesCMRC\8.1\Uniax2D  one element uniaxial tension / compression
echo 0 .\DataExamplesCMRC\9.1\Biax2D   one element biaxial compression
echo 0 .\DataExamplesCMRC\9.2\LSP_1    L-Shaped Panel
echo 0 .\DataExamplesCMRC\9.3\E8-04    Reinforced Plate
echo ---
.\ConFemAll\ConFemAll.exe
pause