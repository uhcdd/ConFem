SWIG is an interface compiler that connects programs written in C and C++ 
with scripting languages such as Perl, Python, Ruby, and Tcl. It works by
taking the declarations found in C/C++ header files and using them to 
generate the wrapper code that scripting languages need to access the 
underlying C/C++ code. 
 
currently used with swigwin-4.1.1, see also https://www.swig.org/download.html
for e.g. ConFemMatC

in working directory of ConFemMatC
	swig.bat -> ConfemMatC_wrap.cxx, ConfemMatC.py
	
Visual Studio
	create a new new project
		empty project c++ windows console
	configure new project
		project name CoonFemMatC - location ....\C_sources\ConFemMatC
	Create
	
	Drag and Drop
		source files ConFemMatC.cpp, ConfemMatC_wrap.cxx
			ConFemMatC only dsyevc3.c, dsyevj3.c, dsyevv3.c
		header files  ConFemMatC.h,  ConFemMatC.i
			ConFemMatC only dsyevc3.h, dsyevj3.h, dsyevv3.h
		
	Release & x64

	Project
		Properties
			General
				Target Name -> _ConFemMatC  - don't forget intial _
				Configuration Type -> Dynamic Library (.dll)
				take (dt. uebernehmen) - dont' forget this!
			Advanced
				Target file extension -> .pyd
				take
			C/C++
				General
					Additional Include Directiories
						-> %localappdata%\Programs\Python\Python311\Lib\site-packages\numpy\core\include
						-> %localappdata%\Programs\Python\Python311\include
					take
			Linker
				General
					Enable Incremental Linking -> NO
					Additional Library Directories
						-> %localappdata%\Programs\Python\Python311\libs
					take
	ConFemMatC Build

from working directory ConFemMatC.py, ConFemMatX\x64\Release\_ConFemMatC.pyd 
	-> ...ConFem/src/ 