Generally each with large 0.1 and small 0.01 m element size

Biax2D				Biaxial compression
					--> ConFem3::Biax
Tension2D			Biaxial tension ?
					? --> CaeFem Bench Spec One2D

TensionCompression2D 		Tension 1-dir with unloading, compression 3-dir 
					--> ConFem::TenCom2P
TensionShear2D			Tension/Compression 3-dir, shear 1-dir 
					--> ConFem::TenComShear
Uniax2D				tension followed by compression / compression only 
					--> ConFem3::UniaxTension
					ok IsoDam
WillamsTest2D			Willams test 
					--> ConFem::WillamsTest