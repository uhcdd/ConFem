exe-Version with hard coded Im=True -- see ConFemAll.py -- allows for batch computing.

A batch-script like ConFEM1.bat starts a powershell-script ConFem.ps1 with a file of 
dataset names -- without suffix like .in.txt -- as input.

The powershell script ConFEM.ps1 reads line by line of the input file and for every
line starts ConFemAll.exe -- which must be provided in the given path. This is done
sequentially, i.e. one dataset after the other is processed by ConFEM.exe.

A number of batch-scrips with different names and different files of dataset names
and different data set names may be startet in parallel for a better utilization 
of processor kernels. 
The four given here serve only as example, the number is arbitrary.