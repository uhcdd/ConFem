exe-Version with hard coded Im=False -- see ConFemAll.py -- allows for prompted computing.

A batch-script like ConFemCMRC.bat starts \verb+ComFemAll.exe+ -- which must be provided in the given path -- 
which is a wrapper for all modules. Echo lines contain examples for the prompt requests. Echo lines may 
be edited for customized requirements. Cut and paste of echo lines avoids annoying typing errors.

Working manually through prompts can be avoided with the windows powershell (PS). It allows that prompts
of programs are fed by files. PS may again be invoked by a windows batch script, see PSRun.bat. This 
feeds the content of - <name> is arbitrary - <name>.in.txt to ConFem. Each line of <nmae>.in.txt but 
the last corresponds to a prompt response. The last line may be an arbitrary character or empty line.
The echo is redirected into a log-file by \verb+PSRun.bat+.
