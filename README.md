# ConeCalculator-for-VMD

The program calculates cone angles for star shaped morphologies in charged molecule systems.

Process:
1. Run the app
2. app will ask you to choose files to load, a topology file and a dcd file. type the full file name with the extension, ex: Na31-2148H2O.dcd
3. pick the frame number you are trying to use. Do not use frames 0, 1, 2 or the final frame. The program can sometimes cause issues with them
4. enter as many residue IDs along the edge of the base of the cone as you see fit (minimum 3 or results wont be accurate)
5. choose the Residue ID of the tip of the cone CORRECTLY (BE SURE TO PICK RIGHT ONE OR PROGRAM MIGHT CRASH)
6. READ STEP 5, VERY IMPORTANT
7. program repeats to step 3 until you are satisfied with the number of cones you need in your program
8. Since the data is not stored externally, if the program crashes, or you decide to end the program all data will be lost. So I suggest running a set of 10 cones at maximum so that long term progress is not lost incase of an input error or program crash

Done.


be sure to install necessary dependancies. Numpy, MDtraj, MatPlotLib

program can be ran on command line, but I prefer you run in either spyder or vscode
