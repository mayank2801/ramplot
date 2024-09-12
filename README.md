# **Ramplot**
Ramplot is a utility designed to help researchers and students visualize and analyze protein structures using Ramachandran plots. Ramachandran plots are a fundamental tool in structural biology, providing a way to visualize the torsion angles of protein backbones and assess the validity of protein structures. Ramplot offers several input options and can generate various types of plots to help understand protein conformations better.
## Introduction
Ramplot is a versatile tool for visualizing different types of steric plots, enabling the study of protein and peptide conformations in specific regions. This utility, available for both online and offline use, consists of various applications that generate a wide range of Ramachandran maps, including 2D and 3D plots for standard and specialized categories such as Glycine (Gly), Valine/Isoleucine (Val/Ile), pre-Proline (pre-Pro), trans-Proline (trans-Pro), cis-Proline (cis-Pro), and a General category representing the remaining 16 amino acids.
You can also access Ramplot as a web server at https://ramplot.in.

## Reference
 
## Standalone
The Standalone version of ramplot is written in python3 and following libraries are necessary for the successful run:
- numpy
- pandas
- biopython
- matplotlib
- mdanalysis
- pytest-shutil
- setuptools

## Installations
```
pip install ramplot

```
or Using Source Code
```
python setup.py install 
```
## USAGE
To know about the available option for the stanadlone, type the following command:
```
ramplot -h
```
To run the example, type the following command:
positional arguments:
  {pdb,trajectory,TorsionAngle,ThreeTorsionAngle}
                        Available commands
    pdb                 Run Ramachandran Plot Using PDB Files
    trajectory          Analysis of residue trajectory
    TorsionAngle        Run Ramachandran Plot Using Custom Torsion Angle CSV File
    ThreeTorsionAngle   Run Ramachandran Plot Using Custom Three Torsion Angle CSV File like PHI,PSI,THETA

To know about the available option for the pdb, type the following command:
```
ramplot pdb -h 
```
options:
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
                        Input PDB Folder Path
  -m MAPTYPE, --MapType MAPTYPE
                        Specify Map Types All :0 2d Only : 2 3D Only :3 (Default: 0)
  -r RESOLUTION, --Resolution RESOLUTION
                        Specify Resolution of Plots (Dafult resolution : 600 )
  -p PLOTFILETYPE, --PlotFileType PLOTFILETYPE
                        Specify Output Plot File formet like png jpeg tif (Default file type : png )
  -o OUTPUT, --Output OUTPUT
                        Specify Output Directory Folder Name or Path

**Input Folder:** It allow users to provide input in the folder path which contains PDB files.

**Output Folder:** Program will save the results in the Output Directory.".

**Maptype:** User should provide Map Types All : 0 2d Only : 2 3D Only : 3 (Default: 0 ).

**Resolution:** User should provide plot resolution in intiger, by default its 600. 

**Plot File Type:** User should provide plot file type, such as, png, tif and jpeg, by default its png.

Example
```
	ramplot pdb -h "Test/PDB/" -m 0 -r 600 -p png -o OutPut	
```
To know about the available option for the trajectory, type the following command:
```
ramplot trajectory -h 
```
options:
  -h, --help            show this help message and exit
  -t INPUTTPR, --InputTPR INPUTTPR
                        Input TPR File Path
  -x INPUTXTC, --InputXTC INPUTXTC
                        Input XTC File Path
  -c INPUTRESIDUES, --InputResidues INPUTRESIDUES
                        Specify Input Residie like ChainResidueNo A101
  -f FRAMEINTERVAL, --FrameInterval FRAMEINTERVAL
                        Specify Frame Interval (Default frame interval : 20
  -m MAPTYPE, --MapType MAPTYPE
                        Specify Map Types All :0 2d Only : 1 3D Only :3
  -r RESOLUTION, --Resolution RESOLUTION
                        Specify Resolution of Plots (Default resolution : 600 )
  -p PLOTFILETYPE, --PlotFileType PLOTFILETYPE
                        Specify Output Plot File formet like png jpeg tif (Dafult file type : png )
  -o OUTPUT, --Output OUTPUT
                        Specify Output Directory Folder Name or Path
						
**Input TPR File:** Provide Molecular Simulation Topology file.

**Input XTC File:** Provide Molecular Simulation Trajectory file.

**Input Residues:** Provide the chain and residue number in the format ChainResidueNumber, such as A101, A105...

**Frame Interval:** You need to provide a valid frame interval as an integer, which the program will use to extract PDB files from the trajectory. For example, if you specify an interval of 10, it will extract PDB files for every 10th frame.

**Output Folder:** Program will save the results in the Output Directory.".
Example
```
	ramplot trajectory -t Test/mdsim.tpr -x Test/mdsim.xtc -c A101,A105 -f 10 -m 0 -r 600 -p png -o OutPut	
```
To know about the available option for the torsion angle, type the following command:
```
ramplot TorsionAngle  -h 
```
options:
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
                        Input Custom Torsion Angle CSV File Path File contains 6 columns with names
                        ID,Chain,Residue,ResidueNo,PHI,PSI . For sample output file you can visit ramplot.in
  -m MAPTYPE, --MapType MAPTYPE
                        Specify Map Types All :0 2d Only : 1 3D Only :3
  -r RESOLUTION, --Resolution RESOLUTION
                        Specify Resolution of Plots (Default resolution : 600 )
  -p PLOTFILETYPE, --PlotFileType PLOTFILETYPE
                        Specify Output Plot File formet like png jpeg tif (Default file type : png )
  -o OUTPUT, --Output OUTPUT
                        Specify Output Directory Folder Name or Path
**Input Folder:**  Input Custom Torsion Angle CSV File Path File contains 6 columns with names ID,Chain,Residue,ResidueNo,PHI,PSI

**Output Folder:** Program will save the results in the Output Directory.
Example
```
	ramplot TorsionAngle -i Test/CustomTorsionAngles.csv -m 0 -r 600 -p png -o OutPut	
```
To know about the available option for the three rotatable bond angle, type the following command:
```
ramplot ThreeTorsionAngle  -h 
```
options:
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
                        Input Custom Torsion Angle CSV File Path File contains 6 columns with names
                        ID,Chain,Residue,ResidueNo,PHI,PSI . For sample output file you can visit ramplot.in
  -r RESOLUTION, --Resolution RESOLUTION
                        Specify Resolution of Plots (Default resolution : 600 )
  -p PLOTFILETYPE, --PlotFileType PLOTFILETYPE
                        Specify Output Plot File formet like png jpeg tif (Default file type : png )
  -o OUTPUT, --Output OUTPUT
                        Specify Output Directory Folder Name or Path
						
**Input Folder:**  Custom Torsion Angle CSV File Path
The file should contain 7 columns with the following headers: ID, Chain, Residue, ResidueNo, PHI, PSI, and THETA.

**Output Folder:** Program will save the results in the Output Directory.
Example
```
	ramplot ThreeTorsionAngle -i Test/CustomTorsionAnglesTheta.csv -m 0 -r 600 -p png -o OutPut	
```

Ramplot Package Files
=======================
It contantain following files, brief descript of these files given below


LICENSE                         : License information

README.md                       : This file provide information about this package

ramplot		                    : Main python program folder

test		                    : Sample data files for test

setup.py						: Python installation file


** For any query or suggestion, you can reach out to us: Mayank Kumar (mayank2801@gmail.com) and Prof. R.S. Rathore (rsrathore@cusb.ac.in). **
** Department of Bioinformatics, School of Earth Biological and Environmental Sciences, Central University of South Bihar, Gaya 824236, India ** "# ramplot" 
"# ramplot" 
