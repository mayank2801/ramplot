# -*- coding: utf-8 -*-
"""
Created on Wed Aug 21 15:58:31 2024

@author: mayank
"""

# # my_package/cli.py

# import argparse
# from ramplot import TorsionAngleCalculation
# from ramplot import FunctionUserTrajectoryInputPhiPsiPlot
# def main():
#     parser = argparse.ArgumentParser(description="A CLI tool for processing data.")
#     parser.add_argument('-i', '--input', required=True, help="Input PDB Folder Path")
#     parser.add_argument('-m', '--MapType', required=True, help="Specify Map Types ")
#     parser.add_argument('-r', '--Resolution', required=True, help="Specify Resolution of Map Types ")
#     parser.add_argument('-p', '--PlotFileType', required=True, help="Specify Resolution of Map TypesPlotFileType like png jpeg")
#     parser.add_argument('-o', '--Output', required=True, help="Specify Output Directory")
#     args = parser.parse_args()    
#     result = TorsionAngleCalculation(args.input,args.Output,args.MapType,args.Resolution,args.PlotFileType)
#     print(result)
#     result =FunctionUserTrajectoryInputPhiPsiPlot(InputTPRFilePath,InputXTCFilePath,OutPutDir,InputResidues,FrameInterval,MapType,PlotResolution,PlotFileType)
    
# if __name__ == "__main__":
#     main()
    
    
    
# my_package/cli.py

import argparse
from ramplot import TorsionAngleCalculation
from ramplot import FunctionUserTrajectoryInputPhiPsiPlot
from ramplot import FunctionUserInputPhiPsiPlot
from ramplot import  FunctionUserInputPhiPsiThetaPlot

def main():
    parser = argparse.ArgumentParser(description="Ramplot")
    
    # Common arguments
    parser.add_argument('-v', '--verbose', action='store_true', help="Increase output verbosity")
    
    # Subcommands
    subparsers = parser.add_subparsers(dest='command', help='Available commands')

    # Subcommand for pdb
    parser_pdb = subparsers.add_parser('pdb', help='Run Ramachandran Plot Using PDB Files')
    parser_pdb.add_argument('-i', '--input', required=True, help="Input PDB Folder Path")
    parser_pdb.add_argument('-m', '--MapType', default=0,  help="Specify Map Types \n All :0 \n 2d Only :  2 \n 3D Only :3 (Default: 0)")
    parser_pdb.add_argument('-r', '--Resolution', default=600,  help="Specify Resolution of Plots (Dafult resolution : 600 )")
    parser_pdb.add_argument('-p', '--PlotFileType',default='png',  help="Specify Output Plot File formet like png jpeg tif (Default file type  : png )" )
    parser_pdb.add_argument('-o', '--Output', required=True, help="Specify Output Directory Folder Name or Path")
    
    # Subcommand for trajectory
    parser_trajectory = subparsers.add_parser('trajectory', help='Analysis of residue trajectory')
    parser_trajectory.add_argument('-t', '--InputTPR', required=True, help="Input TPR File Path")
    parser_trajectory.add_argument('-x', '--InputXTC', required=True, help="Input  XTC File  Path")
    parser_trajectory.add_argument('-c', '--InputResidues', required=True, help="Specify Input Residie like ChainResidueNo A101")
    parser_trajectory.add_argument('-f', '--FrameInterval',default='20',  help="Specify Frame Interval (Default frame interval : 20")
    parser_trajectory.add_argument('-m', '--MapType', default=0,  help="Specify Map Types \n All :0 \n 2d Only :  1 \n 3D Only :3")
    parser_trajectory.add_argument('-r', '--Resolution', default=600,  help="Specify Resolution of Plots (Default resolution : 600 )")
    parser_trajectory.add_argument('-p', '--PlotFileType',default='png',  help="Specify Output Plot File formet like png jpeg tif (Dafult file type  : png )" )
    parser_trajectory.add_argument('-o', '--Output', required=True, help="Specify Output Directory Folder Name or Path")
    
    #Custom Torsion angle Plot 
    parser_Custom_Torsion_Angle = subparsers.add_parser('TorsionAngle', help='Run Ramachandran Plot Using Custom Torsion Angle CSV File')
    parser_Custom_Torsion_Angle.add_argument('-i', '--input', required=True, help="Input Custom Torsion Angle CSV File Path \n File contains 6 columns with names ID,Chain,Residue,ResidueNo,PHI,PSI . For sample output file you can visit ramplot.in")
    parser_Custom_Torsion_Angle.add_argument('-m', '--MapType', default=0, help="Specify Map Types \n All :0 \n 2d Only :  1 \n 3D Only :3")
    parser_Custom_Torsion_Angle.add_argument('-r', '--Resolution', default=600,  help="Specify Resolution of Plots (Default resolution : 600 )")
    parser_Custom_Torsion_Angle.add_argument('-p', '--PlotFileType',default='png',  help="Specify Output Plot File formet like png jpeg tif (Default file type  : png )" )
    parser_Custom_Torsion_Angle.add_argument('-o', '--Output', required=True, help="Specify Output Directory Folder Name or Path")
    # Parse the arguments
    
    
    #Custom Three Torsion angle Plot 
    parser_Custom_Three_Torsion_Angle = subparsers.add_parser('ThreeTorsionAngle', help='Run Ramachandran Plot Using Custom Three Torsion Angle CSV File like PHI,PSI,THETA')
    parser_Custom_Three_Torsion_Angle.add_argument('-i', '--input', required=True, help="Input Custom Torsion Angle CSV File Path \n File contains 6 columns with names ID,Chain,Residue,ResidueNo,PHI,PSI . For sample output file you can visit ramplot.in")
    parser_Custom_Three_Torsion_Angle.add_argument('-r', '--Resolution', default=600,  help="Specify Resolution of Plots (Default resolution : 600 )")
    parser_Custom_Three_Torsion_Angle.add_argument('-p', '--PlotFileType',default='png',  help="Specify Output Plot File formet like png jpeg tif (Default file type  : png )" )
    parser_Custom_Three_Torsion_Angle.add_argument('-o', '--Output', required=True, help="Specify Output Directory Folder Name or Path")
    # Parse the arguments
    args = parser.parse_args()

    # Dispatch the appropriate function based on the subcommand
    if args.command == 'pdb':
        result = TorsionAngleCalculation(args.input,args.Output,args.MapType,args.Resolution,args.PlotFileType)
        # result = function_a(args.input, args.verbose)
    elif args.command == 'trajectory':
        result =FunctionUserTrajectoryInputPhiPsiPlot(args.InputTPR,args.InputXTC,args.Output,args.InputResidues,args.FrameInterval,args.MapType,args.Resolution,args.PlotFileType)
    elif args.command == 'TorsionAngle':
        result =FunctionUserInputPhiPsiPlot(args.input,args.Output,args.MapType,args.Resolution,args.PlotFileType)
    elif args.command == 'ThreeTorsionAngle':
        result =FunctionUserInputPhiPsiThetaPlot(args.input,args.Output,args.Resolution,args.PlotFileType)
    else:
        parser.print_help()
        return

    print(result)

if __name__ == "__main__":
    main()
