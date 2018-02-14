# SelfCollimation
A series of MATLAB scripts that simulate periodic structures using the Plane Wave Expansion Method (PWEM) to obtain electromagnetic band diagrams and isofrequency contours. The scripts also allow the user to design a self collimating periodic structure given the design parameters.

## Included Files:
- convmat.m	: Calculates convolution matrices for 1D, 2D, and 3D problems.
- pwem2d.m 	: Generic PWEM calculator function for periodic supercells.
- BandDiagrams.m: Generates band diagrams, full band diagrams, and isofrequency contours for generic unit cells.

## Sample Output:
|**Description**      | **Image**                                                                                                 |
|:-------------------:|:---------------------------------------------------------------------------------------------------------:|
|Unit Cell            |![Unit Cell](https://github.com/martinezManuelF/SelfCollimation/blob/master/Graphics/UnitCellMarked.png)   |
|E Mode Band Diagram  |![EMBD](https://github.com/martinezManuelF/SelfCollimation/blob/master/Graphics/EModeBandDiagram.png)      |
|E Mode Full Band     |![EMFBD](https://github.com/martinezManuelF/SelfCollimation/blob/master/Graphics/EModeFullBandDiagram.png) |
|Isofrequency 1st     |![1stB](https://github.com/martinezManuelF/SelfCollimation/blob/master/Graphics/EMode1stBandISO.png)       |
|Isofrequency 2nd     |![2ndB](https://github.com/martinezManuelF/SelfCollimation/blob/master/Graphics/EMode2ndBandISO.png)       |
|H Mode Band Diagram  |![HMBD](https://github.com/martinezManuelF/SelfCollimation/blob/master/Graphics/HModeBandDiagram.png)      |
|H Mode Full Band     |![HMFBD](https://github.com/martinezManuelF/SelfCollimation/blob/master/Graphics/HModeFullBandDiagram.png) |
|Isofrequency 1st     |![1stB](https://github.com/martinezManuelF/SelfCollimation/blob/master/Graphics/HMode1stBandISO.png)       |
|Isofrequency 2nd     |![2ndB](https://github.com/martinezManuelF/SelfCollimation/blob/master/Graphics/HMode2ndBandISO.png)       |

Created as part of the course work for EE 5322--21st Century Electromagnetics at the University of Texas at El Paso.
