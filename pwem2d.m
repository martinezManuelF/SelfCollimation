function K0 = pwem2d(DEV,SYM,RES,PQ,MODE)
  % K0 = pwem2d(DEV,GRID,RES,PQ,MODE) Plane Wave Expansion Method
  %
  % INPUT ARGUMENTS
  % ============================================================================
  % DEV         Device structure containing the following
  %   .ER       2D Array describing permittivity of unit cell
  %   .UR       2D Array describing permeability of unit cell
  %   .T1       Reciprocal lattice vector
  %   .T2       Reciprocal lattice vector
  %
  % SYM         Key points of symmetry (for IBZ). Can be used with N points of
  %             symmetry.
  %   .NP       Number of points to be used in band diagram calculation
  %   .POINTS   Array containing all symmetry points
  %
  % RES         Grid resolution
  %   RES(1)    x-axis resolution
  %   RES(2)    y-axis resolution
  %
  % PQ          Number of harmonics along P,Q,and R
  %   .P
  %   .Q
  %
  % MODE        Mode structure containing
  %   .EM       Electromagnetic mode 'E' or 'H'
  %   .BC       1 = Band Diagram, 2 = Complete Band Diagram

  % DETERMINE GRID SIZE
  [Nx,Ny] = size(ER);

  % EXTRACT GRID PARAMETERS
  dx = RES(1);
  dy = RES(2);

  % DETERMINE TOTAL NUMBER OF HARMONICS
  M = PQ.P * PQ.Q;
