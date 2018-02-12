function WN = pwem2d(DEV,SYM,BETA,RES,PQ,MODE)
  % K0 = pwem2d(DEV,GRID,RES,PQ,MODE) Plane Wave Expansion Method
  %
  % INPUT ARGUMENTS
  % ============================================================================
  % DEV         Device structure containing the following
  %   .ER       2D Array describing permittivity of unit cell
  %   .UR       2D Array describing permeability of unit cell
  %   .T1       Reciprocal lattice vector
  %   .T2       Reciprocal lattice vector
  %   .LATTICE  Lattice constant a
  %
  % SYM         Key points of symmetry (for IBZ). Can be used with N points of
  %             symmetry.
  %   .NP       Number of points to be used in band diagram calculation
  %   .POINTS   Array containing all symmetry points
  %
  % BETA        2-by-? Array of BLOCH wave vectors
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
  %
  % OUTPUT ARGUMENTS
  % ============================================================================
  % WN          Normalized frequency

  % DETERMINE GRID SIZE
  [Nx,Ny] = size(ER);

  % EXTRACT GRID PARAMETERS
  dx = RES(1);
  dy = RES(2);

  % DETERMINE TOTAL NUMBER OF HARMONICS
  M = PQ.P * PQ.Q;

  % EXTRACT BLOCH WAVE VECTORS
  bx = BETA(1,:);
  by = BETA(2,:);
  NBX = length(bx);
  NBY = length(by);

  % COMPUTE CONVOLUTION MATRICES
  ERC = convmat(DEV.ER,PQ.P,PQ.Q);
  URC = convmat(DEV.UR,PQ.P,PQ.Q);

  % COMPUTE HARMONIC AXES
  p = [ -floor(PQ.P/2) : +floor(PQ.P/2) ];
  q = [ -floor(PQ.Q/2) : +floor(PQ.Q/2) ];

  % INITIALIZE NORMALIZED FREQUENCY ARRAYS
  if (MODE.BC == 1)
    if (MODE.EM == 'E')
      WNE = zeros(M,NBX);
    elseif (MODE.EM == 'H')
      WNH = zeros(M,NBX);
    end
  elseif (MODE.BC == 2)
  return;
  end

  % SOLVE PROBLEM
  if (MODE.BC == 1)
    for nbeta = 1 : NBX
      KX = bx(nbeta) - 2*pi*p/DEV.LATTICE;
      KY = by(nbeta) - 2*pi*q/DEV.LATTICE;
      [KY,KX] = meshgrid(KY,KX);
      KX = diag(sparse(KX(:)));
      KY = diag(sparse(KY(:)));
        if (MODE.EM == 'E')
          % Build eigen-value problem
          A = KX/URC*KX + KY/URC*KY;
          B = ERC;

          % Solve eigen-value problem
          [K0] = eig(A,B);

          % Record normalized frequencies
          K0 = sort(K0);
          K0 = (DEV.LATTICE/(2*pi))*real(sqrt(K0));
          WNE(:,nbeta) = K0;
        elseif (MODE.EM == 'H')
          % Build eigen-value problem
          A = KX/ERC*KX + KY/ERC*KY;
          B = URC;

          % Solve eigen-value problem
          [K0] = eig(A,B);

          % Record normalized frequencies
          K0 = sort(K0);
          K0 = (DEV.LATTICE/(2*pi))*real(sqrt(K0));
          WNH(:,nbeta) = K0;
        end
      end
    elseif (MODE.BC == 2)
    return;
    end
