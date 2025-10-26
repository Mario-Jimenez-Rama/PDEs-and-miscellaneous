  xl=0; xr=1;    % x domain
  yl=0; yr=1;    % y domain

  J= 6;       % number of points in both x- and y-directions
  N = J;      % number of inner points in one dimension: N-1 (Total points in one dimension: N+1)
  h = (xr-xl)/J;   % mesh size

% build up the coefficient matrix
  nr = (J-1)^2;    % order of the matrix

% Set the left side of the equation (tensor product Laplacian)
  A_2 = diag(2*ones(1,N-1)) + diag((-1)*ones(1,N-2),-1) + diag((-1)*ones(1,N-2),1);
  matA = kron(A_2, eye(N-1)) + kron(eye(N-1), A_2);

% Set up grids
  x = linspace(xl,xr, N+1);
  y = linspace(yl,yr, N+1);
  x_int = x(2:N);
  y_int = y(2:N);

% Set the right side of the equation
  F_buttom = 50*x_int';
  F_right = y_int'.^2 - 101*y_int' + 50;
  F_top = -50*x_int';
  F1 = kron([1; zeros(N-2,1)],F_buttom) + kron([zeros(N-2,1); 1],F_top) + kron(F_right, [zeros(N-2,1);1]);

% Solve the linear system
  u = matA\F1; 

% Reshape long 1D results onto 2D grid:
  uu = zeros(N+1,N+1); uu(2:N,2:N) = reshape(u,N-1,N-1); % rows go through x and columns through y
  [xx,yy] = meshgrid(x,y);

% Set the Dirichlet boundary conditions for the plot
  uu(1:N+1,1) = 50*x;
  uu(1:N+1,N+1) = -50*x;
  uu(1,N+1) = 0;
  uu(N+1,1:N+1) = y.^2 -101*y +50;

% Interpolate to finer grid and plot:
  [xxx,yyy] = meshgrid(0:.04:1,0:.04:1);
  uuu = interp2(xx,yy,uu',xxx,yyy,'cubic'); 
  figure(2), clf, mesh(xxx,yyy,uuu), colormap(hsv), colorbar, FaceColor="flat";
  xlabel x, ylabel y, zlabel u
 
