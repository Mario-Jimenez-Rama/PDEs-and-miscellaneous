  xl=0; xr=2;      % x domain
  yl=0; yr=1;      % y domain

  J= 10;           % number of points in the x-direction      
  h = (xr-xl)/J;   % mesh size in x
  Nx = J;          % number of inner points in one dimension: N-1 (Total points in one dimension: N+1)
  Ny = 0.5*Nx;     % number of inner points in y (Nx has to be even)

% Set the left side of the equation (tensor product Laplacian)
  A_1 = diag(2*ones(1,Nx-1)) + diag((-1)*ones(1,Nx-2),-1) + diag((-1)*ones(1,Nx-2),1); %grande
  A_2 = diag(2*ones(1,Ny-1)) + diag((-1)*ones(1,Ny-2),-1) + diag((-1)*ones(1,Ny-2),1); %pequeña
% Set the the Neumann boundary conditions
  C1 = zeros(Ny - 1); 
  C2 = zeros(Ny - 1);
  C3 = zeros(Nx - 1);
  C1(1,1) = 1;
  C2(Ny-1,Ny-1) = 1;
  C3(1,1) = 1;
  matA = kron(A_2, eye(Nx-1)) + kron(eye(Ny-1), A_1) - kron(C1,eye(Nx-1)) - kron(C2,eye(Nx-1)) -kron(eye(Ny-1),C3);

% Set the right side of the equation
  F_buttom = 50*h;
  F_right = 100;
  F_left = 50*h;
  F_top = -50*h;
  F1 = kron([1; zeros(Ny-2,1)],ones(Nx -1,1)*F_buttom) + kron([zeros(Ny-2,1); 1],ones(Nx -1,1)*F_top) + kron(ones(Ny -1,1), [zeros(Nx-2,1);1]*F_right) + + kron(ones(Ny -1,1), [1;zeros(Nx-2,1)]*F_left);

% Solve the linear system
  u = matA\F1; 

% Set up grids
  x = linspace(xl,xr, Nx+1);
  y = linspace(yl,yr, Ny+1);
  x_int = x(2:Nx);
  y_int = y(2:Ny);

% Reshape long 1D results onto 2D grid:
  uu = zeros(Nx+1,Ny+1); uu(2:Nx,2:Ny) = reshape(u,Nx-1,Ny-1); % rows go through x and columns through y
% Set the boundary conditions for the plot
  uu(Nx+1,1:Ny+1) = 100;
  uu(2:Nx,1) = 50*h + uu(2:Nx,2); 
  uu(2:Nx,Ny+1) = -50*h + uu(2:Nx,Ny);
  uu(1,1:Ny+1) = 50*h + uu(2,1:Ny+1); 
  [xx,yy] = meshgrid(x,y);

% Interpolate to finer grid and plot:
  [xxx,yyy] = meshgrid(0:.04:2,0:.04:1)
  uuu = interp2(xx,yy,uu',xxx,yyy,'cubic');
  figure(2), clf, mesh(xxx,yyy,uuu), colormap(hsv), colorbar, FaceColor="flat";
  xlabel x, ylabel y, zlabel u
