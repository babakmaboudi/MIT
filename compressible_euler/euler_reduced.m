function euler_reduced()
	load('reduced_basis.mat')

	gamma = 5/3;
	L = 1;
	N = 500;
	dx = L/N;
	x = linspace(0,L-dx,N)';

	[rho , rhou] = initial_cond(x);

	d = ones(1,N-1);
	Dx = diag(d,1) + diag(-d,-1); 
	Dx(1,end) = -1;
	Dx(end,1) = 1;
	Dx = Dx/dx/2;
	
	Dxx = diag(d,1) + diag(d,-1) + diag( -2*ones(N,1) );
	Dxx(1,end) = -1;
	Dxx(end,1) = 1;
	Dxx = Dx/dx/dx;

	rho_r = A'*rho;
	rhou_r = A'*rhou;

	MAX_ITER = 1000;
	dt = 0.001;

	[RHOR , RHOUR] = time_stepping(rho_r,rhou_r,Dx,Dxx,gamma,dt,MAX_ITER,dx,A);

function [RHOR , RHOUR] = time_stepping(rho_r,rhou_r,Dx,Dxx,gamma,dt,MAX_ITER,dx,A)
	
	RHOR = rho_r;
	RHOUR = rhou_r;

	for i=1:MAX_ITER

		p = 1/gamma*real((A*rho_r).^gamma);

%		[drho,drhou] = deriv(rho,rhou,p,Dx,Dxx,dt);
%
%		rho = rho + drho;
%		rhou = rhou + drhou;

		[k1r,k1ru] = deriv(rho_r,rhou_r,p,Dx,Dxx,dt,A);
		[k2r,k2ru] = deriv(rho_r+k1r/2,rhou_r+k1ru/2,p,Dx,Dxx,dt,A);
		[k3r,k3ru] = deriv(rho_r+k2r/2,rhou_r+k2ru/2,p,Dx,Dxx,dt,A);
		[k4r,k4ru] = deriv(rho_r+k3r,rhou_r+k3ru,p,Dx,Dxx,dt,A);


		rho_r = rho_r + k1r/6 + k2r/3 + k3r/3 + k4r/6;
		rhou_r = rhou_r + k1ru/6 + k2ru/3 + k3ru/3 + k4ru/6;

%			F = Filter1D(499,250,4);
%			rhout = fftshift( fft(rhou) );
%			rhout = rhout.*F;
%			rhou = real( ifft( fftshift( rhout ) ));

		[i , compute_mass(A*rho_r,dx), compute_mom(A*rhou_r,A*rho_r,dx)]

		RHOR = [ RHOR , rho_r ];
		RHOUR = [ RHOUR , rhou_r ];

		plot(A*rho_r)
		drawnow
	end

function [rho,rhou] = initial_cond(x)
	rho0 = 1;
	drho = 0.01;
	sigma = 0.08;

	e = -(x - 0.5).^2/sigma^2;
	rho = rho0 + drho*exp(e);

	rhou = zeros( length(x) , 1 );

function [drho,drhou] = deriv(rho_r,rhou_r,p,Dx,Dxx,dt,A)
	rho = A*rho_r;
	rhou = A*rhou_r;

	u = rhou./rho;
	rhouu = rhou.*u;

	drho = -0.5*Dx*rhou - 0.5*diag(rho)*Dx*u - 0.5*diag(u)*Dx*rho;
	drhou = -0.5*Dx*rhouu - 0.5*diag(rhou)*Dx*u - 0.5*diag(u)*Dx*rhou - Dx*p;% - rhou;
	
	drho = dt*A'*drho;
	drhou = dt*A'*drhou;

function mass = compute_mass(rho,dx)
	rho0 = 1;
	mass = sum( rho - rho0 )*dx;

function mom = compute_mom(rhou,rho,dx)
	
	u = rhou./rho;
	rhouu = rhou.*u;
	mom = 0.5*sum(rhouu)*dx
