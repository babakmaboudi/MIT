function euler()
	
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

	MAX_ITER = 1000;
	dt = 0.001;
	[RHO , RHOU] = time_stepping(rho,rhou,Dx,Dxx,gamma,dt,MAX_ITER,dx);

	save('snaps.mat','RHO','RHOU')

function [RHO , RHOU] = time_stepping(rho,rhou,Dx,Dxx,gamma,dt,MAX_ITER,dx)
	
	RHO = rho;
	RHOU = rhou;

	for i=1:MAX_ITER
		i

		p = 1/gamma*real(rho.^gamma);

%		[drho,drhou] = deriv(rho,rhou,p,Dx,Dxx,dt);
%
%		rho = rho + drho;
%		rhou = rhou + drhou;

		[k1r,k1ru] = deriv(rho,rhou,p,Dx,Dxx,dt);
		[k2r,k2ru] = deriv(rho+k1r/2,rhou+k1ru/2,p,Dx,Dxx,dt);
		[k3r,k3ru] = deriv(rho+k2r/2,rhou+k2ru/2,p,Dx,Dxx,dt);
		[k4r,k4ru] = deriv(rho+k3r,rhou+k3ru,p,Dx,Dxx,dt);


		rho = rho + k1r/6 + k2r/3 + k3r/3 + k4r/6;
		rhou = rhou + k1ru/6 + k2ru/3 + k3ru/3 + k4ru/6;

%		if(mod(i,10) == 1)
			F = Filter1D(499,250,4);
%			plot(F)
%			pause

%			rhot = fftshift( fft(rho) );
%			rhot = rhot.*F;
%			rho = real( ifft( fftshift( rhot ) ));

			rhout = fftshift( fft(rhou) );
			rhout = rhout.*F;
			rhou = real( ifft( fftshift( rhout ) ));
%		end

%		[i , compute_mass(rho,dx), compute_mom(rhou,rho,dx)]

		RHO = [ RHO , rho ];
		RHOU = [ RHOU , rhou ];

		plot(rho)
		drawnow
	end

function [drho,drhou] = deriv(rho,rhou,p,Dx,Dxx,dt)
	u = rhou./rho;
	rhouu = rhou.*u;

	drho = -0.5*Dx*rhou - 0.5*diag(rho)*Dx*u - 0.5*diag(u)*Dx*rho;
	drhou = -0.5*Dx*rhouu - 0.5*diag(rhou)*Dx*u - 0.5*diag(u)*Dx*rhou - Dx*p;% - rhou;
	
	drho = dt*drho;
	drhou = dt*drhou;


function [rho,rhou] = initial_cond(x)
	rho0 = 1;
	drho = 0.01;
	sigma = 0.08;

	e = -(x - 0.5).^2/sigma^2;
	rho = rho0 + drho*exp(e);

	rhou = zeros( length(x) , 1 );

function filterdiag = Filter1D(N,Nc,s)
	filterdiag = ones(N+1,1);
	alpha = -log(eps);
	% Initialize filter function
	for i=Nc:N
	    filterdiag(i+1) = exp(-alpha*((i-Nc)/(N-Nc))^s);
	end

function mass = compute_mass(rho,dx)
	rho0 = 1;
	mass = sum( rho - rho0 )*dx;

function mom = compute_mom(rhou,rho,dx)
	
	u = rhou./rho;
	rhouu = rhou.*u;
	mom = 0.5*sum(rhouu)*dx
