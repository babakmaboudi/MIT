function generate_reduced_basis()
	load('snaps.mat')
	
	snaps = [RHO , RHOU];
	
	[U,S,V] = svd(snaps);

%	semilogy(diag(S))

	A = U(:,1:50);

	save('reduced_basis.mat','A')
