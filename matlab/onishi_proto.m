%% Prototyping for particle projection - Onishi formula
% Have problems with the fortran-code getting numerical errors, 
% determinant outputs 0 sometimes?


U=dlmread('U_mat.dat');
V=dlmread('V_mat.dat');
DIM=length(U)

% Loop variables:no points, expected part. nbr. & phi-step
N_loop=10*30*30
N_exp=2*diag(V,1).'*diag(V,1)
dPhi=8*(4*pi/((DIM+N_exp)*N_loop));

%Init of phi and onishi vectors
phi=zeros(DIM,1);
oni=zeros(DIM,1);

for j=0:N_loop
    D=eye(DIM).*exp(1i*dPhi);
    X=ctranspose(U)*ctranspose(D)*U + ctranspose(V)*transpose(D)*V;

    phi(j+1)=j*dPhi;
    oni(j+1)=real(sqrt(det(X))*exp(-1i*DIM*dPhi*j));
end
   

plot(phi,oni)