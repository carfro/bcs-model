%% Prototyping for particle projection - Onishi formula
% Have problems with the fortran-code getting numerical errors, 
% determinant outputs 0 sometimes?
close all,clear all%,clf

U       = dlmread('U_mat.dat');
V       = dlmread('V_mat.dat');
DIM 	= length(U);

% Loop variables:nbr of points, expected part. nbr. & phi-step
N_loop 	= 10;
L       = 2*N_loop+1
N_exp 	= 2*sum(nonzeros(diag(V,1).^2))
dPhi 	= 2*pi/(2*N_loop+1);
W       = 6;

%Pfaffian init 
WTW 	= zeros(2*DIM);
nor 	= prod(nonzeros(diag(V,1).^2))

S 	= (-1)^(DIM*(DIM-1)*0.5);
S_norm 	= S/nor;

%Init of phi onishi and pfaffian vectors
phi 	= zeros(L,1);
oniA 	= zeros(L,1);
pfA 	= zeros(2*N_loop+1,1);
diffA 	= zeros(2*N_loop+1,1);

oniR 	= zeros(2*N_loop+1,1);
pfR 	= zeros(2*N_loop+1,1);
diffR 	= zeros(2*N_loop+1,1);

oniI 	= zeros(2*N_loop+1,1);
pfI 	= zeros(2*N_loop+1,1);
diffI 	= zeros(2*N_loop+1,1);

%Init of overlap vectors
ROL 	= zeros(2*N_loop+1,1);
ROL2 	= zeros(2*N_loop+1,1);
NROL 	= zeros(2*N_loop+1,1);
NROL2 	= zeros(2*N_loop+1,1);

Rsum 	= zeros(2*W+1,1);
Rsum2 	= zeros(2*W+1,1);
Nexp 	= zeros(2*W+1,1);
Nexp2 	= zeros(2*W+1,1);
Nvec 	= zeros(2*W+1,1);

traceTest = zeros(2*W+1,1);
oltest = zeros(2*W+1,1);

for j=ceil(N_exp)-W:ceil(N_exp)+W
%for j=ceil(N_exp):ceil(N_exp)
	%j
     ind2 = j-(ceil(N_exp)-W)+1;
% 	ind2=1;
    Nvec(ind2) = j;

    NROL 	= zeros(2*N_loop+1,1);
    NROL2 	= zeros(2*N_loop+1,1);
	
    for i=0:2*N_loop
		ind1=i+1+N_loop;
		D=eye(DIM).*exp(i*1i*dPhi);
		phi(ind1)=i*dPhi;

		%Onishi-formula
		X=ctranspose(U)*ctranspose(D)*U + ctranspose(V)*transpose(D)*V;
        %Xflip=ctranspose(U)*(D)*U + ctranspose(V)*conj(D)*V;

		%oni=sqrt(det(X))*exp(-1i*(1)*dPhi*i*0.5);
		oni=sqrt(det(X))*exp(-1i*(DIM*0.5)*dPhi*i);
		%	    oniA(ind1)=abs(oni);
		%	    oniR(ind1)=real(oni);
		%	    oniI(ind1)=imag(oni);     
		ROL2(ind1)=dPhi*exp(-i*1i*j*dPhi*0.5)*oni;     

		%Pfaffian calculation
		WTW=[V.'*U,V.'*D*conj(V);-1*ctranspose(V)*D.'*V,ctranspose(U)*conj(V)];
		pf=S_norm*pfaffian_householder(WTW);
		%	    pfA(ind1)=abs(pf);
		%	    pfR(ind1)=real(pf);
		%	    pfI(ind1)=imag(pf);
		ROL(ind1)=dPhi*exp(-i*1i*j*dPhi)*pf;

		%	    diffA(ind1)=(abs(oni-pf));
		%	    diffR(ind1)=abs(real(oni))-abs(real(pf));
		%	    diffI(ind1)=abs(imag(oni))-abs(imag(pf));

        X           	= 	D*conj(V)*inv(X.')*V.';
        %Xflip          	= 	conj(V)*inv(Xflip.')*V.'*ctranspose(D);

%         X2              =   (X-ctranspose(Xflip));
%         traceTest(ind1) = sqrt(trace(ctranspose(X2)*X2));
%         oltest(ind1)    = abs(pfaffian_householder(WTW)-pfaffian_householder(WTW.'));
        
        TrcRho          = 	trace(X);

		NROL(ind1)  	= 	ROL(ind1)*TrcRho;
		NROL2(ind1) 	= 	ROL2(ind1)*TrcRho;
        
	end

	Rsum(ind2)= 1/(2*pi)*sum(ROL);	
	Rsum2(ind2)= 1/(2*pi)*sum(ROL2);
	Nexp(ind2)= sum(NROL)/(sum(ROL));
	Nexp2(ind2)= sum(NROL2)/sum(ROL2);

end
% traceTest
% oltest

%figure;
%plot(phi,oniA)
%hold on;
%plot(phi,pfA)
%%plot(phi,diffA)
%legend({'{Abs(Onishi)}','{Abs(Pfaffian)}','{Abs(difference)}'})
%
%figure;
%plot(phi,oniR)
%hold on;
%plot(phi,pfR)
%%plot(phi,diffR)
%legend({'{Re(Onishi)}','{Re(Pfaffian)}','{Re(difference)}'})
%
%figure;
%plot(phi,oniI)
%hold on;
%plot(phi,pfI)
%%plot(phi,diffI)
%legend({'{Im(Onishi)}','{Im(Pfaffian)}','{Im(difference)}'})


figure('Name',      'Rotade overlaps'       ,'NumberTitle','off');
plot(Nvec,real(Rsum))
hold on;
plot(Nvec,real(Rsum2))
legend({'{Re(Overlap_{Pfaffian})}','{Re(Overlap_{Onishi})}'})

% figure('Name',      'Rotade overlaps,imag'       ,'NumberTitle','off');
% plot(Nvec,imag(Rsum))
% hold on;
% plot(Nvec,imag(Rsum2))
% legend({'{Im(Overlap_{Pfaffian})}','{Im(Overlap_{Onishi}})'})

figure('Name',      'Particle number'       ,'NumberTitle','off');
plot(Nvec,Nexp)
hold on;
plot(Nvec,Nexp2)
legend({'{<N>_{Pfaffian}}','{<N>_{Onishi}}'})

% figure;
% plot(linspace(-N_loop,N_loop),traceTest)
% hold on;
% plot(linspace(-Nloop,Nloop),oltest)

