%% Prototyping for particle projection - Onishi formula & pfaffian
close all,clf, clear all

U      = dlmread('U_mat_025d.dat');
V      = dlmread('V_mat_025d.dat');
DIM 	= length(U);
Ul      = dlmread('U_mat_4d.dat');
Vl      = dlmread('V_mat_4d.dat');
Us      = dlmread('U_mat_00625d.dat');
Vs      = dlmread('V_mat_00625d.dat');
% Loop variables:nbr of points, expected part. nbr. & phi-step
N_loop 	= 10*1;
L       = 2*N_loop+1;
N_exp 	= 2*sum(nonzeros(diag(V,1).^2))
dPhi 	= 2*pi/(2*N_loop+1);
W       = 6;

%Pfaffian init 
WTW 	= zeros(2*DIM);
nor 	= prod(nonzeros(diag(V,1).^2))
S_norm 	= (-1)^(DIM*(DIM-1)*0.5)/nor;

nor_l 	= prod(nonzeros(diag(Vl,1).^2))
S_norm_l	= (-1)^(DIM*(DIM-1)*0.5)/nor_l;

nor_s 	= prod(nonzeros(diag(Vs,1).^2))
S_norm_s	= (-1)^(DIM*(DIM-1)*0.5)/nor_s;
%Init of phi onishi and pfaffian vectors
phi 	= zeros(L,1);
oniV 	= zeros(L,1);
pfV 	= zeros(L,1);
pfV2 	= zeros(L,1);
pfV3 	= zeros(L,1);
diffV 	= zeros(L,1);

%Init of overlap vectors
ROL 	= zeros(L,1);
ROL2 	= zeros(L,1);
ROL3 	= zeros(L,1);
NROL 	= zeros(L,1);
NROL2 	= zeros(L,1);

Rsum 	= [];
Rsum2 	= [];
Rsum3 	= [];
Nexp 	= [];
Nexp2 	= [];
Nvec 	= [];

%oltest   = zeros(2*W+1,1);

%figure;
%hold on;

for j=round(N_exp)-W:1:round(N_exp)+W
%for j=ceil(N_exp):ceil(N_exp)
    Nvec = [Nvec;j];

    NROL 	= zeros(L,1);
    NROL2 	= zeros(L,1);
    traceTest   = zeros(L,1);
	
    for i=0:2*N_loop
		ind1 =i+1;
        	D=eye(DIM).*exp(i*1i*dPhi);
		phi(ind1)=i*dPhi;

		%Onishi-formula
		X=ctranspose(U)*ctranspose(D)*U + ctranspose(V)*transpose(D)*V;
		Xl=ctranspose(Ul)*ctranspose(D)*Ul + ctranspose(Vl)*transpose(D)*Vl;
		%%oni=sqrt(det(X))*exp(-1i*(DIM*0.5)*dPhi*i);
		%if j==round(N_exp)-W oniV(ind1)=sqrt(det(X))*exp(-1i*(DIM*0.5)*dPhi*i);, end;
		%ROL2(ind1)=dPhi*exp(-i*1i*j*dPhi*0.5)*oniV(ind1);     

		%Pfaffian calculation
		WTW=[V.'*U,V.'*D*conj(V);-1*ctranspose(V)*D.'*V,ctranspose(U)*conj(V)];
		if j==round(N_exp)-W pfV(ind1)=S_norm*pfaffian_householder(WTW);, end;
		ROL(ind1)=dPhi*exp(-i*1i*j*dPhi)*pfV(ind1);
		
		WTW=[Vl.'*Ul,Vl.'*D*conj(Vl);-1*ctranspose(Vl)*D.'*Vl,ctranspose(Ul)*conj(Vl)];
		if j==round(N_exp)-W pfV2(ind1)=S_norm_l*pfaffian_householder(WTW);, end;
		ROL2(ind1)=dPhi*exp(-i*1i*j*dPhi)*pfV2(ind1);

		WTW=[Vs.'*Us,Vs.'*D*conj(Vs);-1*ctranspose(Vs)*D.'*Vs,ctranspose(Us)*conj(Vs)];
		if j==round(N_exp)-W pfV3(ind1)=S_norm_s*pfaffian_householder(WTW);, end;
		ROL3(ind1)=dPhi*exp(-i*1i*j*dPhi)*pfV3(ind1);

		X           	= 	D*conj(V)*inv(X.')*V.';
		Xl           	= 	D*conj(Vl)*inv(Xl.')*Vl.';
		
		%Xflip         	 = 	conj(V)*inv(Xflip.')*V.'*ctranspose(D);
	        %X2              = 	(X-ctranspose(Xflip));
	        %traceTest(ind1) =  	sqrt(trace(ctranspose(X2)*X2));
	        %oltest(ind1)    =  	abs(pfaffian_householder(WTW)-pfaffian_householder(WTW.'));

	        %traceTest(ind1) =  	trace(X);

		NROL(ind1)  	= 	ROL(ind1)*trace(X);
		NROL2(ind1) 	= 	ROL2(ind1)*trace(Xl);
	end
    %plot(phi,real(traceTest),phi,real(ROL))

    Rsum = [Rsum; 1/(2*pi)*sum(ROL)];
    Rsum2= [Rsum2; 1/(2*pi)*sum(ROL2)];
    Rsum3= [Rsum3; 1/(2*pi)*sum(ROL3)];

    Nexp = [Nexp; sum(NROL)/(sum(ROL))];
    Nexp2= [Nexp2; sum(NROL2)/(sum(ROL2))];
    Nexp(find(Rsum < 10^-4))=0;
end

%figure;
%plot(phi,traceTest,phi,ROL)
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
hold on;
plot(Nvec,real(Rsum),Nvec,real(Rsum2),Nvec,real(Rsum3))
title('$\sum_l <\phi|R_l|\phi>$','interpreter','latex')
legend({'{0.25*\delta_0}','{4*\delta_0}','{0.0625*\delta_0}'})

% figure('Name',      'Rotade overlaps,imag'       ,'NumberTitle','off');
% plot(Nvec,imag(Rsum))
% hold on;
% plot(Nvec,imag(Rsum2))
% legend({'{Im(Overlap_{Pfaffian})}','{Im(Overlap_{Onishi}})'})

figure('Name',      'Particle number'       ,'NumberTitle','off');
plot(Nvec,Nexp)
hold on;
%plot(Nvec,Nexp2)
legend({'{<N>_{Pfaffian}}','{<N>_{Onishi}}'})%,'{<N>_{no trace}}'})

