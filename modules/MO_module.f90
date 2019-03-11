! Nucleon class, w/ attributes N:principal qntmnbr, l:angular mom, j:tot ang mom l±s, ome:(m_J) z-axis projection,
module MO_module 
	implicit none
	integer, parameter :: dp = selected_real_kind(2*kind(1.0))
	integer, parameter :: qp = selected_real_kind(33, 4931)
	real(kind=qp), parameter :: ten_quad = 10._qp

	type, public :: Nucleon
		integer :: N
		integer :: l
		integer :: j
		integer :: ome
		integer :: s
		integer :: q
		real :: E
	end type Nucleon

	public nucleus_creator
contains
	! quicksort - https://gist.github.com/1AdAstra1/6f7785373efe5bb6c254d2e20c78ccc4
	! could be made more general by passing comp-function
	recursive subroutine quicksort_energy(a,first,last)
		!use modified_oscillator
		!implicit none
		type(Nucleon), dimension(:), intent(inout) :: a
		integer,intent(in) :: first, last
		real :: x
		type(Nucleon) :: t
		integer :: i, j

		!last = size(a, 1)

		x = a( (first+last) / 2 )%E
		i = first
		j = last

		do
			do while (a(i)%E < x)
				i=i+1
			end do
			do while (x < a(j)%E)
				j=j-1
			end do
			if (i >= j) exit
			t = a(i);  a(i) = a(j);  a(j) = t
			i=i+1
			j=j-1
		end do

		if (first < i - 1) call quicksort_energy(a,first,i-1)
		if (j + 1 < last) call quicksort_energy(a,j + 1,last)
	end subroutine quicksort_energy

	recursive subroutine quicksort_ome(a,first,last)
		type(Nucleon), dimension(:), intent(inout) :: a
		integer,intent(in) :: first, last
		integer :: x
		type(Nucleon) :: t
		integer :: i, j

		i = first
		j = last

		x = a( (first+last) / 2 )%ome
		!write(*,*) 'X- Pivot element ', x, '| first ', i,'| last ', j, '|'

		do
			do while (a(i)%ome < x .and.  a(i)%E==a(j)%E)
				i=i+1
			end do
			do while (x < a(j)%ome .and.  a(i)%E==a(j)%E)
				j=j-1
			end do
			if (i >= j) exit
			!write(*,*) char(9),'* Putting ' ,a(j)%ome, ' below ', a(i)%ome
			t = a(i);  a(i) = a(j);  a(j) = t
			i=i+1
			j=j-1
		end do

		if (first < i - 1) call quicksort_ome(a,first,i-1)
		if (j + 1 < last) call quicksort_ome(a,j + 1,last)
	end subroutine quicksort_ome


	! Calculates the energy [MeV]
	! KAPPA,MU - values taken from http://www.matfys.lth.se/education/FMF121/The_Nilsson_model_project.pdf …
	! would be better to have shell-specific values, can be found in [Ragnarsson,Nilsson].
	subroutine energy(N,Z,nucl,E)
		!use modified_oscillator
		implicit none

		type(Nucleon),intent(in) :: nucl
		integer, intent(in) :: N,Z
		real, intent(out) :: E
		real :: hw,k,u
		real, dimension(9,2) :: kappa,mu
		integer :: i

		kappa(:,2) = (/0.120,0.120,0.105,0.090,0.065,0.060,0.054,0.054,0.054/) ! Protons
		kappa(:,1) = (/0.120,0.120,0.105,0.090,0.070,0.062,0.062,0.062,0.062/) ! Neutrons
		mu(:,2) = (/0.00,0.00,0.00,0.30,0.57,0.65,0.69,0.69,0.60/) ! Protons
		mu(:,1) = (/0.00,0.00,0.00,0.25,0.39,0.43,0.34,0.26,0.26/) ! Neutrons


		if(nucl%q==0) then
			k=kappa(nucl%N+1,1)
			u=mu(nucl%N+1,1)
		else
			k=kappa(nucl%N+1,2)
			u=mu(nucl%N+1,2)
		end if

		hw = hbar_omega(N,Z,nucl%q)

		if(nucl%s>0) then
			E= hw*( nucl%N +3/2-k*nucl%l-k*u*(nucl%l*(nucl%l+1)-nucl%N*(nucl%N+3)/2) )
		else
			E= hw*( nucl%N +3/2+k*(nucl%l+1)-k*u*(nucl%l*(nucl%l+1)-nucl%N*(nucl%N+3)/2) )
		end if

	end subroutine energy

	! hbar_omega for energy [MeV]
	function hbar_omega(N,Z,q) result(hbOm)
		implicit none
		real :: A

		integer, intent(in) :: N,Z,q
		real :: hbOm

		A=real(N+Z)

		if(q>0) then
			hbom = 41.d0*A**(-1.d0/3.d0)*(1.d0-(1.d0/3.d0)*real(N-Z)/A)
		else
			hbom = 41.d0*A**(-1.d0/3.d0)*(1.d0+(1.d0/3.d0)*real(N-Z)/A)
		end if
	end function hbar_omega

	! Returns list of quantum numbers (Plist:list of particles, see module/class below) for (input) A 'nucleons'
	subroutine qntm_nbrs(n,z,lvls,qCh,Nlist)
		!use modified_oscillator
		implicit none

		integer, intent(in) :: n,z,qCh
		type(Nucleon), dimension(qCh*lvls+(1-qCh)*lvls), intent(out) :: Nlist

		integer :: rem 		! Remaining nucleons to assign
		integer :: ncur 	! Current N under assignment
		integer :: ndeg 	! Nucleons left in current N
		integer :: lcur 	! \ell  currently under assignment
		integer :: j1 		! j1=l+s
		integer :: j2 		! j2=l-s
		integer :: omeCur 	! \Omega currently under assignment
		integer :: i 		! Counter, tracks current nucleon: Nlist(i)
		integer :: n_0,z_0 	! N=8 => z=330,n=330
		integer :: lvls 	! nbr of orbitals/levels to use
		real :: Ener 		! Energy currently under assignment

		n_0=lvls
		z_0=lvls
		ndeg = 0
		ncur = -1
		i = 1
		rem = qCh*z_0+(1-qCh)*n_0

		! Assign quantum numbers to the A nucleons, starting with largest j and omega for each N (ncur)
		wloop:	do while ( rem > 0 )
				if ( ndeg <= 0 ) then
					ncur=ncur+1
					ndeg=(ncur+1)*(ncur+2)
				end if

				if ( ncur==0 ) then
					j1=1 	! j1 = 2*j1 to keep it integer

					do omeCur=-j1, j1, +2
						! ASSIGN QNTM NBRS TO NUCLEON
						Nlist(i)%N 	= ncur
						Nlist(i)%l 	= 0
						Nlist(i)%j 	= j1
						Nlist(i)%ome 	= omeCur
						Nlist(i)%s 	= 1
						Nlist(i)%q 	= qCh

						call energy(n,z,Nlist(i),Ener)

						Nlist(i)%E 	= Ener

						rem = rem - 1
						ndeg = ndeg - 1
						i = i + 1
					end do
				else
					!write(*,*) char(11),'N =',char(9),ncur, char(9), ndeg
					do lcur=ncur, 0, -2

						j1=2*( lcur + 0.5 ) 	! j1 = 2*j1 to keep it integer
						j2=2*( lcur - 0.5 ) 	! j2 = 2*j2 to keep it integer
						!write(*,*) 'l_cur=',char(9), lcur
						!write(*,*) 'j_1 =',char(9),j1
						!write(*,*) 'j_2 =',char(9), j2


						do omeCur=-j1, j1, +2
							! ASSIGN QNTM NBRS TO NUCLEON
							Nlist(i)%N 	= ncur
							Nlist(i)%l 	= lcur
							Nlist(i)%j 	= j1
							Nlist(i)%ome 	= omeCur
							Nlist(i)%s 	= 2*1/2
							Nlist(i)%q 	= qCh

							call energy(n,z,Nlist(i),Ener)

							Nlist(i)%E 	= Ener
							!write(*,*) omeCur

							rem = rem - 1
							ndeg = ndeg - 1
							i = i + 1
							if (rem <= 0) then
								exit wloop
							end if
						end do

						do omeCur=-j2, j2, +2
							! ASSIGN QNTM NBRS TO NUCLEON
							Nlist(i)%N 	= ncur
							Nlist(i)%l 	= lcur
							Nlist(i)%j 	= j2
							Nlist(i)%ome 	= omeCur
							Nlist(i)%s 	= (-2)*1/2
							Nlist(i)%q 	= qCh

							call energy(n,z,Nlist(i),Ener)

							Nlist(i)%E 	= Ener
							!write(*,*) omeCur


							rem = rem - 1
							ndeg = ndeg - 1
							i = i + 1
							if (rem <= 0) then
								exit wloop
							end if
						end do


					if (rem <= 0) then
						exit wloop
					end if

					end do
				end if
			end do wloop


	end subroutine qntm_nbrs

	! Does N and Z
	subroutine nucleus_creator(N,Z,lvls,nucl)
		implicit none
		integer, parameter :: N_tot=330
		integer, intent(in) :: N,Z,lvls
		type(Nucleon), dimension(lvls,2), intent(out) :: nucl
		integer :: i,j

! 		Create the nucleus using N,Z from input
		call qntm_nbrs(N,Z,lvls,+1,nucl(:,2))
		call qntm_nbrs(N,Z,lvls,0,nucl(:,1))
! 		Sort the total nucleus of N=8 (i.e. n,z = 330)
		call quicksort_energy(nucl(:,1),1,lvls)
		call quicksort_energy(nucl(:,2),1,lvls)
		i=1
		do while(i<=lvls)
				j=i+nucl(i,1)%j
				if(j > lvls) j=lvls
				call quicksort_ome(nucl(:,1),i,j)
				i=j+1
		end do
		i=1
		do while(i<=lvls)
				j=i+nucl(i,2)%j
				if(j > lvls) j=lvls
				call quicksort_ome(nucl(:,2),i,j)
				i=j+1
		end do
	end subroutine nucleus_creator

	! Prints the nucleus into terminal and, if opt==1, into files  neutrons/protons.dat
	subroutine print_nucleus(nucl,opt,F1,F2)
		implicit none


		type(Nucleon), dimension(:,:), intent(in) :: nucl
		integer,intent(in) :: opt,F1,F2
		integer :: N,Z,i,prevI,padlen
		integer :: mgNbrs(8)

		N= size(nucl(:,1),1)
		Z= size(nucl(:,2),1)

		if(opt==1) then
			open(unit=F1,file='../data/neutrons.dat')
		       	open(unit=F2,file='../data//protons.dat')
		end if


		prevI = 0
		mgNbrs = (/2,8,20,28,50,82,126,189/)

		100 format (a,i3,a9,i3,i9,i5,i5,a2,i5,a2,i5,a2,f14.2)

		write( *, '(/,31x,a)' ) 'NEUTRONS'
		write( *, '(a22,a5,a5,a6,a8,a6,a18)' ) ' ','N','l','j','m_j','s' , &
			'Energy [MeV]'! 71 positions in total
		padlen=71-4-3

		do i=1,N
			write(*,fmt=100) 'N#:',i,'/#_shell:', i-prevI, &
				nucl(i,1)%N, nucl(i,1)%l, nucl(i,1)%j,'/2', &
				nucl(i,1)%ome,'/2', nucl(i,1)%s ,'/2', nucl(i,1)%E
			if(opt==1) write(F1,*) nucl(i,1)
			if ( ANY( mgNbrs==i ) ) then
				prevI=i
				write( *,  '(a32,a4,i3,a32)') REPEAT( TRIM('-'), floor(padlen/2.0)) ,&
					'Tot=',i,REPEAT( TRIM('-'), ceiling(padlen/2.0))
			end if
		end do

		prevI = 0

		write( *, '(/,31x,a)' ) 'PROTONS'

		write( *, '(a22,a5,a5,a6,a8,a6,a18)' ) ' ','N','l','j','m_j','s' , &
			'Energy [MeV]'! 71 positions in total


		do i=1,Z
			write(*,fmt=100) 'Z#:',i,'/#_shell:', i-prevI, &
				nucl(i,2)%N, nucl(i,2)%l, nucl(i,2)%j,'/2', &
				nucl(i,2)%ome,'/2', nucl(i,2)%s ,'/2', nucl(i,2)%E
			if(opt==1) write(F2,*) nucl(i,2)
			if ( ANY( mgNbrs==i ) ) then
				prevI=i
				write( *,  '(a32,a4,i3,a32)') REPEAT( TRIM('-'), floor(padlen/2.0)) ,&
					'Tot=',i,REPEAT( TRIM('-'), ceiling(padlen/2.0))
			end if
		end do

		if(opt==1) then
			close(F1)
			close(F2)
		end if
	end subroutine print_nucleus

!  ========================== Auxiliary routines ===================================================
!
	! Returns the nucleus of N,Z (neutrons,protons) by extracting from the .dat files created by qntmNbrsMO.f03
	function nucleus_extractor(N,Z) result(nuc)
		integer, intent(in) :: n,z ! nbr of neutrons (n) and protons (z)	
		type(nucleon), dimension(max0(n,z),2) :: nuc
		integer :: i

		! open files from which to extract neutrons,protons
		open(3,file = '../data/neutrons.dat', status='old')
		open(4,file = '../data/protons.dat', status='old')
		
		do i=1,n
			read(3,*) nuc(i,1)
		end do
		close(3)

		do i=1,z
			read(4,*) nuc(i,2)
		end do
		close(4)
	end function nucleus_extractor	

! 	Auxiliary routine: Prints the nucleus 

	subroutine print_n(nucl)
		type(Nucleon), dimension(:,:), intent(in) :: nucl
		integer :: N,Z,i,prevI,padlen
		integer :: mgNbrs(8)

		N= size(nucl(:,1),1)
		Z= size(nucl(:,1),1)


		prevI = 0
		mgNbrs = (/2,8,20,28,50,82,126,189/)

		100 format (a,i3,a9,i3,i9,i5,i5,a2,i5,a2,i5,a2,f14.2) 	

		write( *, '(/,31x,a)' ) 'NEUTRONS' 
		write( *, '(a22,a5,a5,a6,a8,a6,a18)' ) ' ','N','l','j','m_j','s' , &
			'Energy [MeV]'! 71 positions in total 
		padlen=71-4-3

		do i=1,N
			write(*,fmt=100) 'N#:',i,'/#_shell:', i-prevI, &
				nucl(i,1)%N, nucl(i,1)%l, nucl(i,1)%j,'/2', & 
				nucl(i,1)%ome,'/2', nucl(i,1)%s ,'/2', nucl(i,1)%E
			if ( ANY( mgNbrs==i ) ) then 
				prevI=i
				write( *,  '(a32,a4,i3,a32)') REPEAT( TRIM('-'), floor(padlen/2.0)) ,&
					'Tot=',i,REPEAT( TRIM('-'), ceiling(padlen/2.0))
			end if
		end do

		prevI = 0

		write( *, '(/,31x,a)' ) 'PROTONS' 

		write( *, '(a22,a5,a5,a6,a8,a6,a18)' ) ' ','N','l','j','m_j','s' , &
			'Energy [MeV]'! 71 positions in total 


		do i=1,Z
			write(*,fmt=100) 'Z#:',i,'/#_shell:', i-prevI, &
				nucl(i,2)%N, nucl(i,2)%l, nucl(i,2)%j,'/2', & 
				nucl(i,2)%ome,'/2', nucl(i,2)%s ,'/2', nucl(i,2)%E
			if ( ANY( mgNbrs==i ) ) then 
				prevI=i
				write( *,  '(a32,a4,i3,a32)') REPEAT( TRIM('-'), floor(padlen/2.0)) ,&
					'Tot=',i,REPEAT( TRIM('-'), ceiling(padlen/2.0))
			end if
		end do
	end subroutine print_n

end module MO_module

