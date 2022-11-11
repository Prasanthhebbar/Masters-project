!entering the spin states as input and obtaining bras and kets and innerproduct between each set
program spin_algebra
use Matrix
implicit none
integer::b,k,i,sn,n,bd,kd,j,s,z,h,d,e,w,R,A,l,t
complex::Is=(0.0,1.0)
complex::Iz=(0.0,-1.0)
complex::cn1,cn7,c6,c7,c8,c9
real::z1expectation,xexpectation,z2expectation,sm1sp2expectation,sp1sm2expectation,z1z2expectation
real::c1,c2,c3,c4,c5,zexpectation,l1,l2,l3,l4,concurrence,S_A,S_B,S_AB,pl,pm,trace_rho,cn5,cn11
Integer,Dimension(:),allocatable:: bra,ket,splus,sminus,sz
real,Dimension(:),allocatable::x,o1state,gstate,xstate,o2state,o4state,o6state,o8state,ostate,e_val,y,theta1,theta2,cn6,mi
real,Dimension(:,:),allocatable::Hamiltonian,ev,f,cdm,o3state,o5state,o7state,rho_gs,identity_matrix,rho1,rho2,trans_ev
real,Dimension(:,:),allocatable::mat1,mat2,mat3,mat4,mat5,mat6,mat7,mat8,mat9,mat10,mat11,rho_rhot,rho_B
real,Dimension(:,:),allocatable::sy,rho1t,rhot,dm2,dm3,dm4,dm5,dm6,dm7,dm8,dm9,dm10,rho3,rho4,rho5,rho6,rhoprime_gs,rho_A
real,Dimension(:,:,:),allocatable:: ip,im,isz,O,isx
complex,Dimension(:,:),allocatable::U,rho_ham2,rho_ham1,U_dagger,rho_time1,rho_time2,rho_time3,rho_time,matr1,matr16,matr17,R_ij
complex,Dimension(:,:),allocatable::rho_ham3,rho_ham4,rho_ham5,rho_ham6,rho_time4,rho_time5,rho_time6,rho_time7,rho_time8,rho_time9
complex,Dimension(:,:),allocatable::rho_time10,rho_time11,rho_time12
complex,Dimension(:,:),allocatable::matr2,matr3,matr4,matr5,matr6,matr7,matr8,matr9,matr10,matr11,matr12,matr13,matr14,matr15
complex,Dimension(:),allocatable::V1,V2
real:: abserr=1.0e-02,cn2,cn3,cn4,cn8,cn9,cn10,c10,c11,c12


!print*,"enter the number of spin"
!read*,sn
sn=6
allocate(bra(sn))
allocate(o1state(sn))
allocate(o2state(sn))
allocate(o4state(sn))
allocate(o6state(sn))
allocate(o8state(sn))
allocate(theta1(2**sn))
allocate(theta2(2**sn))
allocate(cn6(100))
allocate(mi(100))
allocate(o3state(2**sn,2**sn))
allocate(o5state(2**sn,2**sn))
allocate(o7state(2**sn,2**sn))
allocate(rho_time1(2**sn,2**sn))
allocate(rho_time2(2**sn,2**sn))
allocate(rho_time3(2**sn,2**sn))
allocate(rho_time(2**sn,2**sn))
allocate(rho_ham1(2**sn,2**sn))
allocate(rho_ham2(2**sn,2**sn))
allocate(rho_ham3(2**sn,2**sn))
allocate(rho_ham4(2**sn,2**sn))
allocate(rho_ham5(2**sn,2**sn))
allocate(rho_ham6(2**sn,2**sn))
allocate(rho_time4(2**sn,2**sn))
allocate(rho_time5(2**sn,2**sn))
allocate(rho_time6(2**sn,2**sn))
allocate(rho_time7(2**sn,2**sn))
allocate(rho_time8(2**sn,2**sn))
allocate(rho_time9(2**sn,2**sn))
allocate(rho_time10(2**sn,2**sn))
allocate(rho_time11(2**sn,2**sn))
allocate(rho_time12(2**sn,2**sn))
allocate(U(2**sn,2**sn))
allocate(U_dagger(2**sn,2**sn))
allocate(gstate(sn))
allocate(ket(sn))
allocate(splus(2**sn))
allocate(V1(2**sn))
allocate(V2(2**sn))
allocate(e_val(2**sn))
allocate(ip(sn,2**sn,2**sn))
allocate(im(sn,2**sn,2**sn))
allocate(isz(sn,2**sn,2**sn))
allocate(isx(sn,2**sn,2**sn))
allocate(O(2**sn,2**sn,2**sn))
allocate(identity_matrix(2**sn,2**sn))
allocate(Hamiltonian(2**sn,2**sn))
allocate(f(2**sn,2**sn))
allocate(rho_gs(2**sn,2**sn))
allocate(ev(2**sn,2**sn))
allocate(trans_ev(2**sn,2**sn))
allocate(rho1(2**sn,2**sn))
allocate(rho2(2**sn,2**sn))
allocate(rho3(2**sn,2**sn))
allocate(rho4(2**sn,2**sn))
allocate(rho5(2**sn,2**sn))
allocate(rho6(2**sn,2**sn))
allocate(rhoprime_gs(2**sn,2**sn))
allocate(mat1(2**sn,2**sn))
allocate(mat2(2**sn,2**sn))
allocate(mat3(2**sn,2**sn))
allocate(mat4(2**sn,2**sn))
allocate(mat5(2**sn,2**sn))
allocate(mat6(2**sn,2**sn))
allocate(mat7(2**sn,2**sn))
allocate(mat8(2**sn,2**sn))
allocate(mat9(2**sn,2**sn))
allocate(mat10(2**sn,2**sn))
allocate(mat11(2**sn,2**sn))
allocate(matr1(2**sn,2**sn))
allocate(matr2(2**sn,2**sn))
allocate(matr3(2**sn,2**sn))
allocate(matr4(2**sn,2**sn))
allocate(matr5(2**sn,2**sn))
allocate(matr6(2**sn,2**sn))
allocate(matr7(2**sn,2**sn))
allocate(matr8(2**sn,2**sn))
allocate(matr9(2**sn,2**sn))
allocate(matr10(2**sn,2**sn))
allocate(matr11(2**sn,2**sn))
allocate(matr12(2**sn,2**sn))
allocate(matr13(2**sn,2**sn))
allocate(matr14(2**sn,2**sn))
allocate(matr15(2**sn,2**sn))
allocate(matr16(2**sn,2**sn))
allocate(matr17(2**sn,2**sn))
allocate(R_ij(4,4))
allocate(dm2(2,2))
allocate(dm3(2,2))
allocate(dm4(2,2))
allocate(dm5(2,2))
allocate(dm6(2,2))
allocate(dm7(2,2))
allocate(dm8(2,2))
allocate(dm9(2,2))
allocate(dm10(2,2))
allocate(rho_A(2,2))
allocate(rho_B(2,2))
allocate(cdm(4,4))
allocate(rho1t(4,4))
allocate(rhot(4,4))
allocate(rho_rhot(4,4))
allocate(sy(4,4))
!allocate(TPD(4,4))
!allocate(dm3(2,2))
!allocate(dmc3(8,8))


do s=1,sn

do b=0,2**sn-1
n=sn-1
bd=b

do i=1,n+1
 bra(sn+1-i)= mod(bd,2)
 bd=bd/2
enddo

do k=0,2**sn-1
 n=sn-1
 kd=k
 do i=1,n+1
 ket(sn+1-i)= mod(kd,2)
 kd=kd/2
 enddo
!print*,(bra(i),i=1,n+1)
!print*,(ket(i),i=1,n+1)

splus=ket         !............for producing splus operator
do j=0,2**sn-1
 if (ket(s)==1) then
 splus(s)=0
  
 else
 
 splus=2

 
 endif
enddo

sminus=ket    !...............!............for producing sminus operator
do j=0,2**sn-1
 if (ket(s)==0) then
 sminus(s)=1
  
 else
 
 sminus=2

 
 endif
enddo

sz=ket       !.....................for producing sz operator

if (ket(s)==0) then
z=1
  
else
 
z=-1

 
endif



!print'(32i3)',(splus(i),i=1,sn)
!ket=splus
ip(s,b+1,k+1)=2*p(bra,splus)
im(s,b+1,k+1)=2*p(bra,sminus)
isz(s,b+1,k+1)=z*p(bra,sz)
enddo

enddo
!print*,"spin plus",s
do i=1,2**sn
!print'(32i3)',(ip(s,i,j),j=1,2**sn)!'(32i3)'
enddo

!print*,"spin minus",s
do i=1,2**sn
!print'(32i3)',(im(s,i,j),j=1,2**sn)!'(32i3)'
enddo

!print*,"spin sz",s
do i=1,2**sn
!print'(32i3)',(isz(s,i,j),j=1,2**sn)!'(32i3)'
enddo

enddo

Hamiltonian=0     !.................................combined hamiltonian matrix
n=sn
if (sn==2) n=1
do h=1,n

d=h+1
if(d>sn) d=1
Hamiltonian=Hamiltonian+0.5 *(matmul(ip(h,:,:),im(d,:,:))+matmul(im(h,:,:),ip(d,:,:)))+matmul(isz(h,:,:),isz(d,:,:))
enddo

!print*,"hamiltonian"
 do i=1,2**sn
! print'(32f8.5)',(hamiltonian(i,j),j=1,2**sn)
 enddo


call Jacobi(hamiltonian,ev,abserr,2**sn)   !..........calling jacobi

!print*,"eigen values are"
do i=1,2**sn                         
!print*,hamiltonian(i,i)     !.........printing eigen values
enddo

!print*,minloc(hamiltonian,mask=hamiltonian<0.0)

!print*,hamiltonian(1,1)
!print*,minloc(1,1)


y=minloc(hamiltonian,mask=hamiltonian<0.0)


!print*,y(1)
w=y(1)

!print*,y(1)

!print*,"eigen vectors are"
!do i=1,2**sn
!print'(32f8.5)',(ev(j,i),j=1,2**sn)!'(32i3)'   !......printing eigenvectors
!enddo

!print*,"ground state energy is"
!print*,minval(hamiltonian)

!print*,"eigen vector corresponding to ground state"
!print'(32f8.5)',(ev(:,w))!'(32i3)'



trans_ev=transpose(ev)



gstate=ev(:,w )

do i=1,2**sn
 do j=1,2**sn
rho_gs(i,j) =gstate(i)*gstate(j)
 enddo
enddo


!print*,"density matrix corresponds to ground state"
 !do i=1,2**sn
!print'(32f8.3)',(rho_gs(i,j),j=1,2**sn)
! enddo

do i=1,2**sn
 do j=1,2**sn
 if (i==j) then
 identity_matrix(i,j)=1
 else
 identity_matrix(i,j)=0
 endif
 enddo
enddo

!print*,"identity matrix"
!do i=1,2**sn
!print'(32f8.5)',(identity_matrix(i,j),j=1,2**sn)!'(32i3)'   !......printing identity matrix
!enddo
 
!10 print*,"enter the spin numbers by which we want to  transform the state" 
!read*,i,j

i=3
j=4

!print*,"enter the value of pl & pm"
!read*,pl,pm
pm=0.5


isx(i,:,:)=0.5*(ip(i,:,:)+im(i,:,:))
isx(j,:,:)=0.5*(ip(j,:,:)+im(j,:,:))


rho1=matmul(isx(i,:,:),isz(j,:,:))
rho2=matmul(rho_gs,rho1)
rho3=matmul(isx(i,:,:),rho2)
rho4=pm*matmul(isz(j,:,:),rho3)
rho5=matmul(identity_matrix,rho_gs)
rho6=(1-pm)*matmul(rho5,identity_matrix)

do i=1,2**sn
e_val(i)=hamiltonian(i,i)
enddo

!print*,"eigen values are"
!do i=1,2**sn
!print*,e_val(i)
!enddo


do i=1,2**sn
V1(i)=0.0
enddo

do i=1,2**sn
V2(i)=0.0
enddo

!print*,"enter the time "
!read*,t
do t=1,100
!t=th/10
!t=0
do j=1,2**sn


V1(j)=exp(-Is*e_val(j)*t*0.05)
V2(j)=exp(Is*e_val(j)*t*0.05)

enddo


!print*,"v array elements"
!do i=1,2**sn
!print*,V1(i)
!enddo

do l=1,2**sn
 do j=1,2**sn
  do k=1,2**sn
O(l,j,k)=ev(j,l)*ev(k,l)
   enddo
  enddo
 enddo


!print*,"eigenvector matrix"
!do i=1,2**sn
!print'(32f8.4)',(O(3,i,j),j=1,2**sn)!'(32i3)'
!enddo

do i=1,2**sn
 do j=1,2**sn
  rho_time(i,j)=(0.0,0.0)
 enddo
enddo

do i=1,2**sn
 do j=1,2**sn
  U(i,j)=(0.0,0.0)
 enddo
enddo

do i=1,2**sn
 do j=1,2**sn
  U_dagger(i,j)=(0.0,0.0)
 enddo
enddo




do j=1,2**sn
U=U+(V1(j)*O(j,:,:))
U_dagger=U_dagger+(V2(j)*O(j,:,:))
enddo

!print*,"U"
!do i=1,2**sn
!print'(32f8.5)',(U(i,j),j=1,2**sn)
!enddo


trace_rho=0.0

rhoprime_gs=rho4+rho6


!do i=1,2**sn
! do j=1,2**sn
!  rho_time=rho_time+exp(Iz*(e_val(i)-e_val(j))*t*0.5)*matmul(O(i,:,:),matmul(rhoprime_gs,O(j,:,:)))
! enddo
!enddo




rho_time5=matmul(U,rhoprime_gs)
rho_time=matmul(rho_time5,U_dagger)


!do i=1,2**sn
! trace_rho=trace_rho+rho_time(i,i)
!enddo
!print*,"trace of transformed density matrix"
!print*,trace_rho


!print*,"transformed time dependent density matrix"
!do i=1,2**sn
!print'(32f8.5)',(rho_time(i,j),j=1,2**sn)!'(32i3)'   !......printing identity matrix
!enddo


!print*,"Enter the spin numbers that we need to find the entanglement and mutual information"
!read*,i,j
i=1
j=2

mat1=matmul(isz(i,:,:),isz(j,:,:))
mat2=matmul(isz(i,:,:),im(j,:,:))
mat3=matmul(isz(i,:,:),ip(j,:,:))
mat4=matmul(im(i,:,:),isz(j,:,:))
mat5=matmul(im(i,:,:),im(j,:,:))
mat6=matmul(im(i,:,:),ip(j,:,:))
mat7=matmul(im(i,:,:),isz(j,:,:))
mat8=matmul(ip(i,:,:),isz(j,:,:))
mat9=matmul(ip(i,:,:),im(j,:,:))
mat10=matmul(ip(i,:,:),ip(j,:,:))
mat11=matmul(ip(i,:,:),isz(j,:,:))

matr1=matmul(rho_time,isz(i,:,:))
matr2=matmul(rho_time,isz(j,:,:))
matr3=matmul(rho_time,mat1)
matr4=matmul(rho_time,im(j,:,:))
matr5=matmul(rho_time,mat2)
matr6=matmul(rho_time,mat3)
matr7=matmul(rho_time,mat4)
matr8=matmul(rho_time,mat5)
matr9=matmul(rho_time,mat6)
matr10=matmul(rho_time,mat7)
matr11=matmul(rho_time,mat8)
matr12=matmul(rho_time,mat9)
matr13=matmul(rho_time,mat10)
matr14=matmul(rho_time,mat11)
matr15=matmul(rho_time,im(i,:,:))
matr16=matmul(rho_time,ip(i,:,:))
matr17=matmul(rho_time,ip(j,:,:))

R_ij(1,1)=0.25*(1.0+trace(matr1)+trace(matr2)+trace(matr3))
R_ij(1,2)=0.25*(trace(matr4)+trace(matr5))
R_ij(1,3)=0.25*(trace(matr15)+trace(matr10))
R_ij(1,4)=0.25*trace(matr8)
R_ij(2,1)=0.25*(trace(matr16)+trace(matr6))
R_ij(2,2)=0.25*(1.0+trace(matr1)-trace(matr2)-trace(matr3))
R_ij(2,3)=0.25*trace(matr9)
R_ij(2,4)=0.25*(trace(matr15)-trace(matr7))
R_ij(3,1)=0.25*(trace(matr16)+trace(matr11))
R_ij(3,2)=0.25*trace(matr9)
R_ij(3,3)=0.25*(1.0-trace(matr1)+trace(matr2)-trace(matr3))
R_ij(3,4)=0.25*(trace(matr4)-trace(matr6))
R_ij(4,1)=0.25*trace(matr13)
R_ij(4,2)=0.25*(trace(matr16)-trace(matr14))
R_ij(4,3)=0.25*(trace(matr17)-trace(matr6))
R_ij(4,4)=0.25*(1.0-trace(matr1)-trace(matr2)+trace(matr3))


!print*,"reduced density matrix is"
!do i=1,4
!print'(32f8.5)',(R_ij(i,j),j=1,4)!'(32i3)'
!enddo

cn1=R_ij(2,3)
cn2=abs(cn1)
cn3=R_ij(1,1)*R_ij(4,4)
cn4=sqrt(cn3)
cn5=cn2-cn4

cn7=R_ij(4,1)
cn8=abs(cn7)
cn9=R_ij(2,2)*R_ij(3,3)
cn10=sqrt(cn9)
cn11=cn7-cn10

cn6(t)=2*max(cn5,0.0,cn11)

!print*,"concurrence"
print*,cn6(t)

 

!!! !!!!!!!!!!MUTUAL INFORMATION !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

C1=R_ij(1,1)
C2=R_ij(2,2)
C3=R_ij(3,3)
C4=R_ij(4,4)
C5=R_ij(1,4)
C6=R_ij(4,1)
C7=R_ij(2,3)
C8=R_ij(3,2)

C9=0.5*(C1+C4)+0.5*(SQRT((C1-C4)*(C1-C4)+4*(abs(C6)*abs(C6))))
C10=0.5*(C1+C4)-0.5*(SQRT((C1-C4)*(C1-C4)+4*(abs(C6)*abs(C6))))
C11=0.5*(C2+C3)+0.5*(SQRT((C2-C3)*(C2-C3)+4*(abs(C7)*abs(C7))))
C12=0.5*(C2+C3)-0.5*(SQRT((C2-C3)*(C2-C3)+4*(abs(C7)*abs(C7))))

S_AB=1.443*(-C9*log(C9)-C10*log(C10)-C11*log(C11)-C12*log(C12))


mi(t)=2.0-S_AB

!print*,"mutual information"
print*,mi(t)

enddo
!print*,"do you want to continue?"
!read*,R
!if(R==1)goto 10

OPEN(unit = 10, access = "sequential", action = "write", status = "replace", file = 'concursix34-12-0.5.csv', form = "formatted")
 do t=1,100
   WRITE(10, '(f10.8,",",f10.8,",")') cn6(t),mi(t)
 enddo
CLOSE(10)



end program spin_algebra











