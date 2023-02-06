Module mod_afssh
use matmul_lapack
!! Hammes-Schiffer, Tully, JCP 101, 4657 (1994)
implicit none
real*8, parameter :: clight=2.99792458D10,av=6.0221367D23,hbar=1.05457266D-34
real*8, parameter :: kb=1.3806503d-23
real*8, parameter :: pascal_to_atm=9.86923267d-6,kcal_to_J=6.95d-21
real*8, parameter :: amu2kg=1.66053892d-27
real*8, parameter :: au2kg=9.10938291d-31,au2J=4.35974417d-18,au2m=5.2917721092d-11,au2s=2.418884326505d-17
real*8, parameter :: q_el=2.307075d-28
complex*16,parameter :: iota = (0.d0,1.d0)
real*8 pi,wave_to_J

!! Potential
real*8 A_pot(2),C_pot(2),r_MeCl  !! 1-methyl 2-chloride; solvent-solvent potential;
real*8 a_solute,b_solute,c_solute,Diss_A,d_A,d_B,n_A,n_B !!phenol-trimethylamine; solute-solute interaction
real*8 sig_meth_meth,sig_meth_cl,sig_cl_cl,sig_sol_com
real*8 eps_meth_meth,eps_meth_cl,eps_cl_cl,eps_sol_com
real*8 LJ_cut_meth_meth,LJ_cut_meth_cl,LJ_cut_cl_cl,LJ_cut_sol_com
real*8 dLJ_cut_meth_meth,dLJ_cut_meth_cl,dLJ_cut_cl_cl,dLJ_cut_sol_com
real*8 e_alpha_C(3),e_alpha_i(3),r0,l0  !! solute solvent columbic parameters
real*8 R_c,R_T,Rc_verlet !! cutoffs for LJ and columbic
real*8 box_length(3)
real*8 temperature_init,temperature,gamma_B
real*8 mass_H,mass_MeCl
real*8,allocatable :: mass(:),charge(:)
real*8 charge_H
character*2,allocatable :: atom(:)
real*8 r_cutoff

!! Output/Input
real*8 k2_val
real*8 cnt_frust,cnt_collapse,cnt_init,cnt_term
real*8,allocatable :: pop(:,:),pop_surf(:,:),pop_amp(:,:)

!! Classical
integer nclass,idistribution
real*8 rAH,rAB
real*8,allocatable :: x(:,:),v(:,:),acc(:,:)
real*8,allocatable :: x_old(:,:),v_old(:,:),acc_old(:,:)
real*8,allocatable :: x_in(:,:),v_in(:,:),x_start(:,:)
real*8,allocatable :: distance(:,:)
integer,allocatable :: verlet_list(:,:)
integer iverlet_list
real*8 alpha
integer iforward,flag_terminate,flag_frust,flag_reactant,flag_classical,flag_hop
complex*16,allocatable :: delr(:,:,:),delp(:,:,:)
complex*16,allocatable :: delr_old(:,:,:),delp_old(:,:,:)

!! Quantum
integer nquant,state,nbasis,state_tentative
integer state_old
real*8,allocatable :: si_adiab(:,:),V_k(:),d_ij(:,:,:,:),vdotd(:,:),V_k_old(:)
real*8,allocatable :: Hamil_diab(:,:),delH_dels(:,:,:,:),delH_dels_ad(:,:,:,:)
real*8,allocatable :: pot(:,:),force(:,:,:),force_old(:,:,:)
complex*16,allocatable :: ci(:),ci_old(:)
real*8,allocatable :: si_adiab_prev(:,:)
complex*16,allocatable :: mat(:,:),mat_adiab(:,:)
real*8,allocatable :: hop_prob(:),W_overlap(:,:),hop_prob_net(:)
real*8,allocatable :: KE_DVR(:,:)
real*8 rAH_min,rAH_max,rAH_del,r_exp
real*8,allocatable:: rAH_grid(:)

!! Evolution
integer n_traj,nsteps,nsteps_eq,nstep_write,iwrite,nstep_avg
real*8 dtc,total_time,curr_time,traj_num,tim_eq
real*8 energy,pot_en,KE_en,energy_cutoff,energy_old
real*8 ensq_avg,en_avg
integer nst_av
integer ihop,icollapse,iterminate,ifriction,iaverage
real*8 prob_tun

!! Parallelization
integer iparallel,iflow,iproc,iwait,ifolder

!! Misc
integer nold,cnt_rate
real*8 tim_tot,tim_ev_cl,tim_diag,tim_cl,tim_rattle,tim_pbc,tim_LJ_tot,tim_solv_solv,tim_check,tim_check2
real*8 tim_T_jk
integer,allocatable:: seed(:)
real*8,allocatable:: work(:)
complex*16,allocatable:: cwork(:)
integer,allocatable:: iwork(:),isuppz(:)

contains
!---------------------------------------------------------- 
!---------------------------------------------------------- 

subroutine setup
  implicit none
  character st_ch
  integer i,size_seed,seed2(2)
  real*8 rnd

  pi=dacos(-1.d0)
  wave_to_J=2*pi*clight*hbar
  nold=0

  open(10,file="AFSSH.inp")
  read(10,*) iflow
  read(10,*) iproc
  read(10,*) iparallel
  read(10,*) iwait
  read(10,*) N_traj
  read(10,*) dtc
  read(10,*) tim_eq
  read(10,*) total_time
  read(10,*) iwrite
  read(10,*) nstep_write
  read(10,*) nstep_avg
  read(10,*) iverlet_list
  read(10,*) idistribution
  read(10,*) flag_frust
  read(10,*) energy_cutoff
  read(10,*) nclass
  read(10,*) nquant
  read(10,*) nbasis
  read(10,*) rAH_min
  read(10,*) rAH_max
  read(10,*) iforward
  read(10,*) icollapse
  read(10,*) ifriction
  read(10,*) temperature_init
  read(10,*) seed2
  read(10,*) st_ch
  close(10)
  !----------------------------------------------------------

  if(st_ch.ne.'x') then
    write(6,*) "problem in reading input file"
    stop
  endif
  !---------------------------------------------------------- 

  energy_cutoff=energy_cutoff*wave_to_J

  nsteps=nint(total_time/dtc)+1
  nsteps_eq=nint(tim_eq/dtc)+1

  !-----------------------------------------------------------------  
  i=nsteps/nstep_avg+1
  allocate(pop(nquant,i),pop_surf(nquant,i),pop_amp(nquant,i))
  allocate(x(nclass,3),v(nclass,3),acc(nclass,3))
  allocate(x_old(nclass,3),v_old(nclass,3),acc_old(nclass,3))
  allocate(mass(nclass),charge(nclass),atom(nclass),verlet_list(nclass,nclass))!,distance(nclass,nclass))
  allocate(KE_DVR(nbasis,nbasis),rAH_grid(nbasis))
  allocate(delr(nquant,nclass,3),delp(nquant,nclass,3))
  allocate(delr_old(nquant,nclass,3),delp_old(nquant,nclass,3))
  allocate(si_adiab(nbasis,nquant),ci(nquant),V_k(nquant),V_k_old(nquant))
  allocate(Hamil_diab(nbasis,nbasis),delH_dels(nbasis,nbasis,nclass,3),delH_dels_ad(nquant,nquant,nclass,3))
  allocate(pot(nquant,nquant),force(nquant,nclass,3),force_old(nquant,nclass,3))
  allocate(mat(nbasis,nbasis),mat_adiab(nquant,nquant))
  allocate(d_ij(nquant,nquant,nclass,3),vdotd(nquant,nquant),hop_prob(nquant),W_overlap(nquant,nquant))
  allocate(hop_prob_net(nquant))
  allocate(ci_old(nquant),si_adiab_prev(nbasis,nquant))
  call random_seed(size=size_seed)
  allocate(seed(size_seed))
  do i=1,size_seed/2
    seed(i)=seed2(1)*(2*i+i*i-7)
  enddo
  do i=size_seed/2+1,size_seed
    seed(i)=seed2(2)*(i/2+34-i**3)
  enddo
  call random_seed(put=seed)
  call system_clock(count_rate=cnt_rate)
  !-----------------------------------------------------------------  

  if(iflow==2) then
    open(10,file="ifolder.inp")
    read(10,*) ifolder
    close(10)
    N_traj=N_traj/iparallel
    if(ifolder>1) then
      do i=1,(ifolder-1)*N_traj
        seed=seed+1
        call init_cond
      enddo
      call random_seed(put=seed)
    endif
  else
    ifolder=1
  endif

  tim_tot=0.d0

end subroutine setup
!---------------------------------------------------------- 

subroutine main
  implicit none
  integer i,j,k,n
  real*8 t1,t2
  integer if_reactant

  call files(0)

  call cpu_time(t1)

  call setup_parameters
  call initialize_averages

  do i=1,N_traj
    traj_num=i
    call init_cond
    call equilibrate(if_reactant)
    if(if_reactant==1) then
      cnt_init=cnt_init+1
      iaverage=1;en_avg=0.d0;ensq_avg=0.d0
      call evolve(nsteps)
    endif
  enddo
  call write_average

  call cpu_time(t2);tim_tot=tim_tot+t2-t1
  call files(1)

end subroutine main
!---------------------------------------------------------- 

subroutine files(flag)
  implicit none
  integer,intent(in)::flag

  if(flag==0) then
    open(10,file="output")
    open(11,file="output_cl")
    open(12,file="output_qm")
    open(13,file="output_hop")
    open(14,file="output_overlap")
    open(15,file="output_dec")
    open(50,file="traj_vmd.xyz")

    open(100,file="pop.out")
    open(101,file="cnts.out")
  else
    write(10,*)
    write(10,*)"Total time=",tim_tot
    close(10);close(11);close(12);close(13);close(14);close(15)
    close(100);close(101)
  endif

end subroutine files
!-----------------------------------------------------------------  

subroutine initialize_averages
  implicit none

  cnt_frust=0.d0
  cnt_collapse=0.d0
  cnt_init=0.d0
  cnt_term=0.d0
  pop=0.d0
  pop_surf=0.d0
  pop_amp=0.d0

end subroutine initialize_averages
!-----------------------------------------------------------------  

subroutine init_cond
  implicit none
  real*8 beta
  real*8 fac,sig_x,sig_p,rnd
  integer i,j,n,nc,sgn
  real*8 bond_len
  real*8,allocatable,dimension(:,:) :: rr,ee
  real*8 rij,xij(3)

  !! Boltzmann distribution
  energy=0.d0
  V_k(2)=1.d0

  n=nclass/2
  nc=nint((n/4)**(1/3.d0))

  allocate(rr(n,3),ee(n,3))
  call fcc(rr(:,1),rr(:,2),rr(:,3),ee(:,1),ee(:,2),ee(:,3),n,nc)

  sgn=1
  do i=1,nclass/2
    if(i==1) bond_len=2.7d-10
    if(i>1) bond_len=sgn*r_MeCl
    sgn=sgn*(-1)
    j=2*i-1
    x(j,:)=rr(i,:)-ee(i,:)*bond_len*mass(j+1)/(mass(j)+mass(j+1))!  /2.d0
    j=2*i
    x(j,:)=rr(i,:)+ee(i,:)*bond_len*mass(j-1)/(mass(j-1)+mass(j))     !/2.d0
  enddo

  call pbc_adjust
  state=1
  flag_classical=1
  call Boltzmann_velocities(v,nclass,mass,temperature_init)
  !v=0.d0

  x_in=x;v_in=v

  alpha=0
!  kappa_num_tr=0.d0

  ci=0.d0
  ci(State)=1

  call update_verlet_list
  flag_classical=1
  ihop=1
  iaverage=1
  iterminate=0
  flag_terminate=0

  curr_time=0.d0
  call evaluate_variables(0)
  call evaluate_variables(1)
  call compute_mat_diab

  en_avg=0.d0;ensq_avg=0.d0
  nst_av=0

end subroutine init_cond
!-----------------------------------------------------------------  

subroutine equilibrate(if_reactant)
  implicit none
  integer,intent(out)::if_reactant
  integer i,n,nst
  real*8 frac

  call set_charges(1.d0)

  ifriction=0!1
  iaverage=0
  en_avg=0.d0;ensq_avg=0.d0
  call evolve(nsteps_eq)

!  ifriction=0
!  n=5
!  nst=nint(nsteps_eq/real(n))/2
!
!  do i=1,n
!!    call Boltzmann_velocities(v,nclass,mass,temperature_init)
!    frac=i/real(n)
!    !call set_charges(frac)
!    call set_charges(1.d0)
!    en_avg=0.d0;ensq_avg=0.d0
!    iaverage=0
!    call evolve(nst)
!write(6,*) frac,r_exp*1.d10,temperature
!  enddo
!
  call check_reactant(if_reactant)
  call set_charges(1.d0)

end subroutine equilibrate
!-----------------------------------------------------------------  

subroutine evolve(nsteps)
  implicit none
  integer,intent(in) :: nsteps
  integer i,j,nstep_sm,iflag_coll,i_do_something
  real*8 t1,t2
  integer iterm

  !call cpu_time(t1)

  call write_output(1,1)
  do i=1,nsteps
    call write_output(i,0)
    !if(mod(i,100)==1) call write_vmd
    call average(i)
    call save_old_state
    if(mod(i,iverlet_list)==0) call update_verlet_list
    call evolve_classical(dtc)
    i_do_something=0
    if(ifriction==0) then
      do while(dabs(energy-energy_old)>energy_cutoff.and.i>1) 
        i_do_something=i_do_something+1
        call do_something(i_do_something)
      enddo
    endif
    if(i_do_something.ne.1)call evolve_quantum_small_dtq
    !call evolve_quantum_small_dtq
    if(i_do_something.ne.1.and.ihop==1)call hop
    if(i_do_something.ne.1.and.icollapse==1)call collapse(dtc,iflag_coll)
    if(flag_terminate==1) call traj_terminate(iterm)
      if(iterm==1)exit

    curr_time=curr_time+dtc
  enddo
  call write_output(1,1)

  !call cpu_time(t2)
  !tim_evolve=tim_evolve+t2-t1

end subroutine evolve
!-----------------------------------------------------------------  

subroutine do_something(i_do_something)
  implicit none
  integer,intent(in)::i_do_something
  real*8 dt,dtm(3)
  integer i,nstep

write(6,*) curr_Time*1.d15,i_do_something,dt*1.d15,(energy-energy_old)/wave_to_J

  if(i_do_something==1) then
    call evolve_quantum_small_dtq
    if(flag_hop==1) then
      state=state_tentative
      call evaluate_variables(0)
      v=v_old+0.5*(acc_old+acc)*dtc
      call evaluate_variables(1)
    endif
  else
    dtm=1.d0 !! some large value
    dtm(1)=0.1d0/maxval(vdotd)
    dtm(2)=0.5*dtc*dsqrt(energy_cutoff/dabs(energy-energy_old))
    dtm(3)=dtc

    dt=minval(dtm)

    dt=dt/real(i_do_something)
    nstep=nint(dtc/dt)
    dt=dtc/real(nstep)
    call revert_state

    do i=1,nstep
      call evolve_classical(dt)
    enddo
    !if(dabs(energy-energy_old)<energy_cutoff) write(6,*)curr_time*1.d15,dt*1.d15
  endif


end subroutine do_something
!-----------------------------------------------------------------  

subroutine average(i)
  implicit none
  integer,intent(in) :: i
  integer j
  complex*16 ci_diab(nquant)
  integer if_reactant
  real*8 t1,t2

  !call cpu_time(t1)

  if(iwrite==1) then
    en_avg=en_avg+energy
    ensq_avg=ensq_avg+energy*energy
    nst_av=nst_av+1
  endif

  if(iaverage==1.and.(mod(i,nstep_avg)==1.or.nstep_avg==1)) then
    if(nstep_avg==1) then
      j=i
    else
      j=i/nstep_avg+1
    endif

!    pop(:,j)=pop(:,j)+si_adiab(:,state)**2
!    pop_surf(:,j)=pop_surf(:,j)+si_adiab(:,state)**2
!    pop(:,j)=pop(:,j)+2*realpart(ci(1)*dconjg(ci(2)))*si_adiab(:,1)*si_adiab(:,2)
!    pop_amp(1,j)=pop_amp(1,j)+cdabs(sum(si_adiab(1,:)*ci))**2
!    pop_amp(2,j)=pop_amp(2,j)+cdabs(sum(si_adiab(1,:)*ci))**2

    call check_reactant(if_Reactant)
    if(if_reactant==1)pop(1,j)=pop(1,j)+1

  endif

  !call cpu_time(t2)
  !tim_coll=tim_coll+t2-t1

end subroutine average
!-----------------------------------------------------------------  

subroutine average_end
  implicit none

end subroutine average_end
!-----------------------------------------------------------------  

subroutine check_reactant(if_reactant)
  implicit none
  integer,intent(out)::if_reactant

  if_reactant=0
  if(r_exp<r_cutoff) if_reactant=1

end subroutine check_reactant
!-----------------------------------------------------------------  

subroutine save_old_state
  implicit none

  x_old=x
  v_old=v
  acc_old=acc
  ci_old=ci
  state_old=state
  !ci2_old=ci2
  si_adiab_prev=si_adiab
  V_k_old=V_k
  force_old=force
  energy_old=energy
  delr_old=delr
  delp_old=delp

end subroutine save_old_state
!-----------------------------------------------------------------  

subroutine revert_state
  implicit none

  x=x_old
  v=v_old
  state=state_old
  ci=ci_old
  delr=delr_old
  delp=delp_old
  force=force_old
  !ci2=ci2_old
  call evaluate_variables(0)
  call evaluate_variables(1)

end subroutine revert_state
!-----------------------------------------------------------------  

subroutine evolve_quantum_exponential
  !! Basic quantum propagator
  !! CAUTION AMBER !!
  !! delr and delp not implmented!!
  implicit none
  integer i,j
  real*8 W_overlap(nquant,nquant),si_adiab_old(nbasis,nquant),V_k_old(nquant)

  si_adiab_old=si_adiab
  V_k_old=V_k

!  call evolve_classical(dtc)

  do i=1,nquant
    do j=1,nquant
      W_overlap(i,j)=sum(si_adiab(:,i)*si_adiab_old(:,j))
    enddo
  enddo

  ci=ci*cdexp(-iota*V_k_old*dtc/hbar)
  ci=matmul(W_overlap,ci)

end subroutine evolve_quantum_exponential
!-----------------------------------------------------------------  

subroutine evolve_quantum_small_dtq
  implicit none
  integer i,nstep_qm
  real*8 dtq,dtq1,dtq2
  real*8 V_k_hold(nquant),dVk_dt(nquant)
  real*8 dforce_dt(nquant,nclass)

  call compute_vdotd
  dVk_dt=(V_k-V_k_old)/dtc
  if(icollapse==1) then
    call compute_delH_dels_ad
  endif

  dtq1=0.002/maxval(vdotd)
  dtq2=0.002*hbar/maxval(V_k-sum(V_k)/real(nquant))
  dtq=dtq1
  if(dtq>dtq2)dtq=dtq2

  if(dtq>dtc)dtq=dtc
  nstep_qm=nint(dtc/dtq)
  dtq=dtc/real(nstep_qm)
  hop_prob=0.d0
  hop_prob_net=0.d0
  V_k_hold=V_k
  V_k=V_k_old
  call compute_mat_adiab

  flag_hop=0
  do i=1,nstep_qm
    call compute_hop_prob(dtq)
    if(flag_hop==0)call check_hop(i*dtq)
    call rk4(ci,dtq,dVk_dt)
  enddo

  if(icollapse==1) then
    call verlet_decoherence(dtc,W_overlap,V_k_old,dvk_dt)
  endif

  do i=1,nquant
    if(hop_prob_net(i)<0.d0)hop_prob_net=0.d0
    hop_prob_net(i)=1.d0-dexp(-hop_prob_net(i))
  enddo

end subroutine evolve_quantum_small_dtq
!-----------------------------------------------------------------  

subroutine compute_hop_prob(dtq)
  implicit none
  real*8,intent(in)::dtq
  integer i
  real*8 pr

  do i=1,nquant
    if(i.ne.state) then
      pr=-2*real(ci(i)*dconjg(ci(state)))*vdotd(i,state)
      pr=pr*dtq/cdabs(ci(state))**2
      if(pr<0.d0)pr=0.d0     !!!! CAUTION AMBER CHECK !!!!
      hop_prob(i)=pr
      hop_prob_net(i)=hop_prob_net(i)+pr
    endif
  enddo

end subroutine compute_hop_prob
!-----------------------------------------------------------------  

subroutine check_hop(tim)
  implicit none
  real*8,intent(in)::tim
  integer i
  real*8 rnd,pr

  call random_number(rnd)
  pr=0.d0
  flag_hop=0
  do i=1,nquant
    if(i.ne.state) then
      pr=pr+hop_prob(i)
      if(rnd<pr) then
        state_tentative=i
        flag_hop=1
        exit
      endif
    endif
  enddo

end subroutine check_hop
!-----------------------------------------------------------------  

subroutine rk4(ci,dtq,dVk_dt)
  implicit none
  complex*16,intent(inout)::ci(nquant)
  real*8,intent(in) :: dtq,dVk_dt(nquant)
  complex*16,dimension(1:nquant):: k1,k2,k3,k4

  k1=matmul_lap(mat_adiab,ci)

  V_k=V_k+dVk_dt*dtq/2.d0
  call compute_mat_adiab

  k2=matmul_lap(mat_adiab,ci+0.5*dtq*k1)
  k3=matmul_lap(mat_adiab,ci+0.5*dtq*k2)

  V_k=V_k+dVk_dt*dtq/2.d0
  call compute_mat_adiab

  k4=matmul_lap(mat_adiab,ci+dtq*k3)

  ci=ci+dtq/6.d0*(k1+2*k2+2*k3+k4)

end subroutine rk4
!-----------------------------------------------------------------  

subroutine verlet_decoherence(dt,W_mat,V_k0,dvk_dt)
  implicit none
  real*8,intent(in):: dt,W_mat(nquant,nquant),V_k0(nquant),dvk_dt(nquant)
  real*8 acc_dec(nquant,nclass,3),delf(nquant,nclass,3),temp(nclass,3)
  complex*16 temp_delr(nquant,nclass,3),temp_delp(nquant,nclass,3)
  !complex*16 ci_diab(nquant)
  integer i,j,k

  delF=force_old
  temp=delF(state,:,:)
  do i=1,nquant
    do j=1,nclass
      delF(i,j,:)=delF(i,j,:)-temp(j,:)
      acc_dec(i,j,:)=delF(i,j,:)*cdabs(ci_old(i))**2/mass(j)
    enddo
  enddo

  do i=1,nquant
    do j=1,nclass
      delr(i,j,:)=delr(i,j,:)+delp(i,j,:)/mass*dt+0.5*acc_dec(i,j,:)*dt**2
      delp(i,j,:)=delp(i,j,:)+0.5*mass(j)*acc_dec(i,j,:)*dt
    enddo
  enddo

  !ci_diab=cdexp(iota*V_k0*dt/hbar)*cdexp(0.5*iota*dvk_dt*dt**2/hbar)*ci
  !ci_diab=matmul_lap(W_mat,ci_diab)
  delF=0.d0
  do j=1,nquant
    do k=1,nquant
      delF(j,:,:)=delF(j,:,:)+dabs(W_mat(j,k)**2)*(force(k,:,:)-force(state,:,:))
    enddo
  enddo
  !temp=delF(state,:)
  do i=1,nquant
    do j=1,nclass
      acc_dec(i,j,:)=delF(i,j,:)*cdabs(ci_old(i))**2/mass(j)
    enddo
  enddo

  do i=1,nquant
    do j=1,nclass
      delp(i,j,:)=delp(i,j,:)+0.5*mass(j)*acc_dec(i,j,:)*dt
    enddo
  enddo

  temp_delr=0.d0;temp_delp=0.d0
  do j=1,nquant
    do k=1,nquant
      temp_delr(j,:,:)=temp_delr(j,:,:)+dabs(W_mat(k,j)**2)*delr(k,:,:)
      temp_delp(j,:,:)=temp_delp(j,:,:)+dabs(W_mat(k,j)**2)*delp(k,:,:)
    enddo
  enddo
  delr=temp_delr
  delp=temp_delp

  !do i=1,nclass
  !  delr(:,i)=delr(:,i)-delr(state,i)
  !  delp(:,i)=delp(:,i)-delp(state,i)
  !enddo

end subroutine verlet_decoherence
!-----------------------------------------------------------------  

subroutine evolve_quantum_adiabatic
  implicit none
  complex*16,dimension(1:nquant):: k1,k2,k3,k4

  k1=matmul_lap(mat_adiab,ci)

  call evolve_classical(dtc/2.d0)
  call compute_mat_adiab

  k2=matmul_lap(mat_adiab,ci+0.5*dtc*k1)
  k3=matmul_lap(mat_adiab,ci+0.5*dtc*k2)

  call evolve_classical(dtc/2.d0)
  call compute_mat_adiab

  k4=matmul_lap(mat_adiab,ci+dtc*k3)

  ci=ci+dtc/6.d0*(k1+2*k2+2*k3+k4)

end subroutine evolve_quantum_adiabatic
!-----------------------------------------------------------------  

subroutine evolve_classical(dt)
  !! Velocity Verlet
  implicit none
  real*8,intent(in) :: dt
  real*8 gama_dt,c0,c1,c2
  real*8 delta_r(nclass,3),delta_v(nclass,3)
  real*8 t1,t2

  !call cpu_time(t1)

  if(ifriction==0) then
    !! Step 1
    x_old=x
    x=x+v*dt+0.5*acc*dt*dt
    v=v+0.5*acc*dt
    call rattle(1,dt)
    call pbc_adjust

    acc_old=acc
    call evaluate_variables(0)

    !! Step 2
    v=v+0.5*dt*acc
    call rattle(2,dt)

    call evaluate_variables(1)

  endif

  if(ifriction==1) then
    gama_dt=gamma_B*dt
    c0=dexp(-gama_dt)
    c1=1.d0/gama_dt*(1.d0-c0)
    c2=1.d0/gama_dt*(1.d0-c1)
   
     call stochastic_force(delta_r,delta_v,dt)
     x_old=x
     x=x+c1*dt*v+c2*dt*dt*acc+delta_r
     v=c0*v+(c1-c2)*dt*acc+delta_v
     call rattle(1,dt)
     call pbc_adjust

     acc_old=acc
     call evaluate_variables(0)

     !v=c0*v+(c1-c2)*dt*acc_old+c2*dt*acc+delta_v

     v=v+c2*dt*acc!+delta_v
     call rattle(2,dt)

     call evaluate_variables(1)
  
  endif

  !call cpu_time(t2);tim_ev_cl=tim_ev_cl+t2-t1

end subroutine evolve_classical
!-----------------------------------------------------------------  

subroutine rattle(istep,dt)
  !! Allen and Tildesley, page 96
  !! Constraints solved analytically for Me-Cl bond lengths
  implicit none
  integer,intent(in)::istep
  real*8,intent(in)::dt
  real*8 aa,bb,cc,lambda
  real*8 discr,tmp1,tmp2
  integer i
  real*8 x_Mecl(3),xp_MeCl(3),vp_MeCl(3)
  real*8 t1,t2

  !call cpu_time(t1)

  if(istep==1) then
    !! Step 1
    do i=3,nclass-1,2
      call pbc(x_old(i,:),x_old(i+1,:),x_MeCl,tmp1)
      call pbc(x(i,:),x(i+1,:),xp_MeCl,tmp1)

      aa=r_MeCl**2*dt**4/(4*mass_MeCl**2)
      bb=dt*dt*sum(x_MeCl*xp_MeCl)/mass_MeCl
      cc=sum(xp_MeCl*xp_MeCl)-r_MeCl**2
      discr=bb**2-4*aa*cc
      if(discr<0.d0) then
        write(6,*) "problem in rattle",discr
        stop
      endif
      discr=dsqrt(discr)
      tmp1=(-bb+discr)/(2*aa)
      tmp2=(-bb-discr)/(2*aa)
      if(dabs(tmp1)>dabs(tmp2))lambda=tmp2
      if(dabs(tmp1)<=dabs(tmp2))lambda=tmp1

      x(i,:)=x(i,:)-0.5d0*lambda*x_MeCl*dt*dt/mass(i)
      v(i,:)=v(i,:)-0.5d0*lambda*x_MeCl*dt/mass(i)
      x(i+1,:)=x(i+1,:)+0.5d0*lambda*x_MeCl*dt*dt/mass(i+1)
      v(i+1,:)=v(i+1,:)+0.5d0*lambda*x_MeCl*dt/mass(i+1)
    enddo
  endif

  if(istep==2) then
    !! Step 2
    do i=3,nclass-1,2
      !x_MeCl=x(i+1,:)-x(i,:)
      call pbc(x(i,:),x(i+1,:),x_MeCl,tmp1)
      vp_MeCl=v(i+1,:)-v(i,:)
      lambda=-2*sum(vp_MeCl*x_MeCl)*mass_MeCl/r_MeCl**2
      v(i,:)=v(i,:)-0.5d0*lambda*x_MeCl/mass(i)
      v(i+1,:)=v(i+1,:)+0.5d0*lambda*x_MeCl/mass(i+1)
    enddo
  endif

  !call cpu_time(t2);tim_rattle=tim_rattle+t2-t1

end subroutine rattle
!-----------------------------------------------------------------  

subroutine traj_terminate(iterm)
  implicit none
  integer,intent(out) :: iterm

  iterm=0

end subroutine traj_terminate
!-----------------------------------------------------------------  

subroutine compute_mat_diab
  implicit none
  integer i,j
  real*8 t1,t2

  !call cpu_time(t1)

  mat=0.d0
  do i=1,nbasis
    do j=1,nbasis
      mat(i,j)=-iota/hbar*sum(si_adiab(i,:)*si_adiab(j,:)*V_k(1:nquant))
    enddo
  enddo

  !call cpu_time(t2)
  !tim_mat=tim_mat+t2-t1

end subroutine compute_mat_diab
!-----------------------------------------------------------------  

subroutine compute_mat_adiab
  implicit none
  integer i,j
  real*8 t1,t2
  real*8 V_avg
  
  !call cpu_time(t1)

  mat_adiab=-vdotd
  V_avg=sum(V_k)/real(nquant)
  do i=1,nquant
    mat_adiab(i,i)=mat_adiab(i,i)-iota/hbar*(V_k(i)-V_avg)
  enddo
      
  !call cpu_time(t2)
  !tim_mat=tim_mat+t2-t1
  
end subroutine compute_mat_adiab
!-----------------------------------------------------------------  

subroutine hop
  implicit none
  integer ifrust

  if(flag_hop==1) then
    call velocity_adjust(state_tentative,ifrust)
  endif

end subroutine hop
!-----------------------------------------------------------------  

subroutine velocity_adjust(state_tentative,ifrust)
  implicit none
  integer,intent(in)::state_tentative
  integer,intent(out)::ifrust
  real*8 gij,gama,aa,bb,cc,discr,rij(3),dp(nclass,3),dij,dd(nclass,3),f1,f2,vd
  integer i,j,k,kp

  k=state;kp=state_tentative
  cc=V_k(k)-V_k(kp)
  call compute_dij_2state(k,kp,dd)

  aa=0.d0
  bb=0.d0
  do i=1,nclass
    if(i==1.or.i==2) then
      gij=0.d0
    else if(mod(i,2)==1) then
      !rij=x(i,:)-x(i+1,:)
      call pbc(x(i+1,:),x(i,:),rij,dij)
      gij=sum(rij*(dd(i,:)/mass(i)-dd(i+1,:)/mass(i+1)))
      gij=gij/(sum(rij*rij)*(1.d0/mass(i)+1.d0/mass(i+1)))
    else
      !rij=x(i,:)-x(i-1,:)
      call pbc(x(i-1,:),x(i,:),rij,dij)
      gij=sum(rij*(dd(i,:)/mass(i)-dd(i-1,:)/mass(i-1)))
      gij=gij/(sum(rij*rij)*(1.d0/mass(i)+1.d0/mass(i-1)))
    endif

    dp(i,:)=dd(i,:)-gij*rij

    aa=aa+0.5/mass(i)*sum(dp(i,:)*dp(i,:))
    bb=bb+sum(v(i,:)*dp(i,:))

  enddo

  discr=bb**2+4*aa*cc
  if(discr<0.d0) then
    ifrust=1
    cnt_frust=cnt_frust+1.d0
    if(flag_frust==0)then
      gama=0.d0
      call compute_delH_dels_ad
      f1=sum(force(k,:,:)*dp)
      f2=sum(force(kp,:,:)*dp)
      vd=sum(v*dp)
      if(f1*f2<0.d0.and.vd*f2<0.d0) then
      !if(f1*f2<0.d0) then
        gama=bb/aa
      endif
    endif
    if(flag_frust>0)gama=0.d0
  else
    ifrust=0
    if(bb>=0.d0) gama=(bb-dsqrt(discr))/(2*aa)
    if(bb<0.d0)  gama=(bb+dsqrt(discr))/(2*aa)
    state=state_tentative
    delr=0.d0
    delp=0.d0
  endif

  do i=1,nclass
    v(i,:)=v(i,:)-gama*dp(i,:)/mass(i)
  enddo

  call evaluate_variables(0)
  call evaluate_variables(1)

end subroutine velocity_adjust
!-----------------------------------------------------------------  

subroutine collapse(dt,iflag_coll)
  implicit none
  real*8,intent(in) :: dt
  integer,intent(out) :: iflag_coll
  real*8 rnd,gama_collapse,gama_reset
  complex*16 su1
  real*8 su2
  integer n,i,j

  i=state

  if(icollapse==1) then

    iflag_coll=0
    do n=1,nquant
      if(n.ne.state) then
        gama_reset=0.d0
        do j=1,3
          gama_reset=gama_reset+sum((force(n,:,j)-force(i,:,j))*dble(delr(n,:,j)-delr(state,:,j)))/(2*hbar)
        enddo
        !! CAUTION !! !! Assumes delr(n,n,:) is in direction of v(:) !!
        su1=0.d0;su2=0.d0
        do j=1,3
          su1=su1+cdabs((V_k(i)-V_k(n))*vdotd(i,n)*sum((delr(n,:,j)-delr(state,:,j))*v(:,j)))
          su2=su2+sum(v(:,j)*v(:,j))
        enddo
        su1=su1/su2
        gama_collapse=gama_reset-2/hbar*cdabs(su1)
        gama_collapse=gama_collapse*dt
        gama_reset=-gama_reset*dt
        call random_number(rnd)

        if(rnd<gama_collapse) then
          iflag_coll=1
          cnt_collapse=cnt_collapse+1
          if(icollapse==1) then
            !do j=1,nquant
            !  if(j.ne.n) ci(j)=ci(j)/dsqrt(1-cdabs(ci(n)**2))
            !enddo
            !! Erratum: Landry, Subotnik JCP 137, 229901 (2012)
            ci(i)=ci(i)/cdabs(ci(i))*dsqrt(cdabs(ci(i))**2+cdabs(ci(n))**2)
            ci(n)=0.d0

          endif
        endif
        if(rnd<gama_collapse.or.rnd<gama_reset) then
          if(icollapse==1) then
            delr(n,:,:)=0.d0
            delp(n,:,:)=0.d0
          endif
        endif
      endif
    enddo

  endif

end subroutine collapse
!-----------------------------------------------------------------  

subroutine write_output(n,nflag)
  implicit none
  integer,intent(in)::nflag,n
  integer i
  real*8 r12,x12(3),v12(3)
  real*8 t1,t2
  real*8 phase

  !call cpu_time(t1)

  if(nflag==0) then
    if(iwrite==1) then
      if(mod(n,nstep_write)==1.or.nstep_write==1) then
        call pbc(x(3,:),x(4,:),x12,r12)
        v12=v(4,:)-v(3,:)
        write(10,'(6es17.7,i5)')curr_time*1.d15,energy/wave_to_J,sum(cdabs(ci)**2),r12*1.d10,sum(v12*x12),temperature,state
        write(11,'(f15.5$)')curr_time*1.d15
        write(12,'(6f15.5)')curr_time*1.d15,cdabs(ci(1:2))**2,r_exp*1.d10,V_k/wave_to_J!,cdabs(matmul_lap(si_adiab,ci))**2
        !write(13,'(9es15.5)')curr_time*1.d15,1.d10*delr(1,1,1),delr(2,2,1)*1.d10,delp(1,1,1)/mass,delp(2,2,1)/mass
        do i=1,1!nclass
          write(11,'(7f15.5$)')x(i,:)*1.d10,v(i,:),rAB*1.d10
        enddo
        write(11,*)!;write(12,*)
      endif
    endif
  endif

  if(nflag==1) then
    if(iwrite==0)write(10,'(5es15.5)')traj_num,energy/wave_to_J,sum(cdabs(ci)**2),r_exp*1.d10,temperature
    if(iwrite==1) then
      write(10,*)"traj num=",traj_num
      write(10,*)"standard deviation=",dsqrt((ensq_avg-en_avg**2/dfloat(nst_av))/dfloat(nst_av))/wave_to_J
      write(10,*)"ci**2=",sum(cdabs(ci)**2)
      write(10,*);write(10,*)
      write(11,*);write(11,*)
      write(12,*);write(12,*)
      write(13,*);write(13,*)
    endif
  endif

  !call cpu_time(t2)
  !tim_wr_out=tim_wr_out+t2-t1

end subroutine write_output
!-----------------------------------------------------------------  

subroutine write_average
  implicit none
  integer i,j
  real*8 nf

  nf=dfloat(n_traj)
  cnt_frust=cnt_frust/nf
  cnt_collapse=cnt_collapse/nf

  !pop=pop/nf
  !pop_surf=pop_surf/nf
  !pop_amp=pop_amp/nf

  do i=1,nsteps/nstep_avg+1
    write(100,'(6f15.7)')(i-1)*nstep_avg*dtc*1.d15,pop(1,i),cnt_init
  enddo

  write(101,*) cnt_frust,cnt_collapse

end subroutine write_average
!-----------------------------------------------------------------  

subroutine evaluate_variables(flag)
  implicit none
  integer,intent(in):: flag
  integer i,j,ndof

  if(flag==0) then
    !call compute_distances
    call tise(flag_classical)
    r_exp=sum(si_adiab(:,state)**2*rAH_grid)
  endif

  if(flag==1) then
    KE_en=0.d0
    do i=1,nclass
      KE_en=KE_en+0.5*mass(i)*sum(v(i,:)*v(i,:))
    enddo

    ndof=(3*nclass-6)-(nclass/2-1)
    temperature=KE_en/(ndof*kb)
    energy=pot_en+KE_en
   
  endif

end subroutine evaluate_variables
!-----------------------------------------------------------------  

subroutine tise(flag_classical)
  !! Output - pot_en,acc
  !! Output - V_k,d_ij
  implicit none
  integer,intent(in)::flag_classical
  integer i,j,k
  real*8 Hamil(nbasis,nbasis),ens(nbasis),vect(nbasis,nquant)
  real*8 pot_cl,acc_cl(nclass,3),acc_qm(nclass,3),dpotcl_dx(nclass,3)
  real*8 si_adiab_old(nquant,nbasis)
  real*8 t1,t2

  !call cpu_time(t1)

  call compute_potential(Hamil,delH_dels)
  Hamil_diab=Hamil
  call diag(Hamil,nbasis,ens,vect,nquant)

!  si_adiab_old=si_adiab
  do i=1,nquant
    si_adiab(:,i)=vect(:,i)
    if(sum(si_adiab(:,i)*si_adiab_prev(:,i))<0.d0)si_adiab(:,i)=-si_adiab(:,i)
  enddo

  do i=1,nclass
    do j=1,3
      delH_dels_ad(state,state,i,j)=sum(si_adiab(:,state)*matmul_lap(delH_dels(:,:,i,j),si_adiab(:,state)))
    enddo
  enddo
  !call cpu_time(t2);tim_diag=tim_diag+(t2-t1)

  !call cpu_time(t1)

  if(flag_classical==1) call potential_classical(pot_cl,dpotcl_dx)
  if(flag_classical==0) then
    pot_cl=0.d0
    dpotcl_dx=0.d0
  endif

  do i=1,nclass
    acc_qm(i,:)=-1.d0/mass(i)*delH_dels_ad(state,state,i,:)
    acc_cl(i,:)=-1.d0/mass(i)*dpotcl_dx(i,:)
  enddo
  pot_en=pot_cl+ens(state)
  acc=acc_cl+acc_qm

  V_k=pot_cl+ens(1:nquant)

  !call cpu_time(t2);tim_cl=tim_cl+(t2-t1)

end subroutine tise
!-----------------------------------------------------------------  

subroutine compute_delH_dels_ad
  implicit none
  integer i,j,k

  force=0.d0
  do k=1,nquant
    do i=1,nclass
      do j=1,3
        delH_dels_ad(k,k,i,j)=sum(si_adiab(:,k)*matmul_lap(delH_dels(:,:,i,j),si_adiab(:,k)))
      enddo
    enddo
    force(k,:,:)=-delH_dels_ad(k,k,:,:)
  enddo

end subroutine compute_delH_dels_ad
!-----------------------------------------------------------------  

subroutine compute_dij
  implicit none
  integer i,j,k,kp

  do k=1,nquant-1
    do kp=k+1,nquant
      do i=1,nclass
        do j=1,3
          d_ij(k,kp,i,j)=sum(si_adiab(:,k)*matmul_lap(delH_dels(:,:,i,j),si_adiab(:,kp)))
        enddo
      enddo
      d_ij(k,kp,:,:)=d_ij(k,kp,:,:)/(V_k(kp)-V_k(k))
      d_ij(kp,k,:,:)=-d_ij(k,kp,:,:)
    enddo
  enddo

end subroutine compute_dij
!-----------------------------------------------------------------  

subroutine compute_dij_2state(k,kp,dp)
  implicit none
  integer,intent(in):: k,kp
  real*8,intent(out):: dp(nclass,3)
  real*8 x_sav(nclass)
  integer i,j

  do i=1,nclass
    do j=1,3
      dp(i,j)=sum(si_adiab(:,k)*matmul_lap(delH_dels(:,:,i,j),si_adiab(:,kp)))
    enddo
  enddo
  dp=dp/(V_k(kp)-V_k(k))

end subroutine compute_dij_2state
!-----------------------------------------------------------------  

subroutine compute_vdotd
  ! Meek, Levine, JPCL 5, 2351 (2014). Look at Supp info.
  implicit none
  integer i,j,k
  real*8,dimension(nquant,nquant) :: W,ci_W,si_W
  real*8 A,B,C,D,E
  real*8 Wlj,Wlk

  !Method 1
  !call compute_dij
  !vdotd=0.d0
  !do i=1,nclass
  !  do j=1,3
  !    vdotd=vdotd+v(i,j)*d_ij(:,:,i,j)
  !  enddo
  !enddo

  !Method 2
!  do j=1,nquant
!    do k=1,nquant
!      W(j,k)=sum(si_adiab_prev(:,j)*si_adiab(:,k))
!      ci_W(j,k)=dacos(W(j,k))
!      si_W(j,k)=dasin(W(j,k))
!    enddo
!  enddo
!
!  vdotd=0.d0
!  do k=1,nquant-1
!    do j=k+1,nquant
!      A=-sinx_x(ci_W(j,j)-si_W(j,k))
!      B=sinx_x(ci_W(j,j)+si_W(j,k))
!      C=sinx_x(ci_W(k,k)-si_W(k,j))
!      D=sinx_x(ci_W(k,k)+si_W(k,j))
!      Wlj=dsqrt(1.d0-W(j,j)**2-W(k,j)**2)
!      if(Wlj==0.d0.or.nquant==2) then
!        E=0.d0
!      else
!        Wlk=(-W(j,k)*W(j,j)-W(k,k)*W(k,j))/Wlj
!        E=2*dasin(Wlj)/(dasin(Wlj)**2-dasin(Wlk)**2)
!        E=E*(Wlj*Wlk*dasin(Wlj)+dasin(Wlk)*(dsqrt((1-Wlj**2)*(1-Wlk**2))-1.d0))
!      endif
!      vdotd(k,j)=0.5/dtc*(ci_W(j,j)*(A+B)+si_W(k,j)*(C+D)+E)
!      vdotd(j,k)=-vdotd(k,j)
!    enddo
!  enddo

  !Method 3
  do i=1,nquant
    do j=1,nquant
      W_overlap(i,j)=sum(si_adiab_prev(:,i)*si_adiab(:,j))
    enddo
  enddo

  call orthoganalize(W_overlap,nquant)
  call logm(W_overlap,vdotd,nquant)
  vdotd=vdotd/dtc

end subroutine compute_vdotd
!-----------------------------------------------------------------  

subroutine orthoganalize(mat,n)
  integer,intent(in)::n
  real*8,intent(inout)::mat(n,n)
  real*8 S_mat(n,n)

  S_mat=matmul_lap(transpose(mat),mat)
  call inverse_squareroot(S_mat,n)
  mat=matmul_lap(mat,S_mat)

end subroutine orthoganalize
!-----------------------------------------------------------------  

subroutine setup_parameters
  implicit none
  integer i
  real*8 Aij,Cij

  !-----------------------------------------------------------------  
  rAH_del=(rAH_max-rAH_min)/real(nbasis-1)
  mass_H=1.d0*amu2kg
  do i=1,nbasis
    rAH_grid(i)=rAH_min+(i-1)*rAH_del
  enddo
  call compute_KE_matrix_dvr(KE_DVR,nbasis,rAH_del,mass_H)
  !-----------------------------------------------------------------  

  A_pot(1)=dsqrt(7.95d6*kcal_to_J*1.d-120)
  A_pot(2)=dsqrt(5.25d6*kcal_to_J*1.d-120)
  C_pot(1)=dsqrt(2750.d0*kcal_to_J*1.d-60)
  C_pot(2)=dsqrt(2950.d0*kcal_to_J*1.d-60)

  r_MeCl=1.781d-10

  a_solute=11.2d10
  b_solute=7.1d13*kcal_to_J
  c_solute=0.776
  Diss_A=110.d0*kcal_to_J
  d_A=0.95d-10
  d_b=0.97d-10
  n_A=9.26d10
  n_B=11.42d10

  r0=1.43d-10;l0=0.125d-10
  e_alpha_C(1)=-0.5d0
  e_alpha_C(2)=0.5d0
  e_alpha_C(3)=0.d0
  e_alpha_i(1)=-1.d0
  e_alpha_i(2)=0.5d0
  e_alpha_i(3)=0.5d0

  box_length=28.d-10

  r_c=13.8d-10
  r_t=0.95d0*r_c

  mass(1)=93.d0*amu2kg
  mass(2)=59.d0*amu2kg
  atom(1)='HB'
  atom(2)='HB'

  sig_sol_com=3.5d-10;eps_sol_com=200.d0*kb

  Aij=A_pot(1)*A_pot(1)
  Cij=C_pot(1)*C_pot(1)
  sig_meth_meth=(Aij/Cij)**(1/6.d0)
  eps_meth_meth=Cij**2/(4*Aij)

  Aij=A_pot(1)*A_pot(2)
  Cij=C_pot(1)*C_pot(2)
  sig_meth_cl=(Aij/Cij)**(1/6.d0)
  eps_meth_cl=Cij**2/(4*Aij)

  Aij=A_pot(2)*A_pot(2)
  Cij=C_pot(2)*C_pot(2)
  sig_cl_cl=(Aij/Cij)**(1/6.d0)
  eps_cl_cl=Cij**2/(4*Aij)

  do i=3,nclass,2
    mass(i)=15.d0*amu2kg
    charge(i)=0.25d0
    atom(i)='Me'
  enddo
  do i=4,nclass,2
    mass(i)=35.5d0*amu2kg
    charge(i)=-0.25d0
    atom(i)='Cl'
  enddo

  mass_MeCl=1.d0/(1.d0/mass(3)+1.d0/mass(4))

  call LJ(LJ_cut_meth_meth,dLJ_cut_meth_meth,R_c,eps_meth_meth,sig_meth_meth)
  call LJ(LJ_cut_meth_cl,dLJ_cut_meth_cl,R_c,eps_meth_cl,sig_meth_cl)
  call LJ(LJ_cut_cl_cl,dLJ_cut_cl_cl,R_c,eps_cl_cl,sig_cl_cl)
  call LJ(LJ_cut_sol_com,dLJ_cut_sol_com,R_c,eps_sol_com,sig_sol_com)

  Rc_verlet=14.d-10*dsqrt(3.d0)
  r_cutoff=1.25d-10

  gamma_B=1000.d0*2*pi*clight

end subroutine setup_parameters
!-----------------------------------------------------------------  

subroutine potential_classical(pot_cl,dv_dx)
  implicit none
  real*8,intent(out) :: pot_cl,dv_dx(nclass,3)
  real*8 pot_LJ,dv_dx_LJ(nclass,3)
  real*8 pot_coul,dv_dx_coul(nclass,3)
  real*8 xij(3),rij
  integer i,j,k
  integer c1,c2,d1,d2

  pot_LJ=0.d0;pot_coul=0.d0
  dv_dx_LJ=0.d0;dv_dx_coul=0.d0

!  call system_clock(d1)
  do i=1,nclass-1
    do j=i+1,nclass
!    if(verlet_list(i,1)>0) then
!      do k=2,verlet_list(i,1)+1
!        j=verlet_list(i,k)
        call pbc(x(i,:),x(j,:),xij,rij)
        if(rij<r_C) then
          call LJ_total(pot_LJ,dv_dx_LJ,xij,rij,i,j)
        endif
        !if(i>2)call coul_solv_solv(pot_coul,dv_dx_coul,xij,rij,i,j)
      enddo
!    endif
  enddo
!  call system_clock(d2)
!  tim_LJ=tim_LJ+(d2-d1)/real(cnt_rate)

!  call system_clock(c1)
  do i=3,nclass-2,2
    do j=i+2,nclass,2
      call coul_molecular_solvent(pot_coul,dv_dx_coul,i,j)
    enddo
  enddo
!  call system_clock(c2)
!  tim_cl=tim_cl+(c2-c1)/real(cnt_rate)


  pot_cl=pot_LJ+pot_coul
  dv_dx=dv_dx_LJ+dv_dx_coul


end subroutine potential_classical
!-----------------------------------------------------------------  

subroutine pot_AHB(pot,dpot_dx,rAH)
  implicit none
  real*8,intent(out):: pot,dpot_dx(nclass,3)
  real*8,intent(in):: rAH
  real*8 morse1,dM1_dr
  real*8 morse2,dM2_dr
  real*8 dV_drAB,dV_drAH,x_H(3),x12(3),xAH(3)
  real*8 tmp

  call pbc(x(1,:),x(2,:),x12,rAB)

  x_H=x(1,:)+rAH*x12/rAB
  call pbc(x(1,:),x_H,xAH,tmp)
  call morse(morse1,dM1_dr,rAH,Diss_A,n_A,d_A)
  call morse(morse2,dM2_dr,(RAB-rAH),c_solute*Diss_A,n_B,d_B)
  tmp=b_solute*dexp(-a_solute*rAB)
  dV_drAB=tmp*(-a_solute)+dM2_dr
  dV_drAH=dM1_dr-dM2_dr

  pot=morse1+morse2+tmp
  dpot_dx=0.d0
  dpot_dx(1,:)=dv_drAB*(-x12)/rAB!+dv_drAH*(-xAH)/rAH
  dpot_dx(2,:)=dv_drAB*(x12)/rAB

!!! CAREFUL
!!! pot_AHB(|r_AB|,r_AH), while treating r_AH, and r_A(3) INDEPENDENTLY!!
!!! Hence dV/dr_A(3) = dV/d|r_AB| d|rAB|/dr_A(3)

end subroutine pot_AHB
!-----------------------------------------------------------------

subroutine LJ_total(pot,dv_dx,xij,rij,i,j)
  implicit none
  real*8,intent(inout)::pot,dv_dx(nclass,3)
  integer,intent(in):: i,j
  real*8,intent(in):: rij,xij(3)
  real*8 LJ_ij,dpot_drij,vec(3)
  real*8 t1,t2

  !call cpu_time(t1)

!  pot=0.d0
!  dv_dx=0.d0
!  do i=1,nclass-1
!    do j=i+1,nclass
!      call pbc(x(i,:),x(j,:),xij,rij)

!      if(rij<R_c) then
        vec=0.d0
        if(atom(i).eq.'Me'.and.atom(j).eq.'Me') then
          call LJ(LJ_ij,dpot_drij,rij,eps_meth_meth,sig_meth_meth)
          pot=pot+LJ_ij-LJ_cut_meth_meth-dLJ_cut_meth_meth*(rij-R_c)
          vec(:)=(dpot_drij-dLJ_cut_meth_meth)*(-xij)/rij
        endif

        if(atom(i).eq.'Me'.and.atom(j).eq.'Cl'.and.j>i+1) then
          call LJ(LJ_ij,dpot_drij,rij,eps_meth_cl,sig_meth_cl)
          pot=pot+LJ_ij-LJ_cut_meth_cl-dLJ_cut_meth_cl*(rij-R_c)
          vec(:)=(dpot_drij-dLJ_cut_meth_cl)*(-xij)/rij
        endif

        if(atom(i).eq.'Cl'.and.atom(j).eq.'Me') then
          call LJ(LJ_ij,dpot_drij,rij,eps_meth_cl,sig_meth_cl)
          pot=pot+LJ_ij-LJ_cut_meth_cl-dLJ_cut_meth_cl*(rij-R_c)
          vec(:)=(dpot_drij-dLJ_cut_meth_cl)*(-xij)/rij
        endif

        if(atom(i).eq.'Cl'.and.atom(j).eq.'Cl') then
          call LJ(LJ_ij,dpot_drij,rij,eps_cl_cl,sig_cl_cl)
          pot=pot+LJ_ij-LJ_cut_cl_cl-dLJ_cut_cl_cl*(rij-R_c)
          vec(:)=(dpot_drij-dLJ_cut_cl_cl)*(-xij)/rij
        endif

        if(atom(i).eq.'HB'.and.(atom(j).eq.'Me'.or.atom(j).eq.'Cl')) then
          call LJ(LJ_ij,dpot_drij,rij,eps_sol_com,sig_sol_com)
          pot=pot+LJ_ij-LJ_cut_sol_com-dLJ_cut_sol_com*(rij-R_c)
          vec(:)=(dpot_drij-dLJ_cut_sol_com)*(-xij)/rij
        endif

        dv_dx(i,:)=dv_dx(i,:)+vec
        dv_dx(j,:)=dv_dx(j,:)-vec

!      endif
!    enddo
!  enddo

  !call cpu_time(t2)
  !tim_LJ_tot=tim_LJ_tot+t2-t1

end subroutine LJ_total
!-----------------------------------------------------------------  

subroutine coul_solv_solv(pot,dv_dx,xij,rij,i,j)
  implicit none
  real*8,intent(inout)::pot,dv_dx(nclass,3)
  integer,intent(in):: i,j
  real*8,intent(in)::rij,xij(3)
  real*8 pot_ij,dpot_drij
  real*8 T,dT_drij,vec(3)
  real*8 t1,t2

  !call cpu_time(t1)

  !pot=0.d0
  !dv_dx=0.d0
  !do i=3,nclass-1
  !  do j=i+1,nclass
      if(mod(i,2)==0.or.j>i+1) then
        !rij=distance(i,j)
        !call pbc(x(i,:),x(j,:),xij,rij)
        call coul(pot_ij,dpot_drij,rij,charge(i),charge(j))
        call steinhauser(T,dT_drij,rij,r_t,r_c)
        pot=pot+pot_ij*T
        vec=(dpot_drij*T+pot_ij*dT_drij)*(-xij)/rij
        dv_dx(i,:)=dv_dx(i,:)+vec
        dv_dx(j,:)=dv_dx(j,:)-vec
      endif
  !  enddo
  !enddo

  !call cpu_time(t2)
  !tim_solv_solv=tim_solv_solv+t2-t1

end subroutine coul_solv_solv
!-----------------------------------------------------------------  

subroutine coul_solv_solute(pot,dv_dx,xij,rij,i,j,k)
  implicit none
  real*8,intent(inout)::pot(nbasis),dv_dx(nbasis,nclass,3)
  real*8,intent(in)::xij(3),rij
  integer,intent(in):: i,j,k
  real*8 pot_ij,dpot_drij
  real*8 T,dT_drij,vec(3)

  !do j=3,nclass
  !  do i=1,2
      !rij=distance(i,j)
      !call pbc(x(i,:),x(j,:),xij,rij)
      call coul(pot_ij,dpot_drij,rij,charge(i),charge(j))
      call steinhauser(T,dT_drij,rij,r_t,r_c)
      pot(k)=pot(k)+pot_ij*T
      vec=(dpot_drij*T+pot_ij*dT_drij)*(-xij)/rij
      dv_dx(k,i,:)=dv_dx(k,i,:)+vec
      dv_dx(k,j,:)=dv_dx(k,j,:)-vec
  !  enddo

end subroutine coul_solv_solute
!-----------------------------------------------------------------  

subroutine coul_molecular(pot,dv_dx,i,j,k)
  implicit none
  real*8,intent(inout)::pot(nbasis),dv_dx(nbasis,nclass,3)
  integer,intent(in):: i,j,k
  real*8 xij(3),rij
  real*8 xii1(3),xjj1(3)
  real*8 xij_st(3),rij_st
  integer k1,k2
  real*8 pot_ij,dpot_drij,pot_tot
  real*8 T,dT_drij,vec(3),drij_dx(nclass,3)

  call center(x(i,:),x(i+1,:),xii1)
  call center(x(j,:),x(j+1,:),xjj1)
  call pbc(xii1,xjj1,xij_st,rij_st)
  call steinhauser(T,dT_drij,rij_st,r_t,r_c)

  if(rij_st<r_c) then
    drij_dx=0.d0
    drij_dx(i,:)=-0.5*xij_st/rij_st
    drij_dx(i+1,:)=-0.5*xij_st/rij_st
    drij_dx(j,:)=0.5*xij_st/rij_st
    drij_dx(j+1,:)=0.5*xij_st/rij_st

    pot_tot=0.d0
    do k1=i,i+1
      do k2=j,j+1
        call pbc(x(k1,:),x(k2,:),xij,rij)
        call coul(pot_ij,dpot_drij,rij,charge(k1),charge(k2))
        pot(k)=pot(k)+pot_ij*T
        pot_tot=pot_tot+pot_ij
        vec=dpot_drij*(-xij)/rij*T
        dv_dx(k,k1,:)=dv_dx(k,k1,:)+vec
        dv_dx(k,k2,:)=dv_dx(k,k2,:)-vec
      enddo
    enddo

    vec=pot_tot*dT_drij*(-0.5*xij_st/rij_st)
    dv_dx(k,i,:)   = dv_dx(k,i,:)   + vec
    dv_dx(k,i+1,:) = dv_dx(k,i+1,:) + vec
    dv_dx(k,j,:)   = dv_dx(k,j,:)   - vec
    dv_dx(k,j+1,:) = dv_dx(k,j+1,:) - vec
  endif

end subroutine coul_molecular
!-----------------------------------------------------------------  

subroutine coul_molecular_solvent(pot,dv_dx,i,j)
  implicit none
  real*8,intent(inout)::pot,dv_dx(nclass,3)
  integer,intent(in):: i,j
  real*8 xij(3),rij
  real*8 xii1(3),xjj1(3)
  real*8 xij_st(3),rij_st
  integer k1,k2
  real*8 pot_ij,dpot_drij,pot_tot
  real*8 T,dT_drij,vec(3),drij_dx(nclass,3)

  call center(x(i,:),x(i+1,:),xii1)
  call center(x(j,:),x(j+1,:),xjj1)
  call pbc(xii1,xjj1,xij_st,rij_st)
  call steinhauser(T,dT_drij,rij_st,r_t,r_c)

  if(rij_st<r_c) then
    drij_dx=0.d0
    drij_dx(i,:)=-0.5*xij_st/rij_st
    drij_dx(i+1,:)=-0.5*xij_st/rij_st
    drij_dx(j,:)=0.5*xij_st/rij_st
    drij_dx(j+1,:)=0.5*xij_st/rij_st

    pot_tot=0.d0
    do k1=i,i+1
      do k2=j,j+1
        call pbc(x(k1,:),x(k2,:),xij,rij)
        call coul(pot_ij,dpot_drij,rij,charge(k1),charge(k2))
        pot=pot+pot_ij*T
        pot_tot=pot_tot+pot_ij
        vec=dpot_drij*(-xij)/rij*T
        dv_dx(k1,:)=dv_dx(k1,:)+vec
        dv_dx(k2,:)=dv_dx(k2,:)-vec
      enddo
    enddo

    vec=pot_tot*dT_drij*(-0.5*xij_st/rij_st)
    dv_dx(i,:)   = dv_dx(i,:)   + vec
    dv_dx(i+1,:) = dv_dx(i+1,:) + vec
    dv_dx(j,:)   = dv_dx(j,:)   - vec
    dv_dx(j+1,:) = dv_dx(j+1,:) - vec
  endif

end subroutine coul_molecular_solvent
!-----------------------------------------------------------------  

subroutine coul_H_solute(pot,dv_dx,rAH,j,x_AB,k)
  implicit none
  real*8,intent(inout)::pot(nbasis),dv_dx(nbasis,nclass,3)
  real*8,intent(in)::rAH,x_AB(3)
  integer,intent(in):: j,k
  real*8 pot_ij,dpot_drij,dpot_drhx,pot_tot
  real*8 T,dT_drij,vec(3),X_H(3),xij(3),rij
  real*8 xii1(3),xjj1(3)
  real*8 xij_st(3),rij_st,dv_drhx
  integer k2

  !call pbc(x(1,:),x(2,:),xij,rij)

  x_H=x(1,:)+rAH*(x_AB)/rAB

  call center(x(j,:),x(j+1,:),xjj1)
  call pbc(x_H,xjj1,xij_st,rij_st)
  call steinhauser(T,dT_drij,rij_st,r_t,r_c)

  if(rij_st<r_c) then

    pot_tot=0.d0
    dv_drhx=0.d0
    do k2=j,j+1
      call pbc(x_H,x(k2,:),xij,rij)
      call coul(pot_ij,dpot_drij,rij,charge_H,charge(k2))
      pot(k)=pot(k)+pot_ij*T
      pot_tot=pot_tot+pot_ij
      vec=dpot_drij*(-xij)/rij*T
      dv_dx(k,k2,:)=dv_dx(k,k2,:)-vec
      dv_drhx=dpot_drij*T
      vec=xij/rij
      dv_dx(k,1,:)=dv_dx(k,1,:)-dv_drhx*vec*(1.d0-rAH/rAB)+dv_drhx*x_AB*(rAH/(rAB**3*rij))*sum(-xij*x_AB)
      dv_dx(k,2,:)=dv_dx(k,2,:)-dv_drhx*vec*(rAH/rAB)-dv_drhx*x_AB*(rAH/(rAB**3*rij))*sum(-xij*x_AB)
    enddo

    vec=pot_tot*dT_drij*(-0.5*xij_st/rij_st)
    dv_dx(k,j,:)   = dv_dx(k,j,:)   - vec
    dv_dx(k,j+1,:) = dv_dx(k,j+1,:) - vec
    dv_drhx=pot_tot*dT_drij
    vec=xij_st/rij_st
    dv_dx(k,1,:)=dv_dx(k,1,:)-dv_drhx*vec*(1.d0-rAH/rAB)+dv_drhx*x_AB*(rAH/(rAB**3*rij_st))*sum(-xij_st*x_AB)
    dv_dx(k,2,:)=dv_dx(k,2,:)-dv_drhx*vec*(rAH/rAB)-dv_drhx*x_AB*(rAH/(rAB**3*rij_st))*sum(-xij_st*x_AB)
  endif

end subroutine coul_H_solute
!-----------------------------------------------------------------  

subroutine morse(pot,dpot_dr,r,Diss_A,n_A,d_A)
  implicit none
  real*8,intent(out)::pot,dpot_dr
  real*8,intent(in)::r,diss_A,n_A,d_A
  real*8 tmp

  tmp=Diss_A*dexp(-n_A*(r-d_A)**2/(2*r))

  pot=Diss_A-tmp
  dpot_dr=tmp*n_A*0.5*(1-(d_A/r)**2)

end subroutine morse
!-----------------------------------------------------------------  

subroutine LJ(pot,dpot_drij,rij,eps,sig)
  implicit none
  real*8,intent(in)::rij,eps,sig
  real*8,intent(out)::pot,dpot_drij
  real*8 rij_inv,sig_rij_6

  rij_inv=1.d0/rij
  sig_rij_6=(sig*rij_inv)**6

  pot=4*eps*sig_rij_6*(sig_rij_6-1.d0)
  dpot_drij=-24*eps*rij_inv*sig_rij_6*(2*sig_rij_6-1.d0)

end subroutine LJ
!-----------------------------------------------------------------  

subroutine coul(pot,dpot_drij,rij,q1,q2)
  implicit none
  real*8,intent(out)::pot,dpot_drij
  real*8,intent(in)::rij,q1,q2

  pot=q1*q2*q_el/rij
  dpot_drij=-pot/rij

end subroutine coul
!-----------------------------------------------------------------  

subroutine steinhauser(T,dt_drij,rij,rt,rc)
  implicit none
  real*8,intent(out)::T,dt_drij
  real*8,intent(in)::rij,rt,rc

  if(rij<=rt) then
    T=1.d0;dT_drij=0.d0
  endif

  if(rij>=rc) then
    T=0.d0;dT_drij=0.d0
  endif

  if(rij>rt.and.rij<rc) then
    T=1-(rij-rt)**2*(3*rc-rt-2*rij)/(rc-rt)**3
    dt_drij=6*(rij-rt)*(rij-rc)/(rc-rt)**3
  endif

end subroutine steinhauser
!-----------------------------------------------------------------  

subroutine compute_charges(rAH)
  implicit none
  real*8,intent(in)::rAH
  real*8 f,r

  f(r)=0.5d0*(1+(r-r0)/dsqrt((r-r0)**2+l0**2))

  charge(1)=e_alpha_c(1)+f(rAH)*(e_alpha_i(1)-e_alpha_c(1))
  charge(2)=e_alpha_c(3)+f(rAH)*(e_alpha_i(3)-e_alpha_c(3))
  charge_H=e_alpha_c(2)+f(rAH)*(e_alpha_i(2)-e_alpha_c(2))

end subroutine compute_charges
!-----------------------------------------------------------------  

subroutine compute_potential(H_diab,delV_dels)
  implicit none
  real*8,intent(out) :: H_diab(nbasis,nbasis),delV_dels(nbasis,nbasis,nclass,3)
  integer i,j,k
  real*8 V_AHB,dv1_dx(nclass,3)
  real*8 pot_coul(nbasis),dv2_dx(nbasis,nclass,3)
  real*8 xij(3),rij,x_AB(3)

  H_diab=KE_DVR
  delV_dels=0.d0

  call pbc(x(1,:),x(2,:),x_AB,rAB)

  pot_coul=0.d0;dv2_dx=0.d0

  do j=3,nclass,2
    do i=1,2,2
      call pbc(x(i,:),x(j,:),xij,rij)

      do k=1,nbasis
        rAH=rAH_min+(k-1)*rAH_del
        call compute_charges(rAH)
!        call coul_solv_solute(pot_coul,dv2_dx,xij,rij,i,j,k)
!        if(i==2)call coul_H_solute(pot_coul,dv2_dx,rAH,j,x_AB,k)
        call coul_molecular(pot_coul,dv2_dx,i,j,k)
        call coul_H_solute(pot_coul,dv2_dx,rAH,j,x_AB,k)

        !call pot_AHB(V_AHB,dv1_dx,rAH)
        !H_diab(k,k)=H_diab(k,k)+V_AHB+pot_coul
        !delV_dels(k,k,:,:)=dv1_dx+dv2_dx
      enddo
    enddo
  enddo

  do k=1,nbasis
    rAH=rAH_min+(k-1)*rAH_del
    call pot_AHB(V_AHB,dv1_dx,rAH)
    H_diab(k,k)=H_diab(k,k)+V_AHB+pot_coul(k)
    delV_dels(k,k,:,:)=dv1_dx+dv2_dx(k,:,:)
  enddo

end subroutine compute_potential
!-----------------------------------------------------------------  

subroutine compute_distances
  implicit none
!  integer i,j
  real*8 vec(3)
!
!  do i=1,nclass-1
!    do j=i+1,nclass
!      call pbc(x(i,:),x(j,:),vec,distance(i,j))
!    enddo
!  enddo
!
  call pbc(x(1,:),x(2,:),vec,rAB)
!
end subroutine compute_distances
!-----------------------------------------------------------------  

subroutine pbc(x1,x2,r12,length)
  implicit none
  real*8,intent(in)::x1(3),x2(3)
  real*8,intent(out)::r12(3),length
  integer i
  real*8 tmp
  integer c1,c2

  !call system_clock(c1)

  do i=1,3
    tmp=x2(i)-x1(i)
    if(tmp>box_length(i)/2.d0) tmp=tmp-box_length(i)
    if(tmp<-box_length(i)/2.d0) tmp=tmp+box_length(i)
    r12(i)=tmp
  enddo
  length=dsqrt(sum(r12*r12))

  !call system_clock(c2)
  !!tim_pbc=tim_pbc+(c2-c1)/real(cnt_rate)

end subroutine pbc
!-----------------------------------------------------------------  

subroutine center(x1,x2,cent)
  implicit none
  real*8,intent(in)::x1(3),x2(3)
  real*8,intent(out)::cent(3)
  integer i
  real*8 tmp,su
  integer c1,c2

  !call system_clock(c1)

  do i=1,3
    tmp=dabs(x2(i)-x1(i))
    su=x2(i)+x1(i)
    if(tmp>box_length(i)/2.d0) su=su-box_length(i)
    cent(i)=su/2.d0
  enddo

  !call system_clock(c2)
  !!tim_pbc=tim_pbc+(c2-c1)/real(cnt_rate)

end subroutine center
!-----------------------------------------------------------------  

subroutine pbc_adjust
  implicit none
  integer i,j

  do i=1,nclass
    do j=1,3
      if(x(i,j)<0.d0)x(i,j)=x(i,j)+box_length(j)
      if(x(i,j)>box_length(j))x(i,j)=x(i,j)-box_length(j)
    enddo
  enddo

end subroutine pbc_adjust
!-----------------------------------------------------------------  

subroutine update_verlet_list
  implicit none
  integer i,j,nn
  real*8 xij(3),rij

  do i=1,nclass-1
    nn=0
    do j=i+1,nclass
      call pbc(x(i,:),x(j,:),xij,rij)
      if(rij<Rc_verlet) then
        nn=nn+1
        verlet_list(i,nn+1)=j
      endif
    enddo
    verlet_list(i,1)=nn
  enddo


end subroutine update_verlet_list
!-----------------------------------------------------------------  

SUBROUTINE FCC(rx,ry,rz,ex,ey,ez,n,nc)

! Allen Tildeseley page 169
! Website: ftp://ftp.dl.ac.uk/ccp5/ALLEN_TILDESLEY/F.23

!********************************************************************************
!** FICHE F.23.  ROUTINE TO SET UP ALPHA FCC LATTICE OF LINEAR MOLECULES       **
!** This FORTRAN code is intended to illustrate points made in the text.       **
!** To our knowledge it works correctly.  However it is the responsibility of  **
!** the user to test it, if it is to be used in a research application.        **
!********************************************************************************


!C    *******************************************************************
!C    ** SETS UP THE ALPHA FCC LATTICE FOR N LINEAR MOLECULES.         **
!C    **                                                               **
!C    ** THE SIMULATION BOX IS A UNIT CUBE CENTRED AT THE ORIGIN.      **
!C    ** N SHOULD BE AN INTEGER OF THE FORM ( 4 * ( NC ** 3 ) ),       **
!C    ** WHERE NC IS THE NUMBER OF FCC UNIT CELLS IN EACH DIRECTION.   **
!C    ** SEE FIGURE 5.10 FOR A DIAGRAM OF THE LATTICE AND A            **
!C    ** DEFINITION OF THE FOUR ORIENTATIONAL SUBLATTICES.             **
!C    **                                                               **
!C    ** PRINCIPAL VARIABLES:                                          **
!C    **                                                               **
!C    ** INTEGER N                    NUMBER OF MOLECULES              **
!C    ** REAL    RX(N),RY(N),RZ(N)    MOLECULAR POSITIONS              **
!C    ** REAL    EX(N),EY(N),EZ(N)    UNIT VECTORS GIVING ORIENTATIONS **
!C    ** REAL    RROOT3               1.0 / SQRT ( 3.0 )               **
!C    *******************************************************************

  INTEGER,intent(in) :: n,NC
  REAL*8,intent(out),dimension(n)::  RX,RY, RZ, EX, EY, EZ
  REAL*8      RROOT3

!  PARAMETER ( RROOT3 = 0.5773503 )

  REAL*8      CELL, CELL2
  INTEGER     I, IX, IY, IZ, IREF, M

!C    ** CALCULATE THE SIDE OF THE UNIT CELL **

  rroot3=1.d0/dsqrt(3.d0)

  CELL  = box_length(1) / REAL ( NC )
  CELL2 = 0.5 * CELL

!C    ** BUILD THE UNIT CELL **

!C    ** SUBLATTICE A **

  RX(1) =  0.0
  RY(1) =  0.0
  RZ(1) =  0.0
  EX(1) =  RROOT3
  EY(1) =  RROOT3
  EZ(1) =  RROOT3

!C    ** SUBLATTICE B **

  RX(2) =  CELL2
  RY(2) =  CELL2
  RZ(2) =  0.0
  EX(2) =  RROOT3
  EY(2) = -RROOT3
  EZ(2) = -RROOT3

!C    ** SUBLATTICE C **

  RX(3) =  0.0
  RY(3) =  CELL2
  RZ(3) =  CELL2
  EX(3) = -RROOT3
  EY(3) =  RROOT3
  EZ(3) = -RROOT3

!C    ** SUBLATTICE D **

  RX(4) =  CELL2
  RY(4) =  0.0
  RZ(4) =  CELL2
  EX(4) = -RROOT3
  EY(4) = -RROOT3
  EZ(4) =  RROOT3

!C    ** CONSTRUCT THE LATTICE FROM THE UNIT CELL **

  M = 0

  DO 99 IZ = 1, NC

     DO 98 IY = 1, NC

        DO 97 IX = 1, NC

           DO 96 IREF = 1, 4

              RX(IREF+M) = RX(IREF) + CELL * REAL ( IX - 1 )
              RY(IREF+M) = RY(IREF) + CELL * REAL ( IY - 1 )
              RZ(IREF+M) = RZ(IREF) + CELL * REAL ( IZ - 1 )

              EX(IREF+M) = EX(IREF)
              EY(IREF+M) = EY(IREF)
              EZ(IREF+M) = EZ(IREF)

96               CONTINUE

           M = M + 4

97            CONTINUE

98         CONTINUE

99      CONTINUE

  RETURN
END subroutine fcc
!-----------------------------------------------------------------  

subroutine Boltzmann_velocities(v,ncl,mass,temperature)
  !! Boltzmann velocities
  implicit none
  integer,intent(in)::ncl
  real*8,intent(in)::temperature,mass(ncl)
  real*8,intent(out)::v(ncl,3)
  integer i,j,ndof
  real*8 sig_p,rnd,fac,vcom

  do i=1,ncl
    sig_p=dsqrt(kb*temperature*mass(i))
    do j=1,3
      call gaussian_random_number(rnd)
      v(i,j)=dabs(1.d0/mass(i)*(rnd*sig_p))
    enddo
  enddo

  !! COM velocity
  do j=1,3
    vcom=sum(v(:,j))/real(ncl)
    v(:,j)=v(:,j)-vcom
  enddo

  call rattle(2,dtc)

  !! rescaling
  call evaluate_variables(0)
  call evaluate_variables(1)
  ndof=(3*ncl-6)-(ncl/2-1)  !! includes bond constraints
  fac=dsqrt(ndof*kb*temperature/(KE_en))
  v=v*fac

  call evaluate_variables(1)

end subroutine Boltzmann_velocities
!-----------------------------------------------------------------  

subroutine check_acceleration
  implicit none
  integer i,j,nflag
  real*8 delx,en_old,acc_sav(nclass,3)
  real*8 q0,rnd

  delx=1.d-17
!  state=1
  flag_classical=1


  write(13,'(i5)') nclass
  write(13,*) 0.d0
  write(13,'(A,3f15.7)')'A ',x(1,:)*1.d10
  write(13,'(A,3f15.7)')'B ',x(2,:)*1.d10
  do i=3,nclass,2
    write(13,'(A,3f15.7)')'Me ',x(i,:)*1.d10!,layer_num(i)
    write(13,'(A,3f15.7)')'Cl ',x(i+1,:)*1.d10!,layer_num(i)
    !if(i>1.and.layer_num(i)-layer_num(i-1)>0) write(200,*) 
    !if(i>1.and.layer_num(i)-layer_num(i-1)>0) write(200,*) 
  enddo


  !do i=1,nclass
  !  do j=1,3
  !    call random_number(rnd)
  !    x(i,j)=(rnd*2-1.d0)*1.d-10
  !  enddo
  !enddo

  call evaluate_variables(0)
  en_old=pot_en;acc_sav=acc

  write(6,*) "delx=",delx
  write(6,*)

write(6,*) x(1,:)
write(6,*) x(2,:)
write(6,*) state

  do i=1,nclass
    do j=1,3
      x(i,j)=x(i,j)+delx
      call evaluate_variables(0)
      acc(i,j)=-(pot_en-en_old)/delx/mass(i)
      write(6,*)"Analytical acceleration =",acc_sav(i,j)
      write(6,*)"Numerical acceleration  =",acc(i,j)
      write(6,*)"Error =",(acc(i,j)-acc_sav(i,j))/acc(i,j)*100.d0
      write(6,*)
      x(i,j)=x(i,j)-delx
    enddo
   if(i==3)read(*,*)
  enddo

  stop

end subroutine check_acceleration
!---------------------------------------------------------- 

function inverse(mat)
  complex*16,intent(in) :: mat(2,2)
  complex*16 inverse(2,2),det
  complex*16 a,b,c,d

  a=mat(1,1);b=mat(1,2)
  c=mat(2,1);d=mat(2,2)

  det=a*d-b*c
  inverse(1,1)=d
  inverse(1,2)=-b
  inverse(2,1)=-c
  inverse(2,2)=a

  inverse=inverse/det

end function inverse
!-----------------------------------------------------------------  

function commute(A,B,iflag)
  integer,intent(in) :: iflag
  real*8,intent(in) :: A(:,:)
  complex*16,intent(in) :: B(:,:)
  complex*16 commute(nquant,nquant)
  real*8 delA(nquant,nquant)
  complex*16 tmp
  integer j,k

  if(iflag==0) commute=matmul_lap(A,B)-matmul_lap(B,A)

  if(iflag==1) then
    !! Assume A is diagonal
    do j=1,nquant
      do k=1,nquant
        commute(j,k)=B(j,k)*(A(j,j)-A(k,k))
      enddo
    enddo
  endif

  if(iflag==2) then
    !! Assume A is tridiagonal, with a_ii=0, and a_ij=-a_ji (a is assumed to be d_ij)
    do j=1,nquant
      do k=1,nquant
        tmp=0.d0
        if(j<nquant) tmp=tmp+A(j,j+1)*B(j+1,k)
        if(j>1) tmp=tmp-A(j-1,j)*B(j-1,k)
        if(k>1) tmp=tmp-A(k-1,k)*B(j,k-1)
        if(k<nquant) tmp=tmp+A(k,k+1)*B(j,k+1)
      enddo
    enddo
  endif

end function commute
!-----------------------------------------------------------------  

function anti_commute(A,B,iflag)
  integer,intent(in) :: iflag
  real*8,intent(in) :: A(:,:)
  complex*16,intent(in) :: B(:,:)
  complex*16 anti_commute(nquant,nquant)
  real*8 delA(nquant,nquant)
  integer i,j

  if(iflag==0) anti_commute=matmul_lap(A,B)+matmul_lap(B,A)

  if(iflag==1) then
    !! Assume A is diagonal
    do i=1,nquant
      do j=1,nquant
       anti_commute(i,j)=B(i,j)*(A(i,i)+A(j,j))
      enddo
    enddo
  endif

end function anti_commute
!-----------------------------------------------------------------  

subroutine draw_pes
  implicit none
  integer i,j,n
  real*8 pot,dpot_dx(nclass,3),rAH

  state=1
  open(10,file="pes.out")
  do i=1,100
    rAH=0.8d-10+1.d-10*i/100.d0
    call pot_AHB(pot,dpot_dx,rAH)
    write(10,*) rAH*1.d10,pot/wave_to_J
  enddo
  close(10)

write(6,*) r_exp
call draw_wavefn

  stop

end subroutine draw_pes
!-----------------------------------------------------------------  

subroutine draw_wavefn
  implicit none
  integer i
  real*8 V1,dv1_dx(nclass)
  real*8 V2,dv2_dx(nclass)

  !write(23,*)x(2)*1.d10
  do i=1,nbasis
    !call pot_q(V1,dv1_dx,r_grid(i))
    !call pot_coup(V2,dv2_dx,r_grid(i))
    !write(22,*)r_grid(i),(V1+V2)/wave_to_J
    write(23,*)rAH_grid(i)*1.d10,si_adiab(i,1:2)
  enddo
  !write(22,*);write(22,*)
  write(23,*);write(23,*)

end subroutine draw_wavefn
!-----------------------------------------------------------------  

function gamma_fn(z)
  !http://www1.uprh.edu/rbaretti/GammaFunction7dic2010.htm
  implicit none
  integer nstep,i
  real*8 ti,tf,dt,t
  complex*16 sum,gamma_fn,z,f

  data ti,tf,nstep /0.d0,30.d0,20001/

  f(t)=t**(z)*dexp(-t)

  dt=(tf-ti)/dfloat(nstep)
!  sum=dt**z/z -dt**(z+1.d0)/(z+1.d0)
  sum=0.d0

!  do i=2,nstep,2
  do i=1,nstep
    t=ti+dt*dfloat(i)
!    sum=sum+(dt/3.d0)*(f(t-dt)+ 4.d0*f(t) +f(t+dt))
    sum=sum+f(t)
  enddo

  gamma_fn=sum*dt/z

end function gamma_fn
!---------------------------------------------------------- 

pure function arg(z)
  implicit none
  real*8 arg
  complex*16,intent(in) :: z
  real*8 zr,zi

  zr=dble(z)
  zi=dimag(z)
  arg=datan2(zi,zr)
!  arg=datan(zi/zr)

end function arg
!----------------------------------------------------------   

subroutine diag(mat,n,eigen_value,eigen_vect,m_values)
  !! Diaganalizing matrix using dsyevr. First m_values eigen values and eigenvectors computed.
  !! The module's common variables should contain:

  !! Initialize nold=0 

  !! nold makes sure that everytime value of n changes, work and iwork are re-allocated for optimal performance.
  !! mat is destroyed after use.

  implicit none
  integer,intent(in) :: n,m_values
  real*8,intent(out) :: eigen_value(n),eigen_vect(n,m_values)
  real*8,intent(inout) :: mat(n,n)
  real*8 vl,vu,abstol
  integer il,iu,info,m,AllocateStatus
  integer lwork,liwork

  vl=0.d0;vu=0.d0   !! not referenced
  il=1;iu=m_values
  abstol=0.d0
  info=0

  if(nold.ne.n .or. .not.allocated(work) .or. .not.allocated(iwork) .or. .not.allocated(isuppz)) then
  !if(nold.ne.n) then
    lwork=-1;liwork=-1
    if(allocated(isuppz))deallocate(isuppz)
    if(allocated(work))deallocate(work)
    if(allocated(iwork))deallocate(iwork)
    allocate(isuppz(2*m_values),work(n),iwork(n))
    call dsyevr('V','I','U',n,mat,n,vl,vu,il,iu,abstol,m,eigen_value,eigen_vect,n,isuppz,work,lwork,iwork,liwork,info)
    lwork=nint(work(1)); liwork=iwork(1)
    deallocate(work,iwork)
    allocate(work(lwork),STAT=AllocateStatus)
    if(allocatestatus.ne.0) write(6,*)"problem in diag, allocation"
    allocate(iwork(liwork),STAT=AllocateStatus)
    if(allocatestatus.ne.0) write(6,*)"problem in diag, allocation"
    nold=n
  endif

  lwork=size(work)
  liwork=size(iwork)

  call dsyevr('V','I','U',n,mat,n,vl,vu,il,iu,abstol,m,eigen_value,eigen_vect,n,isuppz,work,lwork,iwork,liwork,info)
  if(info.ne.0) then
    write(6,*) "problem in diagonalization",info
    stop
  endif

end subroutine diag
!---------------------------------------------------------- 

subroutine logm(mat,log_mat,n)
  !! http://arxiv.org/pdf/1203.6151v4.pdf
  implicit none
  integer,intent(in):: n
  real*8,intent(in):: mat(n,n)
  real*8,intent(out):: log_mat(n,n)
  integer i
  complex*16 T(n,n),en(n),vect(n,n)
  complex*16 dd(n,n)

  call schur(mat,T,n,en,vect,nold,cwork)

  dd=0.d0
  do i=1,n
    dd(i,i)=cdlog(t(i,i)/cdabs(t(i,i)))
  enddo

  log_mat=matmul_lap(vect,matmul_lap(dd,conjg(transpose(vect))))

end subroutine logm
!-----------------------------------------------------------------  

subroutine inverse_squareroot(mat,n)
  !! http://arxiv.org/pdf/1203.6151v4.pdf
  implicit none
  integer,intent(in):: n
  real*8,intent(inout):: mat(n,n)
  integer i
  complex*16 T(n,n),en(n),vect(n,n)
  complex*16 dd(n,n)

  call schur(mat,T,n,en,vect,nold,cwork)

  dd=0.d0
  do i=1,n
    dd(i,i)=1.d0/t(i,i)**0.5d0
  enddo

  mat=matmul_lap(vect,matmul_lap(dd,conjg(transpose(vect))))

end subroutine inverse_squareroot
!-----------------------------------------------------------------  

subroutine schur(mat,T,n,eigen_value,eigen_vect,nold,cwork)
  !! Diaganalizing matrix using dsyevr. First m_values eigen values and eigenvectors computed.
  !! The module's common variables should contain:

  !! Initialize nold=0 

  !! nold makes sure that everytime value of n changes, work and iwork are re-allocated for optimal performance.
  !! mat is destroyed after use.

  implicit none
  integer,intent(in) :: n
  integer,intent(inout) :: nold
  complex*16,intent(out) :: eigen_value(n),eigen_vect(n,n)
  real*8,intent(in) :: mat(n,n)
  complex*16,intent(out) :: T(n,n)
  complex*16,allocatable,intent(inout):: cwork(:)
  real*8 rwork(n)
  complex*16 mat_c(n,n)

  integer lwork
  logical:: select
  logical bwork(n)
  integer sdim,info,AllocateStatus

  T=mat

  info=0
  sdim=0

  if(nold.ne.n .or. .not.allocated(cwork)) then
  !if(nold.ne.n) then
    lwork=-1
    if(allocated(cwork))deallocate(cwork)
    allocate(cwork(n))
    call zgees('V','N',SELECT,N,T,n,SDIM,eigen_value,eigen_vect,n,cWORK,LWORK,rwork,BWORK,INFO)
    lwork=int(cwork(1))
    deallocate(cwork)
    allocate(cwork(lwork),STAT=AllocateStatus)
    if(allocatestatus.ne.0) write(6,*)"problem in diag, allocation"
    nold=n
  endif

  lwork=size(cwork)
  call zgees('V','N',SELECT,N,T,n,SDIM,eigen_value,eigen_vect,n,cWORK,LWORK,rwork,BWORK,INFO)
  if(info.ne.0) then
    write(6,*) "problem in diagonalization",info
    stop
  endif

end subroutine schur
!---------------------------------------------------------- 

REAL FUNCTION determinant(matrix, n)
    !!http://web.hku.hk/~gdli/UsefulFiles/Example-Fortran-program.html
    IMPLICIT NONE
    REAL*8, DIMENSION(n,n) :: matrix
    INTEGER, INTENT(IN) :: n
    REAL*8 :: m, temp
    INTEGER :: i, j, k, l
    LOGICAL :: DetExists = .TRUE.
    l = 1
    !Convert to upper triangular form
    DO k = 1, n-1
        IF (matrix(k,k) == 0) THEN
            DetExists = .FALSE.
            DO i = k+1, n
                IF (matrix(i,k) /= 0) THEN
                    DO j = 1, n
                        temp = matrix(i,j)
                        matrix(i,j)= matrix(k,j)
                        matrix(k,j) = temp
                    END DO
                    DetExists = .TRUE.
                    l=-l
                    EXIT
                ENDIF
            END DO
            IF (DetExists .EQV. .FALSE.) THEN
                determinant= 0
                return
            END IF
        ENDIF
        DO j = k+1, n
            m = matrix(j,k)/matrix(k,k)
            DO i = k+1, n
                matrix(j,i) = matrix(j,i) - m*matrix(k,i)
            END DO
        END DO
    END DO
    
    !Calculate determinant by finding product of diagonal elements
    determinant= l
    DO i = 1, n
        determinant= determinant* matrix(i,i)
    END DO
    
END FUNCTION determinant
!-----------------------------------------------------------------  

subroutine state_of_the_system(iflag)
  !! My dear particles and wavefunctions
  !! This subroutine writes/reads the state of the sytem
  implicit none
  integer,intent(in)::iflag

  open(50,file='state.out')
  if(iflag==0) then
    write(50,*)state
    write(50,*)ci
    write(50,*)x
    write(50,*)v
  endif
  if(iflag==1) then
    read(50,*)state
    read(50,*)ci
    read(50,*)x
    read(50,*)v
    call evaluate_variables(0)
    call evaluate_variables(1)
  endif
  close(50)

end subroutine state_of_the_system
!-----------------------------------------------------------------  

subroutine set_charges(frac)
  implicit none
  real*8,intent(in)::frac
  integer i

  if(frac==0.d0) charge=0

  if(frac>0.d0) then
    do i=3,nclass,2
      charge(i)=0.25d0*frac
    enddo
    do i=4,nclass,2
      charge(i)=-0.25d0*frac
    enddo
  endif

end subroutine set_charges
!-----------------------------------------------------------------  

subroutine write_vmd
  implicit none
  integer i

  write(50,'(i5)') nclass
  write(50,*) curr_time*1.d15
  write(50,'(A,3f15.7)')'A ',x(1,:)*1.d10
  write(50,'(A,3f15.7)')'N ',x(2,:)*1.d10
  do i=3,nclass,2
    write(50,'(A,3f15.7)')'Me ',x(i,:)*1.d10!,layer_num(i)
    write(50,'(A,3f15.7)')'Cl ',x(i+1,:)*1.d10!,layer_num(i)
    !if(i>1.and.layer_num(i)-layer_num(i-1)>0) write(200,*) 
    !if(i>1.and.layer_num(i)-layer_num(i-1)>0) write(200,*) 
  enddo

end subroutine write_vmd
!-----------------------------------------------------------------  

subroutine stochastic_force(delr,delv,dt)
  implicit none
  real*8,intent(in)::dt
  real*8,intent(out) :: delr(nclass,3),delv(nclass,3)!f(nclass)
  integer i,j
  real*8 rnd1,rnd2,sig_r,sig_v,sig_rv,gdt

  gdt=gamma_B*dt

  do i=1,nclass
    do j=1,3
      sig_r=dt*dsqrt(kb*temperature_init/mass(i) *1.d0/gdt*(2-1.d0/gdt*(3-4*dexp(-gdt)+dexp(-2*gdt))))
      sig_v=dsqrt(kb*temperature_init/mass(i)*(1-dexp(-2*gdt)))
      sig_rv=(dt*kb*temperature_init/mass(i)* 1.d0/gdt *(1-dexp(-gdt))**2)/(sig_r*sig_v)  !! correlation coeffecient

      call gaussian_random_number(rnd1)
      call gaussian_random_number(rnd2)
      delr(i,j)=sig_r*rnd1
      delv(i,j)=sig_v*(sig_rv*rnd1+dsqrt(1-sig_rv**2)*rnd2)
    enddo
  enddo

  do j=1,3
    delr(:,j)=delr(:,j)-sum(delr(:,j))/dfloat(nclass)
    delv(:,j)=delv(:,j)-sum(delv(:,j))/dfloat(nclass)
  enddo

end subroutine stochastic_force
!-----------------------------------------------------------------  

End Module mod_afssh
