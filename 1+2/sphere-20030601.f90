!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!! sphere.f90 : main engine and calculate right hand side
!


module sphere_module
  use sparse_matrix_module
  use amg_solver_module
  use triangulation_module
  use h1_fem_space_module
  use mpi

  implicit none

  integer,parameter::nz=511
  type(triangulation)::tri
  type(h1_fem_space),SAVE::fem_space
  type(h1_fem_function), dimension(:), allocatable::f
  real, dimension(0:nz)::zcord
  real, dimension(0:nz)::vg
  real, dimension(0:nz+1)::v
  real, dimension(1:nz)::a
  real, dimension(1:nz)::b
  real, dimension(1:nz)::rhs1
  real, dimension(1:nz)::rhs2
  real, dimension(1:3)::w
  real, dimension(1:102)::work
  real, dimension(0:nz,1:3,1:3)::uu
  real, dimension(0:nz,1:3,1:3,1:3,1:3)::uuuu 
  real, dimension(0:nz,1:3,1:3)::tau
  real, dimension(0:nz):: N1
  real, dimension(0:nz):: N2
  real, dimension(0:nz):: shearstress
  real, dimension(0:nz)::dens
  real, dimension(0:nz)::conv_g_dens
  real, dimension(0:nz,1:3,1:3)::conv_g_uu
  real, dimension(:),allocatable::Pot
  real, dimension(:),allocatable::tmp1
  real, dimension(:),allocatable:: sendbuffer
  real, dimension(:),allocatable:: recvbuffer
  integer,dimension(:),allocatable:: displs
  integer,dimension(:),allocatable:: recvcounts
  real,dimension(1:2*(nz+1)):: rmollif,temp,conv_tmp
  real time,dt,tlast,dz,CFL2,bigU,De,Re,gamma,Pressure,epsilon,d
  integer myrank, nprocs, ierr, njobs, ljob, rjob, npack, SRC, DEST, TAG
  INTEGER STATUS(MPI_STATUS_SIZE)

contains

  subroutine run(filename)
    character(*), intent(in)::filename

    type(sparse_matrix)::mass_matrix
    type(sparse_matrix)::stiff_matrix
    type(sparse_matrix)::M
    type(amg_solver)::solver
    real, dimension(:), allocatable::rhs
    real f_l1_norm, error
    real,dimension(1:3)::u
    integer i,j,k,l,p,q,z,iter,itmp,left,right

    De = 1.
    Re = 1.
    gamma = 5./9.
    bigU = 5.33*1.5
    Pressure = -40.
    epsilon = 0.01
    dz = 1.0/nz
    dt = 0.001
    tlast = 1000.
    CFL2=dt/dz/dz

    !Begin init mpi environment
    call mpi_init(ierr)
    call mpi_comm_rank(mpi_comm_world, myrank, ierr)
    call mpi_comm_size(mpi_comm_world, nprocs, ierr)
    !End init mpi environment

    !Begin assign jobs
    l = mod(nz+1,nprocs)
    k = (nz+1-l)/nprocs
    if (myrank.lt.l) then
       njobs = k+1
       rjob = (myrank+1)*(k+1)
    else
       njobs = k
       rjob = (nz+1)-(nprocs-(myrank+1))*k
    end if
    rjob = rjob-1
    ljob = rjob-njobs+1
    allocate(f(ljob-1:rjob+1))
    allocate(sendbuffer(1:10*njobs))
    allocate(recvbuffer(0:10*(nz)-1))
    allocate(displs(0:nprocs-1))
    allocate(recvcounts(0:nprocs-1))
    allocate(Pot(ljob-1:rjob+1))
    allocate(tmp1(ljob:rjob+1))
    !End assign jobs

    !Begin initilize finite element and solver
    call triangulation_read_mesh(tri, filename)
    call h1_fem_space_initialize(fem_space, tri)
    do i=ljob-1,rjob+1
       call h1_fem_function_initialize(f(i), fem_space)
    end do
    call initial_value()
    call h1_fem_space_get_mass_matrix(fem_space, mass_matrix)
    call h1_fem_space_get_stiff_matrix(fem_space, stiff_matrix)
    call sparse_matrix_initialize(M, fem_space%spH)
    M%matrix_entry = (1./dt)*mass_matrix%matrix_entry &
         + 1./De*stiff_matrix%matrix_entry    
    call amg_solver_initialize(solver, M)
    !End initilize finite element and solver

    allocate(rhs(1:tri%n_node))

    if (myrank.eq.0) then 
       open(unit=11, file="v-5.plt", action="write")
       open(unit=12, file="dens-5.plt", action="write")
       open(unit=13, file="N1-5.plt", action="write")
       open(unit=14, file="N2-5.plt", action="write")
       open(unit=15, file="shearstress-5.plt", action="write")
       open(unit=16, file="Reb22-5.plt", action="write")
       open(unit=17, file="Imb22-5.plt", action="write")
    end if

    if (myrank.eq.0) print *, "start iteration"
    time = 0.
    do iter=1,200000
       !==========================compute and dump data=================!
       call get_dens_uu_tau()
       !Begin packing, gathering, and unpacking the data
       npack = 10
       l = mod(nz+1,nprocs)
       k = (nz+1-l)/nprocs
       do i=0,nprocs-1
          if (i.lt.l) then
             recvcounts(i) = npack*(k+1)
             displs(i) = npack*i*(k+1) 
          else
             recvcounts(i) = npack*k
             displs(i) = npack*((nz+1)-(nprocs-i)*k)
          end if
       end do
       do i=1,njobs
          sendbuffer(i) = dens(ljob-1+i)
          sendbuffer(i+njobs) = uu(ljob-1+i,1,1)
          sendbuffer(i+2*njobs) = uu(ljob-1+i,2,2)
          sendbuffer(i+3*njobs) = uu(ljob-1+i,3,3)
          sendbuffer(i+4*njobs) = uu(ljob-1+i,1,2)
          sendbuffer(i+5*njobs) = uu(ljob-1+i,1,3)
          sendbuffer(i+6*njobs) = uu(ljob-1+i,2,3)
          sendbuffer(i+7*njobs) = tau(ljob-1+i,1,1)-tau(ljob-1+i,3,3)
          sendbuffer(i+8*njobs) = tau(ljob-1+i,3,3)-tau(ljob-1+i,2,2)
          sendbuffer(i+9*njobs) = tau(ljob-1+i,1,3)
       end do
       call MPI_Allgatherv(sendbuffer, npack*njobs, MPI_REAL8, &
            recvbuffer, recvcounts, displs, MPI_REAL8, &
            MPI_Comm_World, ierr)
       do i=0,nprocs-1
          if (i.lt.l) then
             itmp = k+1
          else
             itmp = k
          end if
          do j=0,itmp-1
             dens(displs(i)/npack+j) = recvbuffer(displs(i)+j)
             uu(displs(i)/npack+j,1,1) = recvbuffer(displs(i)+itmp+j)
             uu(displs(i)/npack+j,2,2) = recvbuffer(displs(i)+2*itmp+j)
             uu(displs(i)/npack+j,3,3) = recvbuffer(displs(i)+3*itmp+j)
             uu(displs(i)/npack+j,1,2) = recvbuffer(displs(i)+4*itmp+j)
             uu(displs(i)/npack+j,1,3) = recvbuffer(displs(i)+5*itmp+j)
             uu(displs(i)/npack+j,2,3) = recvbuffer(displs(i)+6*itmp+j)
             N1(displs(i)/npack+j) = recvbuffer(displs(i)+7*itmp+j)
             N2(displs(i)/npack+j) = recvbuffer(displs(i)+8*itmp+j)
             shearstress(displs(i)/npack+j) =  recvbuffer(displs(i)+9*itmp+j)
          end do
       end do
       do i=0,nz
          do p=1,3
             do q=1,p-1
                uu(i,p,q) = uu(i,q,p)
             end do
          end do
       end do
       !End packing, gathering, and unpacking the data

       !Begin dumping data to disk
       if (mod(iter,10).eq.0) then 
          if (myrank.eq.0) then          
             write(unit=11, fmt=*) "ZONE"
             write(unit=12, fmt=*) "ZONE"
             write(unit=13, fmt=*) "ZONE"
             write(unit=14, fmt=*) "ZONE"
             write(unit=15, fmt=*) "ZONE"
             write(unit=16, fmt=*) "ZONE"
             write(unit=17, fmt=*) "ZONE"
             do i=0,nz
                write(unit=11, fmt=*) zcord(i), (v(i)+v(i+1))/2.
                write(unit=12, fmt=*) zcord(i), dens(i)
                write(unit=13, fmt=*) zcord(i), N1(i)
                write(unit=14, fmt=*) zcord(i), N2(i)
                write(unit=15, fmt=*) zcord(i), shearstress(i)
                write(unit=16, fmt=*) zcord(i), uu(i,1,1)-uu(i,2,2)
                write(unit=17, fmt=*) zcord(i), uu(i,1,2)
             end do
          end if
       end if
       !End dumping data to disk

       !=========================Update the flow========================!
       !!implicit scheme
       !do i=1,nz
       !   a(i) = -(1.-gamma)/Re*CFL2
       !   c(i) = a(i)
       !   b(i) = 1.-2.*a(i)
       !   d(i) = v(i)-dt*Pressure+dt*gamma/De/Re*(shearstress(i)-shearstress(i-1))/dz
       !end do
       !b(1) = 1.-3.*a(1)
       !b(nz) = 1.-3.*a(nz)
       !call tridag(a,b,c,d,v(1:nz-1),nz-1)
       
       do i=1,nz
          rhs1(i) = v(i)-dt*Pressure+dt*Gamma/De/Re*(shearstress(i)-shearstress(i-1))/dz
       end do

       !Cholesky decomposition for symetric tridiagonal matrix    
       d = (1.-Gamma)/Re*CFL2
       b(1) = 0.
       a(1) = sqrt(1.+3.*d)
       do i=2,nz-1
          b(i) = -d/a(i-1)
          a(i) = sqrt(1.+2.*d-b(i)*b(i))
       end do

       b(nz) = -d/a(nz-1)
       a(nz) = sqrt(1.+3.*d-b(nz)*b(nz))

       rhs2(1) = rhs1(1)/a(1)
       do i=2,nz
          rhs2(i) = (rhs1(i)-b(i)*rhs2(i-1))/a(i)
       end do

       v(nz) = rhs2(nz)/a(nz)
       do i=nz-1,1,-1
          v(i) = (rhs2(i)-b(i+1)*v(i+1))/a(i)
       end do

       v(0)=-v(1)
       v(nz+1)=-v(nz)

       do i=0,nz
          vg(i) = (v(i+1)-v(i))/dz
       end do

       !===============Half step of finite element method===============!
       call get_conv_g_uu()
       do i = ljob,rjob
          call get_right_hand_side(rhs, i)
          call sparse_matrix_vector_mult_ps(mass_matrix, f(i)%value, 1./dt, rhs)
          call amg_solver_solve(solver, f(i)%value, rhs)
       end do

       !=============Half step of finite difference method==============!       
       call get_dens_uu()
       !Begin packing, gathering, and unpacking the data
       npack = 7
       l = mod(nz+1,nprocs)
       k = (nz+1-l)/nprocs
       do i=0,nprocs-1
          if (i.lt.l) then
             recvcounts(i) = npack*(k+1)
             displs(i) = npack*i*(k+1) 
          else
             recvcounts(i) = npack*k
             displs(i) = npack*((nz+1)-(nprocs-i)*k)
          end if
       end do
       do i=1,njobs
          sendbuffer(i) = dens(ljob-1+i)
          sendbuffer(i+njobs) = uu(ljob-1+i,1,1)
          sendbuffer(i+2*njobs) = uu(ljob-1+i,2,2)
          sendbuffer(i+3*njobs) = uu(ljob-1+i,3,3)
          sendbuffer(i+4*njobs) = uu(ljob-1+i,1,2)
          sendbuffer(i+5*njobs) = uu(ljob-1+i,1,3)
          sendbuffer(i+6*njobs) = uu(ljob-1+i,2,3)
       end do
       call MPI_Allgatherv(sendbuffer, npack*njobs, MPI_REAL8, &
            recvbuffer, recvcounts, displs, MPI_REAL8, &
            MPI_Comm_World, ierr)
       do i=0,nprocs-1
          if (i.lt.l) then
             itmp = k+1
          else
             itmp = k
          end if
          do j=0,itmp-1
             dens(displs(i)/npack+j) = recvbuffer(displs(i)+j)
             uu(displs(i)/npack+j,1,1) = recvbuffer(displs(i)+itmp+j)
             uu(displs(i)/npack+j,2,2) = recvbuffer(displs(i)+2*itmp+j)
             uu(displs(i)/npack+j,3,3) = recvbuffer(displs(i)+3*itmp+j)
             uu(displs(i)/npack+j,1,2) = recvbuffer(displs(i)+4*itmp+j)
             uu(displs(i)/npack+j,1,3) = recvbuffer(displs(i)+5*itmp+j)
             uu(displs(i)/npack+j,2,3) = recvbuffer(displs(i)+6*itmp+j)
          end do
       end do
       do i=0,nz
          do p=1,3
             do q=1,p-1
                uu(i,p,q) = uu(i,q,p)
             end do
          end do
       end do
       !End packing, gathering, and unpacking the data

       !Begin send and recv data to neighboring processors
       SRC = MYRANK-1
       IF ( SRC .LT. 0 ) SRC = NPROCS - 1
       DEST = MYRANK + 1
       IF ( DEST .GE. NPROCS ) DEST = 0

       TAG = 111
       CALL MPI_SENDRECV( f(rjob)%value, tri%n_node, MPI_REAL8, &
            DEST, TAG,  f(ljob-1)%value, tri%n_node, MPI_REAL8, &
            SRC, TAG, MPI_COMM_WORLD, STATUS, IERR )

       SRC = MYRANK + 1
       IF ( SRC .GE. NPROCS ) SRC = 0
       DEST = MYRANK-1
       IF ( DEST .LT. 0 ) DEST = NPROCS - 1

       TAG = 222
       CALL MPI_SENDRECV( f(ljob)%value, tri%n_node, MPI_REAL8, &
            DEST, TAG,  f(rjob+1)%value, tri%n_node, MPI_REAL8, &
            SRC, TAG, MPI_COMM_WORLD, STATUS, IERR )
       !End send and recv data to neighboring processors

       !Begin half step of finite difference method
       call get_conv_g_uu()
       call get_conv_g_dens()
       do j=1,tri%n_node
          u = tri%node(j)%x
          left = ljob-1
          if (left.lt.0) then
             left = 0
             tmp1(0) = 0.
          end if
          right = rjob+1
          if (right.gt.nz) then
             right = nz
             tmp1(nz+1) = 0.
          end if
          Pot=0.
          do z=left,right
             do p=1,3
                do q=1,3
                   Pot(z)=Pot(z)-bigU*u(p)*u(q)*conv_g_uu(z,p,q)
                end do
             end do
             Pot(z)=Pot(z)+bigU*conv_g_dens(z)
          end do
          do z=left+1,right
             tmp1(z) = (f(z)%value(j)-f(z-1)%value(j))/dz+&
                  (Pot(z)-Pot(z-1))/dz*(f(z)%value(j)+f(z-1)%value(j))/2.
          end do
          do z=ljob,rjob
             f(z)%value(j)=f(z)%value(j)+dt*epsilon*epsilon/De*&
                  (tmp1(z+1)-tmp1(z))/dz
          end do
       end do
       !End half step of finite difference method
       time = time + dt
       if (time > tlast) then
          exit
       end if

    end do

    close(unit=11)
    close(unit=12)
    close(unit=13)
    close(unit=14)
    close(unit=15)
    close(unit=16)
    close(unit=17)
    call sparse_matrix_destroy(mass_matrix)
    call sparse_matrix_destroy(stiff_matrix)
    call sparse_matrix_destroy(M)
    call amg_solver_destroy(solver)
    call h1_fem_space_destroy(fem_space)
    do i = ljob-1, rjob+1
       call h1_fem_function_destroy(f(i))
    end do
    call triangulation_destroy(tri)
    deallocate(sendbuffer)
    deallocate(recvbuffer)
    deallocate(recvcounts)
    deallocate(displs)
    deallocate(Pot)
    deallocate(tmp1)
    deallocate(rhs)
    deallocate(f)
    call mpi_finalize(ierr)
  end subroutine run

  subroutine get_right_hand_side(rhs, z)
    real, dimension(:), intent(out)::rhs
    integer, intent(in)::z

    integer::i, j, k, l, p, q
    real, dimension(1:3,1:fem_space%n_quadrature_point)::quadrature_point
    real, dimension(1:fem_space%n_quadrature_point)::Jxw
    real, dimension(1:fem_space%n_quadrature_point)::f_on_quad_point
    real, dimension(1:3,1:3)::basis_gradient
    real, dimension(1:3)::grad_V, temp
    real r
    real, dimension(1:3)::cross_x_Kx

    rhs = 0.
    do i=1,tri%n_element
       call h1_fem_space_get_element_info(fem_space, i, quadrature_point, Jxw)
       call h1_fem_function_on_quad_point(f(z), i, f_on_quad_point)
       do j=1,fem_space%n_quadrature_point
          r = sqrt(dot_product(quadrature_point(:,j),quadrature_point(:,j)))
          temp = matmul((1./r)*quadrature_point(:,j), conv_g_uu(z,:,:))
          call cross_product((-2./r)*quadrature_point(:,j), temp, grad_V)
          call h1_fem_space_basis_gradient(fem_space, i, quadrature_point(:,j), basis_gradient)          
          cross_x_Kx(1) = 0.
          cross_x_Kx(2) = vg(z)*quadrature_point(3,j)*quadrature_point(3,j)
          cross_x_Kx(3) = -vg(z)*quadrature_point(3,j)*quadrature_point(2,j)
          do k=1,3
             l = tri%element(i)%node(k)
             rhs(l) = rhs(l) &
                  -1./De*bigU*Jxw(j)*f_on_quad_point(j)* &
                  dot_product(grad_V,basis_gradient(k,:)) &
                  +Jxw(j)*f_on_quad_point(j)* &
                  dot_product(cross_x_Kx,basis_gradient(k,:)) 
          end do
       end do
    end do
  end subroutine get_right_hand_side

  subroutine get_dens_uu()
    real, dimension(1:3,1:fem_space%n_quadrature_point)::quadrature_point
    real, dimension(1:fem_space%n_quadrature_point)::Jxw
    real, dimension(1:fem_space%n_quadrature_point)::f_on_quad_point
    real r
    integer z,i,j,k,l

    dens = 0.
    uu = 0.
    do z=ljob,rjob
       do i=1,tri%n_element
          call h1_fem_space_get_element_info(fem_space, i, quadrature_point,Jxw)
          call h1_fem_function_on_quad_point(f(z), i, f_on_quad_point)
          dens(z) = dens(z) + dot_product(Jxw, f_on_quad_point)
          do j=1,fem_space%n_quadrature_point
             r = sqrt(dot_product(quadrature_point(:,j),quadrature_point(:,j)))
             do k=1,3
                do l=k,3
                   uu(z,k,l) = uu(z,k,l) + Jxw(j)* &
                        (quadrature_point(k,j)/r)*(quadrature_point(l,j)/r)* &
                        f_on_quad_point(j)
                end do
                do l=1,k-1
                   uu(z,k,l) = uu(z,l,k)
                end do
             end do
          end do
       end do
    end do
  end subroutine get_dens_uu

  subroutine get_dens_uu_tau()
    real, dimension(1:3,1:fem_space%n_quadrature_point)::quadrature_point
    real, dimension(1:fem_space%n_quadrature_point)::Jxw
    real, dimension(1:fem_space%n_quadrature_point)::f_on_quad_point
    real  r
    integer z,i,j,k,l,p,q

    dens = 0.
    uu = 0.
    uuuu = 0.
    tau=0.
    do z=ljob,rjob
       do i=1,tri%n_element
          call h1_fem_space_get_element_info(fem_space,i,quadrature_point,Jxw)
          call h1_fem_function_on_quad_point(f(z), i, f_on_quad_point)
          dens(z) = dens(z) + dot_product(Jxw, f_on_quad_point)
          do j=1,fem_space%n_quadrature_point
             r = sqrt(dot_product(quadrature_point(:,j),quadrature_point(:,j)))
             do k=1,3
                do l=k,3
                   uu(z,k,l) = uu(z,k,l) + Jxw(j)* &
                        (quadrature_point(k,j)/r)*(quadrature_point(l,j)/r)* &
                        f_on_quad_point(j)
                   do p=1,3
                      do q=p,3
                         uuuu(z,p,q,k,l)=uuuu(z,p,q,k,l)+Jxw(j)* &
                              (quadrature_point(k,j)/r)*(quadrature_point(l,j)/r)* &                    
                              (quadrature_point(p,j)/r)*(quadrature_point(q,j)/r)* &                     
                              f_on_quad_point(j)                      
                      end do
                      do q=1,p-1
                         uuuu(z,p,q,k,l)=uuuu(z,q,p,k,l)
                      end do
                   end do
                end do
                do l=1,k-1
                   uu(z,k,l) = uu(z,l,k)
                   uuuu(z,:,:,k,l)=uuuu(z,:,:,l,k)
                end do
             end do
          end do
       end do

       do p=1,3
          do q=p,3
             tau(z,p,q)=3.*uu(z,p,q)+0.5*De*vg(z)*uuuu(z,p,q,3,1)
             do k=1,3
                do l=1,3
                   tau(z,p,q)=tau(z,p,q)+2.*bigU*uu(z,k,l)*uuuu(z,p,q,k,l)
                end do
             end do
             do k=1,3
                tau(z,p,q)=tau(z,p,q)-2.*bigU*uu(z,p,k)*uu(z,q,k)
             end do
          end do
          do q=1,p-1
             tau(z,p,q)=tau(z,q,p)
          end do
       end do
       do p=1,3
          tau(z,p,p) = tau(z,p,p)-1.
       end do
    end do
  end subroutine get_dens_uu_tau

  subroutine get_conv_g_dens()
    integer:: i

    do i=1,nz+1
       temp(i)=0.
       temp(i+nz+1)=dens(i-1)
    end do

    call convol(rmollif,temp,2*(nz+1),conv_tmp)

    do i=0,nz
       conv_g_dens(i)=conv_tmp(i+1)/nz 
    end do

  end subroutine get_conv_g_dens

  subroutine get_conv_g_uu()
    integer:: i,p,q
    do p=1,3
       do q=p,3
          do i=1,nz+1
             temp(i)=0.
             temp(i+nz+1)=uu(i-1,p,q)
          end do

          call convol(rmollif,temp,2*(nz+1),conv_tmp)

          do i=0,nz
             conv_g_uu(i,p,q)=conv_tmp(i+1)/nz 
          end do
       end do
       
       do q=1,p-1
          do i=0,nz
             conv_g_uu(i,p,q)=conv_g_uu(i,q,p)
          end do
       end do
    end do
  end subroutine get_conv_g_uu


  !Compute the mollifier kernel and its differentiate
  subroutine kernel()
    integer,parameter::ndiv=800
    real,dimension(ndiv)::tarr
    real,dimension(ndiv)::xcor1
    real,dimension(1:2*(nz+1))::xcor2
    real  dx, tmp, rinteg
    integer i

    dx=2./ndiv
    do i=1,ndiv
       xcor1(i) = (i-0.5)*dx-1.
    end do

    do i=1,ndiv
       tarr(i) = exp(1./(xcor1(i)**2-1))
    end do

    !compute the constant of the mollifier by midpoint formula.

    rinteg = 0.
    do i=1,ndiv
       rinteg = rinteg + tarr(i)
    end do


    rinteg = 1./(rinteg*dx)

    do i=1,2*(nz+1)
       xcor2(i) = (i-2)*dz-1.
    end do

    do i=1,2*(nz+1)
       if(xcor2(i) <= -epsilon) then
          rmollif(i) = 0.
       else if(xcor2(i) >= epsilon) then
          rmollif(i) = 0.
       else 
          tmp=1./((xcor2(i)/epsilon)**2-1)
          rmollif(i) = rinteg*exp(tmp)/epsilon
       endif
    end do
  end subroutine kernel

  subroutine initial_value()
    integer::i,j
    real pi, x, phi
    real,dimension(1:3)::u

    pi = 4.*atan(1.)
    phi = pi/8.
    do i=ljob,rjob
       do j=1,tri%n_node
          u = tri%node(j)%x
          f(i)%value(j) = 1./(4.*pi)+1./1.6e1*sqrt(5./pi)*&
               (3.*(u(2)*sin(phi)+u(3)*cos(phi))*(u(2)*sin(phi)+u(3)*cos(phi))-1.)       
       end do
       x = h1_fem_function_l1_norm(f(i))
       f(i)%value = (1./x)*f(i)%value
    end do

    v = 0.
    do i=0,nz
       zcord(i) = i*dz
       vg(i) = (v(i+1)-v(i))/dz
    end do

    call kernel()
  end subroutine initial_value

end module sphere_module

!
! end of file
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!
