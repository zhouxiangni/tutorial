!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! sphere.f90 : main engine and calculate right hand side
!


module sphere_module
  use sparse_matrix_module
  use amg_solver_module
  use triangulation_module
  use h1_fem_space_module

  implicit none

  type(triangulation)::tri
  type(h1_fem_space)::fem_space
  type(h1_fem_function)::f
  real::time
  real::delta_t
  real::bigU 
  real::v_z
  real::De   !Deborah number 

contains

  subroutine run(filename)
    character(*), intent(in)::filename

    integer::fid = 10
    type(sparse_matrix)::mass_matrix
    type(sparse_matrix)::stiff_matrix
    type(sparse_matrix)::M
    type(amg_solver)::solver
    real, dimension(:), allocatable::rhs
    real, dimension(1:3,1:3)::uu
    real, dimension(1:3,1:3,1:3,1:3)::uuuu 
    real, dimension(1:3,1:3)::tau
    real, dimension(1:3)::w
    real, dimension(1:102)::work
    real::f_l1_norm, error
    real::v_z, S_y, S_z, S
    integer:: i,k,l,p,q,iter,ierr

    De = 1.0e0
    v_z = 3.0e0
    bigU = 6.0e0*1.5e0
    delta_t = 1.0e-2

    call triangulation_read_mesh(tri, filename)
    call h1_fem_space_initialize(fem_space, tri)
    call h1_fem_function_initialize(f, fem_space)
    call initial_value()
    call h1_fem_space_get_mass_matrix(fem_space, mass_matrix)
    call h1_fem_space_get_stiff_matrix(fem_space, stiff_matrix)
    call sparse_matrix_initialize(M, fem_space%spH)
    M%matrix_entry = (1.0e0/delta_t)*mass_matrix%matrix_entry + (1.0e0/De)*stiff_matrix%matrix_entry
    call amg_solver_initialize(solver, M)
    allocate(rhs(1:tri%n_node))

    open(unit=fid+1, file="U6-G3-point.plt", action="write")
    write(unit=fid+1, fmt=*) "ZONE"
    close(unit=fid+1)

    time = 0.0e0
    do iter=1,200000
       call get_right_hand_side(rhs, uu, uuuu, v_z)
       call sparse_matrix_vector_mult_ps(mass_matrix, f%value, 1.0e0/delta_t, rhs)
       call amg_solver_solve(solver, f%value, rhs)
       
       !f_l1_norm = h1_fem_function_l1_norm(f)
       !f%value = (1.0e0/f_l1_norm)*f%value
       !print 1000, f_l1_norm
1000   format('|f|=',E25.19)
       open(unit=fid, file="solution-G3.dat", action="write")
       write(unit=fid, fmt=*) f%value
       close(unit=fid)

       tau=0.0e0
       do p=1,3
          do q=1,3
             tau(p,q)=3.0e0*uu(p,q)+0.5e0*De*v_z*uuuu(p,q,3,1)
             do k=1,3
                do l=1,3
                   tau(p,q)=tau(p,q)+2.0e0*bigU*uu(k,l)*uuuu(p,q,k,l)
                end do
             end do
             do k=1,3
                tau(p,q)=tau(p,q)-2.0e0*bigU*uu(p,k)*uu(q,k)
             end do
          end do
       end do
       do p=1,3
          tau(p,p) = tau(p,p)-1.0e0
       end do
       do k=1,3
          uu(k,k) = uu(k,k)-1.0e0/3.0e0
       end do
       S = 0.0e0
       do k=1,3
          do l=1,3
             S = S+uu(k,l)*uu(l,k)
          end do
       end do
       S = sqrt(1.5e0*S)
       call dsyev('V','U',3,uu,3,w,work,102,ierr)
       uu(:,3)=uu(:,3)/sqrt(dot_product(uu(:,3),uu(:,3)))
       !S_y=acos(abs(uu(2,3)))
       !S_z=acos(abs(uu(3,3)))
       
       open(unit=fid+1, file="U6-G3-point.plt", position='append', action="write")
       write(unit=fid+1, fmt='(F8.3,2X,3(E19.10,2X))') time,uu(1,3),uu(2,3),uu(3,3)!S,S_y,S_z!tau(1,3),tau(3,3)-tau(1,1),tau(2,2)-tau(3,3) 
       close(unit=fid+1)
       time = time + delta_t
    end do

    call sparse_matrix_destroy(mass_matrix)
    call sparse_matrix_destroy(stiff_matrix)
    call sparse_matrix_destroy(M)
    call amg_solver_destroy(solver)
    call h1_fem_space_destroy(fem_space)
    call h1_fem_function_destroy(f)
    call triangulation_destroy(tri)
    deallocate(rhs)
  end subroutine run

   subroutine get_right_hand_side(rhs, uu, uuuu, v_z)
    real, dimension(:), intent(out)::rhs
    real, dimension(1:3,1:3), intent(out)::uu
    real, dimension(1:3,1:3,1:3,1:3), intent(out)::uuuu
    real, intent(in)::v_z

    integer::i, j, k, l, p, q
    real, dimension(1:3,1:fem_space%n_quadrature_point)::quadrature_point
    real, dimension(1:fem_space%n_quadrature_point)::Jxw
    real, dimension(1:fem_space%n_quadrature_point)::f_on_quad_point
    real, dimension(1:3,1:3)::basis_gradient
    real, dimension(1:3)::grad_U, temp
    real::r
    real, dimension(1:3)::cross_x_Kx

    uu = 0.0e0
    uuuu = 0.0e0
    do i=1,tri%n_element
       call h1_fem_space_get_element_info(fem_space, i, quadrature_point, Jxw)
       call h1_fem_function_on_quad_point(f, i, f_on_quad_point)
       do j=1,fem_space%n_quadrature_point
          r = sqrt(dot_product(quadrature_point(:,j),quadrature_point(:,j)))
          do k=1,3
             do l=1,3
                uu(k,l) = uu(k,l) + Jxw(j)* &
                     (quadrature_point(k,j)/r)*(quadrature_point(l,j)/r)* &
                     f_on_quad_point(j)
                do p=1,3
                   do q=p,3
                      uuuu(p,q,k,l)=uuuu(p,q,k,l)+Jxw(j)* &
                           (quadrature_point(k,j)/r)*(quadrature_point(l,j)/r)* &
                           (quadrature_point(p,j)/r)*(quadrature_point(q,j)/r)* &
                           f_on_quad_point(j)
                   end do
                   do q=1,p-1
                      uuuu(p,q,k,l)=uuuu(q,p,k,l)
                   end do
                end do
             end do
          end do
       end do
    end do

    rhs = 0.0e0
    do i=1,tri%n_element
       call h1_fem_space_get_element_info(fem_space, i, quadrature_point, Jxw)
       call h1_fem_function_on_quad_point(f, i, f_on_quad_point)
       do j=1,fem_space%n_quadrature_point
          r = sqrt(dot_product(quadrature_point(:,j),quadrature_point(:,j)))
          temp = matmul((1.0e0/r)*quadrature_point(:,j), uu)
          call cross_product((-2.0e0/r)*quadrature_point(:,j), temp, grad_U)
          call h1_fem_space_basis_gradient(fem_space, i, quadrature_point(:,j), basis_gradient)
          cross_x_Kx(1) = 0.0e0
          cross_x_Kx(2) = v_z*quadrature_point(3,j)*quadrature_point(3,j)
          cross_x_Kx(3) = -v_z*quadrature_point(3,j)*quadrature_point(2,j)
          do k=1,3
             l = tri%element(i)%node(k)
             rhs(l) = rhs(l) &
                  -1.0e0/De*bigU*Jxw(j)*f_on_quad_point(j)* &
                  dot_product(grad_U,basis_gradient(k,:)) &
                  +Jxw(j)*f_on_quad_point(j)* &
                  dot_product(cross_x_Kx,basis_gradient(k,:)) 
          end do
       end do
    end do
  end subroutine get_right_hand_side


  subroutine initial_value()

    integer::i
    real::pi, x, phi
    real,dimension(1:3)::u

    pi = 4.0e0*atan(1.0e0)
    phi = pi/8.0e0
    do i=1,tri%n_node
       u = tri%node(i)%x
       f%value(i) = 1.0e0/(4.0e0*pi)+1.0e0/1.6e1*sqrt(5.0e0/pi)* &
          (3.0e0*(u(2)*sin(phi)+u(3)*cos(phi))*(u(2)*sin(phi)+u(3)*cos(phi))-1.0e0)
    end do
    !f%value=1.0e0/(4.0e0*pi)
    x = h1_fem_function_l1_norm(f)
    f%value = (1.0e0/x)*f%value
  end subroutine initial_value

end module sphere_module

!
! end of file
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

