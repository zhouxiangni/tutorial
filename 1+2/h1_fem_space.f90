!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! h1_fem_space.f90 : piecewise linear finite element space on sphere
!

module h1_fem_space_module
  use sparse_matrix_module
  use triangulation_module

  implicit none

  type, public::h1_fem_space
     type(triangulation), pointer::triangulation=>NULL()
     integer::n_quadrature_point
     real, dimension(:,:), pointer::quadrature_point=>NULL()
     real, dimension(:), pointer::weight=>NULL()
     type(sparsity_pattern)::spH
  end type h1_fem_space

  type, public::h1_fem_function
     type(h1_fem_space), pointer::fem_space=>NULL()
     real, dimension(:), pointer::value=>NULL()
  end type h1_fem_function

contains

  subroutine h1_fem_space_basis_function(sp, el, point, value)
    type(h1_fem_space), intent(in)::sp
    integer, intent(in)::el
    real, dimension(1:3), intent(in)::point
    real, dimension(1:3), intent(out)::value

    integer, dimension(1:3)::i
    real::D

    i = sp%triangulation%element(el)%node
    D = det3(sp%triangulation%node(i(1))%x, sp%triangulation%node(i(2))%x, sp%triangulation%node(i(3))%x)
    value(1) = det3(point, sp%triangulation%node(i(2))%x, sp%triangulation%node(i(3))%x)/D
    value(2) = det3(sp%triangulation%node(i(1))%x, point, sp%triangulation%node(i(3))%x)/D
    value(3) = det3(sp%triangulation%node(i(1))%x, sp%triangulation%node(i(2))%x, point)/D
  end subroutine h1_fem_space_basis_function

  subroutine h1_fem_space_basis_gradient(sp, el, point, value)
    type(h1_fem_space), intent(in)::sp
    integer, intent(in)::el
    real, dimension(1:3), intent(in)::point
    real, dimension(1:3,1:3), intent(out)::value

    integer, dimension(1:3)::i
    integer::j
    real:: D
    real, dimension(1:3)::l
    real, dimension(1:3)::v

    i = sp%triangulation%element(el)%node
    call cross_product(sp%triangulation%node(i(2))%x, sp%triangulation%node(i(3))%x, value(1,:))
    call cross_product(sp%triangulation%node(i(3))%x, sp%triangulation%node(i(1))%x, value(2,:))
    call cross_product(sp%triangulation%node(i(1))%x, sp%triangulation%node(i(2))%x, value(3,:))
    D = det3(sp%triangulation%node(i(1))%x, sp%triangulation%node(i(2))%x, sp%triangulation%node(i(3))%x)
    l(1) = dot_product(point,value(1,:))/D
    l(2) = dot_product(point,value(2,:))/D
    l(3) = dot_product(point,value(3,:))/D
    call cross_product(point,value(1,:),value(1,:))
    call cross_product(point,value(2,:),value(2,:))
    call cross_product(point,value(3,:),value(3,:))
    v = value(1,:)+value(2,:)+value(3,:)
    value(1,:) = (value(1,:)-l(1)*v)/D
    value(2,:) = (value(2,:)-l(2)*v)/D
    value(3,:) = (value(3,:)-l(3)*v)/D
    !print *, "debug: sumRlambda=",  value(1,:)+ value(2,:) + value(3,:)
  end subroutine h1_fem_space_basis_gradient

  subroutine h1_fem_space_initialize(sp, tri)
    type(h1_fem_space), intent(inout)::sp
    type(triangulation), target, intent(in)::tri

    integer::i, j, k, l, m, n_max_coupling_node
    integer, dimension(1:tri%n_node)::n_coupling_node

    call h1_fem_space_destroy(sp)
    sp%triangulation => tri
    
    !sp%n_quadrature_point = 10

    !allocate(sp%quadrature_point(1:3,1:sp%n_quadrature_point))
    !allocate(sp%weight(1:sp%n_quadrature_point))
    !sp%quadrature_point(:,1) = (/1.0e0, 0.0e0, 0.0e0/)
    !sp%quadrature_point(:,2) = (/0.0e0, 1.0e0, 0.0e0/)
    !sp%quadrature_point(:,3) = (/0.0e0, 0.0e0, 1.0e0/)
    !sp%quadrature_point(:,4) = (/0.0e0, 1.0e0/3.0e0, 2.0e0/3.0e0/)
    !sp%quadrature_point(:,5) = (/1.0e0/3.0e0, 0.0e0, 2.0e0/3.0e0/)
    !sp%quadrature_point(:,6) = (/1.0e0/3.0e0, 2.0e0/3.0e0, 0.0e0/)
    !sp%quadrature_point(:,7) = (/0.0e0, 2.0e0/3.0e0, 1.0e0/3.0e0/)
    !sp%quadrature_point(:,8) = (/2.0e0/3.0e0, 0.0e0, 1.0e0/3.0e0/)
    !sp%quadrature_point(:,9) = (/2.0e0/3.0e0, 1.0e0/3.0e0, 0.0e0/)
    !sp%quadrature_point(:,10) = (/1.0e0/3.0e0, 1.0e0/3.0e0, 1.0e0/3.0e0/)
    !sp%weight = (/1.0e0/3.0e1, 1.0e0/3.0e1, 1.0e0/3.0e1, &
    !     3.0e0/4.0e1, 3.0e0/4.0e1, 3.0e0/4.0e1, &
    !     3.0e0/4.0e1, 3.0e0/4.0e1, 3.0e0/4.0e1, &
    !     2.7e1/6.0e1/)

    sp%n_quadrature_point = 7

    allocate(sp%quadrature_point(1:3,1:sp%n_quadrature_point))
    allocate(sp%weight(1:sp%n_quadrature_point))
    sp%quadrature_point(:,1) = (/1.0e0, 0.0e0, 0.0e0/)
    sp%quadrature_point(:,2) = (/0.0e0, 1.0e0, 0.0e0/)
    sp%quadrature_point(:,3) = (/0.0e0, 0.0e0, 1.0e0/)
    sp%quadrature_point(:,4) = (/0.0e0, 0.5e0, 0.5e0/)
    sp%quadrature_point(:,5) = (/0.5e0, 0.0e0, 0.5e0/)
    sp%quadrature_point(:,6) = (/0.5e0, 0.5e0, 0.0e0/)
    sp%quadrature_point(:,7) = (/1.0e0/3.0e0, 1.0e0/3.0e0, 1.0e0-2.0e0/3.0e0/)
    sp%weight = (/3.0e0/6.0e1, 3.0e0/6.0e1, 3.0e0/6.0e1, &
       8.0e0/6.0e1, 8.0e0/6.0e1, 8.0e0/6.0e1, &
       2.7e1/6.0e1/)

    n_coupling_node = 0
    do i=1,tri%n_element
       n_coupling_node(tri%element(i)%node) = n_coupling_node(tri%element(i)%node)+3
    end do
    n_max_coupling_node = maxval(n_coupling_node)
    call sparsity_pattern_initialize(sp%spH, tri%n_node, tri%n_node, n_max_coupling_node)
    do i=1,tri%n_element
       do j=1,3
          k = tri%element(i)%node(j)
          call sparsity_pattern_add_entries(sp%spH, (/k, k, k/), tri%element(i)%node)
       end do
    end do
    call sparsity_pattern_compress(sp%spH)
  end subroutine h1_fem_space_initialize

  subroutine h1_fem_space_destroy(sp)
    type(h1_fem_space), intent(inout)::sp

    if (associated(sp%quadrature_point)) then
       deallocate(sp%quadrature_point)
    end if
    if (associated(sp%weight)) then
       deallocate(sp%weight)
    end if
    call sparsity_pattern_destroy(sp%spH)
  end subroutine h1_fem_space_destroy

  subroutine h1_fem_space_get_element_info(sp, el, quadrature_point, Jxw)
    type(h1_fem_space), intent(in)::sp
    integer, intent(in)::el
    real, dimension(1:3,1:sp%n_quadrature_point), intent(out)::quadrature_point
    real, dimension(1:sp%n_quadrature_point), intent(out)::Jxw

    integer::i, j

    do i=1,sp%n_quadrature_point
       quadrature_point(:,i) = 0.0e0
       do j=1,3
          quadrature_point(:,i) = quadrature_point(:,i) + sp%quadrature_point(j,i)* &
               sp%triangulation%node(sp%triangulation%element(el)%node(j))%x
       end do
    end do
    do i=1,sp%n_quadrature_point
       Jxw(i) = 0.5e0*sp%weight(i)*jacobian_determinant(sp%triangulation, el, quadrature_point(:,i))
    end do
  end subroutine h1_fem_space_get_element_info

  subroutine h1_fem_space_get_mass_matrix(sp, mat)
    type(h1_fem_space), intent(in)::sp
    type(sparse_matrix), intent(inout)::mat

    integer::i, j, k, j1, k1
    real, dimension(1:3, 1:sp%n_quadrature_point)::quadrature_point
    real, dimension(1:sp%n_quadrature_point)::Jxw
    real, dimension(1:3)::basis_function
    real, dimension(1:3,1:3)::element_matrix

    call sparse_matrix_initialize(mat, sp%spH)
    do i=1,sp%triangulation%n_element
       element_matrix = 0.0e0
       call h1_fem_space_get_element_info(sp, i, quadrature_point, Jxw)
       do j=1,sp%n_quadrature_point
          call h1_fem_space_basis_function(sp, i, quadrature_point(:,j), basis_function)
          do j1=1,3
             do k1=1,3
                element_matrix(j1, k1) = element_matrix(j1, k1) + Jxw(j)* &
                     basis_function(j1)*basis_function(k1)
             end do
          end do
       end do
       do j=1,3
          j1 = sp%triangulation%element(i)%node(j)
          do k=1,3
             k1 = sp%triangulation%element(i)%node(k)
             call sparse_matrix_add_entry(mat, j1, k1, element_matrix(j,k))
          end do
       end do
    end do
  end subroutine h1_fem_space_get_mass_matrix

  subroutine h1_fem_space_get_stiff_matrix(sp, mat)
    type(h1_fem_space), intent(in)::sp
    type(sparse_matrix), intent(inout)::mat

    integer::i, j, k, j1, k1
    real, dimension(1:3, 1:sp%n_quadrature_point)::quadrature_point
    real, dimension(1:sp%n_quadrature_point)::Jxw
    real, dimension(1:3, 1:3)::basis_gradient
    real, dimension(1:3, 1:3)::element_matrix

    call sparse_matrix_initialize(mat, sp%spH)
    do i=1,sp%triangulation%n_element
       element_matrix = 0.0e0
       call h1_fem_space_get_element_info(sp, i, quadrature_point, Jxw)
       do j=1,sp%n_quadrature_point
          call h1_fem_space_basis_gradient(sp, i, quadrature_point(:,j), basis_gradient)
          do j1=1,3
             do k1=1,3
                element_matrix(j1, k1) = element_matrix(j1, k1) + Jxw(j)* &
                     dot_product(basis_gradient(j1,:), basis_gradient(k1,:))
             end do
          end do
       end do
       do j=1,3
          j1 = sp%triangulation%element(i)%node(j)
          do k=1,3
             k1 = sp%triangulation%element(i)%node(k)
             call sparse_matrix_add_entry(mat, j1, k1, element_matrix(j,k))
          end do
       end do
    end do
  end subroutine h1_fem_space_get_stiff_matrix

  subroutine h1_fem_function_initialize(f, sp)
    type(h1_fem_function), intent(inout)::f
    type(h1_fem_space), target, intent(in)::sp

    if (.not. associated(f%fem_space, sp)) then
    	call h1_fem_function_destroy(f)
    	f%fem_space => sp
    	allocate(f%value(1:sp%triangulation%n_node))
    end if
    f%value = 0.0e0
  end subroutine h1_fem_function_initialize

  subroutine h1_fem_function_destroy(f)
    type(h1_fem_function), intent(inout)::f

    if (associated(f%value)) then
       deallocate(f%value)
    end if
    nullify(f%fem_space)
  end subroutine h1_fem_function_destroy
  
  subroutine h1_fem_function_on_quad_point(f, el, value)
    type(h1_fem_function), intent(in)::f
    integer, intent(in)::el
    real, dimension(1:f%fem_space%n_quadrature_point), intent(out)::value

    integer::i, j, k
    real, dimension(1:3,1:f%fem_space%n_quadrature_point)::quadrature_point
    real, dimension(1:f%fem_space%n_quadrature_point)::Jxw
    real, dimension(1:3)::basis_function

    value = 0.0e0
    call h1_fem_space_get_element_info(f%fem_space, el, quadrature_point, Jxw)
    do i=1,f%fem_space%n_quadrature_point
       call h1_fem_space_basis_function(f%fem_space, el, quadrature_point(:,i), basis_function)
       do j=1,3
          k = f%fem_space%triangulation%element(el)%node(j)
          value(i) = value(i) + f%value(k)*basis_function(j)
       end do
    end do
  end subroutine h1_fem_function_on_quad_point

  function h1_fem_function_l1_norm(f) result(l1_norm)
    type(h1_fem_function), intent(in)::f
    real::l1_norm

    integer::i
    real, dimension(1:3,1:f%fem_space%n_quadrature_point)::quadrature_point
    real, dimension(1:f%fem_space%n_quadrature_point)::Jxw
    real, dimension(1:f%fem_space%n_quadrature_point)::f_on_quad_point

    l1_norm = 0.0e0
    do i=1,f%fem_space%triangulation%n_element
       call h1_fem_space_get_element_info(f%fem_space, i, quadrature_point, Jxw)
       call h1_fem_function_on_quad_point(f, i, f_on_quad_point)
       l1_norm = l1_norm + dot_product(Jxw, abs(f_on_quad_point))
    end do
  end function h1_fem_function_l1_norm

  function h1_fem_function_l2_norm(f) result(l2_norm)
    type(h1_fem_function), intent(in)::f
    real::l2_norm

    integer::i
    real, dimension(1:3,1:f%fem_space%n_quadrature_point)::quadrature_point
    real, dimension(1:f%fem_space%n_quadrature_point)::Jxw
    real, dimension(1:f%fem_space%n_quadrature_point)::f_on_quad_point

    l2_norm = 0.0e0
    do i=1,f%fem_space%triangulation%n_element
       call h1_fem_space_get_element_info(f%fem_space, i, quadrature_point, Jxw)
       call h1_fem_function_on_quad_point(f, i, f_on_quad_point)
       l2_norm = l2_norm + dot_product(Jxw, f_on_quad_point**2)
    end do
    l2_norm = sqrt(l2_norm)
  end function h1_fem_function_l2_norm
  
  real function jacobian_determinant(tri, el, point)
    type(triangulation), intent(in)::tri
    integer, intent(in)::el
    real, dimension(1:3), intent(in)::point

    real::r, d, sin_theta, cos_theta
    real, dimension(1:3)::v1, v2

    r = sqrt(dot_product(point, point))
    call cross_product(tri%node(tri%element(el)%node(2))%x-tri%node(tri%element(el)%node(1))%x, &
         tri%node(tri%element(el)%node(3))%x-tri%node(tri%element(el)%node(1))%x, v1)
    d = sqrt(dot_product(v1, v1))
    call cross_product(point, v1, v2)
    sin_theta = sqrt(dot_product(v2, v2))/(r*d)
    cos_theta = sqrt(1.0e0 - sin_theta*sin_theta)
    jacobian_determinant = d/(r*r*cos_theta)
  end function jacobian_determinant

end module h1_fem_space_module
	
!
! end of file
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!

