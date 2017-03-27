!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! triangulation of simplex in 2d
!

module triangulation_module

	implicit none
	
	type, public::node
		integer::id
		real, dimension(1:3)::x
	end type node
	
	type, public::face
		integer::id
		integer, dimension(1:2)::node
		integer, dimension(1:2)::neighbour
	end type face
	
	type, public::element
		integer, dimension(1:3)::node
		integer, dimension(1:3)::neighbour
		integer, dimension(1:3)::face
	end type element
	
	type, public::edge
		integer, dimension(1:2)::vertex
	end type edge
	
	type, public::triangulation
		integer::n_node
		integer::n_face
		integer::n_element
		type(node), dimension(:), pointer::node
		type(face), dimension(:), pointer::face
		type(element), dimension(:), pointer::element
	end type triangulation
	
	contains

	subroutine triangulation_destroy(tri)
		type(triangulation), intent(inout)::tri
		
		if (associated(tri%node)) then
			deallocate(tri%node)
		end if
		if (associated(tri%face)) then
			deallocate(tri%face)
		end if
		if (associated(tri%element)) then
			deallocate(tri%element)
		end if
	end subroutine triangulation_destroy

	subroutine triangulation_read_mesh(tri, filename)
		type(triangulation), intent(inout)::tri
		character(*), intent(in)::filename
		
		integer::i, j, k, fid
		
		call triangulation_destroy(tri)
		fid = 10
		
		! read node data
		print *, "in triangulation_read_mesh: Reading node data ... ..."
		open(unit=fid, file=filename//".n", action="read")
		read(unit=fid, fmt=*) tri%n_node, tri%n_element, tri%n_face
		allocate(tri%node(1:tri%n_node))
		allocate(tri%face(1:tri%n_face))
		allocate(tri%element(1:tri%n_element))
		do i=1,tri%n_node
			read(unit=fid, fmt=*) j, tri%node(i)%x
		end do
		close(unit=fid)
		
		! read face data
		print *, "in triangulation_read_mesh: Reading face data ... ..."
		open(unit=fid, file=filename//".s", action="read")
		read(unit=fid, fmt=*) i
		if (i /= tri%n_face) then
			print *, "in triangulation_read_mesh: number of face not coincide"
			read *
		end if
		do i=1,tri%n_face
			read(unit=fid, fmt=*) j, tri%face(i)%node, tri%face(i)%neighbour
			tri%face(i)%node = tri%face(i)%node + 1
			tri%face(i)%neighbour = tri%face(i)%neighbour + 1
		end do
		close(unit=fid)
		
		! read element data
		print *, "in triangulation_read_mesh: Reading element data ... ..."
		open(unit=fid, file=filename//".e", action="read")
		read(unit=fid, fmt=*) i, j, k
		if (i /= tri%n_element) then
			print *, "in triangulation_read_mesh: number of element not coincide"
			read *
		end if
		if (j /= tri%n_node) then
			print *, "in triangulation_read_mesh: number of node not coincide"
			read *
		end if
		if (k /= tri%n_face) then
			print *, "in triangulation_read_mesh: number of face not coincide"
			read *
		end if
		do i=1,tri%n_element
			read(unit=fid, fmt=*) j, tri%element(i)%node, tri%element(i)%neighbour, &
				tri%element(i)%face
			tri%element(i)%node = tri%element(i)%node + 1
			tri%element(i)%neighbour = tri%element(i)%neighbour + 1
			tri%element(i)%face = tri%element(i)%face + 1
		end do
		close(unit=fid)
		
		print *, "OK!", tri%n_node, "nodes,", tri%n_face," faces,", tri%n_element, "elements."
	end subroutine triangulation_read_mesh
	
	real function triangulation_get_element_area(tri, el)
		type(triangulation), intent(in)::tri
		integer, intent(in)::el
		
		real, dimension(1:3)::v
		
		call cross_product(tri%node(tri%element(el)%node(2))%x-tri%node(tri%element(el)%node(1))%x, &
			tri%node(tri%element(el)%node(3))%x-tri%node(tri%element(el)%node(1))%x, v)
		triangulation_get_element_area = 0.5e0*sqrt(dot_product(v, v))!modified by LuoChong, add sqrt
	end function triangulation_get_element_area

	real function det3(x1, x2, x3)
		real, dimension(1:3), intent(in)::x1, x2, x3
		
		det3 = x1(1)*(x2(2)*x3(3) - x2(3)*x3(2)) &
			+ x1(2)*(x2(3)*x3(1) - x2(1)*x3(3)) &
			+ x1(3)*(x2(1)*x3(2) - x2(2)*x3(1))
	end function det3

	subroutine cross_product(x1, x2, x3)
		real, dimension(1:3), intent(in)::x1, x2
		real, dimension(1:3), intent(out)::x3

		x3(1) = x1(2)*x2(3) - x1(3)*x2(2)
		x3(2) = x1(3)*x2(1) - x1(1)*x2(3)
		x3(3) = x1(1)*x2(2) - x1(2)*x2(1)
	end subroutine cross_product

end module triangulation_module

!
! end of file
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

