! ***********************************************************************
!
!   Copyright (C) 2010-2019  The MESA Team
!
!   this file is part of mesa.
!
!   mesa is free software; you can redistribute it and/or modify
!   it under the terms of the gnu general library public license as published
!   by the free software foundation; either version 2 of the license, or
!   (at your option) any later version.
!
!   mesa is distributed in the hope that it will be useful, 
!   but without any warranty; without even the implied warranty of
!   merchantability or fitness for a particular purpose.  see the
!   gnu library general public license for more details.
!
!   you should have received a copy of the gnu library general public license
!   along with this software; if not, write to the free software
!   foundation, inc., 59 temple place, suite 330, boston, ma 02111-1307 usa
!
! ***********************************************************************
 
      module run_star_extras

      use star_lib
      use star_def
      use const_def
      use math_lib
      use auto_diff
      
      implicit none
      
      contains

        subroutine axion_other_energy(id, ierr)
                use star_def
                use auto_diff
                use const_def, only: Rsun
                integer, intent(in) :: id
                integer, intent(out) :: ierr
                type (star_info), pointer :: s
                integer :: k

                integer :: nz
                integer :: i
                integer :: i_na23
                real(dp) :: g ! effective coupling
                real(dp) :: eps0, T0, mu0 ! baseline values
                real(dp), allocatable :: x_na23(:), T(:), dm(:), E(:) ! intermediate arrays

         call star_ptr(id, s, ierr)

         nz = s% nz

         allocate(x_na23(nz))
         allocate(T(nz))
         allocate(dm(nz))
         allocate(E(nz))

         ierr = 0
         i_na23 = 15 ! index of species

         eps0 = 8.6E+27 ! ergs / g s; baseline
         T0 = 5.106E+9 ! Kelvins; baseline
         mu0 = 1.5 ! unitless; chemical potential

         g = s% x_ctrl(99) ! unitless; effective coupling

         x_na23(1:nz) = s% xa(i_na23, 1:nz)
         T(1:nz) = s% T(1:nz)
         dm(1:nz) = s% dm(1:nz)
         E(1:nz) = eps0 * x_na23(1:nz) * g*g / ( exp(T0 / T(1:nz)) + 1.5 )

         if (ierr /= 0) return

         s% xtra(1) = sum(E(1:nz) * dm(1:nz)) ! total axion power

         write (*,*) "effective coupling:"
         write (*,*) g

         write (*,*) "total axion power (ergs/s):"
         write (*,*) s% xtra(1)

         if (s% x_logical_ctrl(99)) then ! if we want axions to affect evolution
                 write (*,*) "axion cooling enabled"
                 s% extra_heat(1:nz) %val = -E(1:nz)
                 ! note that extra_heat is type(auto_diff_real_star_order1) so includes partials.
         end if

        !do, i=1, nz
        !write (*,*) x_na23(i)
        !enddo

 end subroutine axion_other_energy
      
      subroutine extras_controls(id, ierr)
         integer, intent(in) :: id
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return

         s% other_energy => axion_other_energy
         
         s% extras_startup => extras_startup
         s% extras_check_model => extras_check_model
         s% extras_finish_step => extras_finish_step
         s% extras_after_evolve => extras_after_evolve
         s% how_many_extra_history_columns => how_many_extra_history_columns
         s% data_for_extra_history_columns => data_for_extra_history_columns
         s% how_many_extra_profile_columns => how_many_extra_profile_columns
         s% data_for_extra_profile_columns => data_for_extra_profile_columns  

         s% how_many_extra_history_header_items => how_many_extra_history_header_items
         s% data_for_extra_history_header_items => data_for_extra_history_header_items
         s% how_many_extra_profile_header_items => how_many_extra_profile_header_items
         s% data_for_extra_profile_header_items => data_for_extra_profile_header_items

      end subroutine extras_controls
      
      
      subroutine extras_startup(id, restart, ierr)
         integer, intent(in) :: id
         logical, intent(in) :: restart
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         include 'formats'
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
      end subroutine extras_startup
      
      
      subroutine extras_after_evolve(id, ierr)
         integer, intent(in) :: id
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         real(dp) :: dt
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
      end subroutine extras_after_evolve
      

      ! returns either keep_going, retry, or terminate.
      integer function extras_check_model(id)
         integer, intent(in) :: id
         integer :: ierr
         
         real(dp), parameter :: Blocker_scaling_factor_after_TP = 5d0
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         
         include 'formats'
         
         extras_check_model = keep_going

         ! retry if axion emission changes too much
         if ( abs( s% xtra(1) * s% dt ) > s% x_ctrl(98) * s% x_ctrl(99) ) then
                 extras_check_model = retry
                 write(*, *) 'retry: axion emitted exceeds x_ctrl(98) x_ctrl(99) ergs'
                 write(*, *) s% dt
                 return
         end if

         
         !if (abs(s% Blocker_scaling_factor - Blocker_scaling_factor_after_TP) < 1d-8) return
         
         !if (s% center_he4 < 1d-4 .and. &
         !      any(s% burn_he_conv_region(1:s% num_conv_boundaries)) .and. &
         !      s% he_core_mass - s% c_core_mass <= s% TP_he_shell_max) then
         !   !write(*,1) 'set Blocker_scaling_factor = Blocker_scaling_factor_after_TP', Blocker_scaling_factor_after_TP
         !   !s% Blocker_scaling_factor = Blocker_scaling_factor_after_TP
         !   write(*,*) '1st thermal pulse'
         !   extras_check_model = terminate
         !   s% termination_code = t_extras_finish_step
         !end if
         
      end function extras_check_model


      integer function how_many_extra_history_columns(id)
         integer, intent(in) :: id
         integer :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         how_many_extra_history_columns = 0
      end function how_many_extra_history_columns
      
      
      subroutine data_for_extra_history_columns(id, n, names, vals, ierr)
         integer, intent(in) :: id, n
         character (len=maxlen_history_column_name) :: names(n)
         real(dp) :: vals(n)
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
      end subroutine data_for_extra_history_columns

      
      integer function how_many_extra_profile_columns(id)
         use star_def, only: star_info
         integer, intent(in) :: id
         integer :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         how_many_extra_profile_columns = 0
      end function how_many_extra_profile_columns
      
      subroutine data_for_extra_profile_columns(id, n, nz, names, vals, ierr)
         use star_def, only: star_info, maxlen_profile_column_name
         use const_def, only: dp
         integer, intent(in) :: id, n, nz
         character (len=maxlen_profile_column_name) :: names(n)
         real(dp) :: vals(nz,n)
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         integer :: k
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return

         ! note: do NOT add the extra names to profile_columns.list
         ! the profile_columns.list is only for the built-in profile column options.
         ! it must not include the new column names you are adding here.

         ! here is an example for adding a profile column
         if (n /= 1) stop 'data_for_extra_profile_columns'
         names(1) = 'axion'
         do k = 1, nz
            vals(k,1) = s% extra_heat(k) %val
         end do

      end subroutine data_for_extra_profile_columns

      integer function how_many_extra_profile_header_items(id)
         integer, intent(in) :: id
         integer :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         how_many_extra_profile_header_items = 2
      end function how_many_extra_profile_header_items


      subroutine data_for_extra_profile_header_items(id, n, names, vals, ierr)
         integer, intent(in) :: id, n
         character (len=maxlen_profile_column_name) :: names(n)
         real(dp) :: vals(n)
         type(star_info), pointer :: s
         integer, intent(out) :: ierr
         ierr = 0
         call star_ptr(id,s,ierr)
         if(ierr/=0) return

         ! here is an example for adding an extra profile header item
         ! also set how_many_extra_profile_header_items
         names(1) = 'g_eff'
         vals(1) = s% x_ctrl(99)

         names(2) = 'lum_axion_surf'
         vals(2) = s% xtra(1)

      end subroutine data_for_extra_profile_header_items

      integer function how_many_extra_history_header_items(id)
         integer, intent(in) :: id
         integer :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         how_many_extra_history_header_items = 0
      end function how_many_extra_history_header_items


      subroutine data_for_extra_history_header_items(id, n, names, vals, ierr)
         integer, intent(in) :: id, n
         character (len=maxlen_history_column_name) :: names(n)
         real(dp) :: vals(n)
         type(star_info), pointer :: s
         integer, intent(out) :: ierr
         ierr = 0
         call star_ptr(id,s,ierr)
         if(ierr/=0) return

         ! here is an example for adding an extra history header item
         ! also set how_many_extra_history_header_items
         ! names(1) = 'mixing_length_alpha'
         ! vals(1) = s% mixing_length_alpha

      end subroutine data_for_extra_history_header_items
      

      ! returns either keep_going or terminate.
      integer function extras_finish_step(id)
         integer, intent(in) :: id
         integer :: ierr
         type (star_info), pointer :: s
         include 'formats'
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         extras_finish_step = keep_going
         if (s% x_integer_ctrl(1) == 2) then ! part2
            if (s% L(1)/Lsun < s% x_ctrl(1) .and. &
                s% Teff > s% x_ctrl(2) .and. &
                safe_log10(s% power_he_burn) < s% x_ctrl(3)) then
               write(*,1) 'L/Lsun < limit', s% L(1)/Lsun, s% x_ctrl(1)
               write(*,1) 'Teff > limit', s% Teff, s% x_ctrl(2)
               write(*,1) 'log_LHe < limit', &
                  safe_log10(s% power_he_burn), s% x_ctrl(3)
               extras_finish_step = terminate
            end if
         end if
      end function extras_finish_step
      
      

      end module run_star_extras
      
