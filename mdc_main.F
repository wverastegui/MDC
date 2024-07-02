      module m_mdc_main
     
      private
      public :: mdc_main

      CONTAINS

      subroutine mdc_main(istep)

C-----Load modules

      use precision, only      : dp        
      use siesta_options
      use sparse_matrices, only: maxnh, listh, listhptr, numh
      use sparse_matrices, only: H, S
      use atomlist, only       : no_u
      use files, only          : slabel
      use files, only          : filesOut_t   ! derived type for output file names
      use files, only          : label_length
      use m_spin, only         : nspin
      use siesta_geom, only    : na_u, isa
      use parallel, only       : Node
      use m_steps
      use fdf

      implicit none

C-----Define variables

      integer, intent(in)                 :: istep
      character(len=label_length+4), save :: fname
      character(len=label_length+4)       :: paste
      character(len=90 )                  :: err_crit
      integer,                       save :: unitmdc
      real(dp), save                      :: tunncond
      logical,  save                      :: writmdc
      external          io_assign, io_close, paste
      integer                             :: mdstep
C**************************************************************************
C    Call MDC routine
C**************************************************************************
          
          call tunn_coupl(na_u,isa,no_u,nspin,maxnh,numh,listhptr,listh,
     &                    H,S,tunncond,writmdc,err_crit,istep)

C**************************************************************************
          fname = paste(slabel,'.MDC')

          call io_assign(unitmdc)

          ! start writing MDC file
          open (unitmdc, file=fname, form = 'formatted', 
     .          position='append', status='unknown')

          if (Node .eq. 0) then
              ! write the column start (only for istep 1)
              if(istep .eq. inicoor) then
                 write(unitmdc,'(6a)') 
     .                 '# Step  ','   Splitting^2'
              endif
              ! write the tunneling coupling               
              if(writmdc .eq. .true.) then
                 write(unitmdc,'(i6,2x,E16.8)')
     .                 istep, tunncond
              else
                 write(unitmdc,'(a,i6,2x,a)')'#',istep,
     .                 err_crit
              endif          
          endif 

          call io_close(unitmdc)

      END subroutine mdc_main
      END module m_mdc_main
      
