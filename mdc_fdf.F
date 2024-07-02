      module m_mdc_fdf
      implicit none 
      private
      public :: mdc_fdf

      contains
      function dirname(number)
      integer    :: number
      character(len=60)  :: dirname               
      write (dirname, '(i6.4)')  number
      end function 


      subroutine mdc_fdf(istep)
C-----Load modules
      
      use precision, only      : dp        
      use siesta_options
      use files, only          : slabel
      use files, only          : filesOut_t   ! derived type for output file names
      use files, only          : label_length
      use siesta_geom, only    : na_u, isa, xa, ucell
      use parallel, only       : Node, IOnode
      use m_steps
      use fdf
      use atmfuncs, only       : nofis, symfio, cnfigfio, labelfis
      use atmfuncs, only       : atmpopfio, izofis 
      use atomlist, only       : iaorb, iphorb
      use parallelsubs 
      use basis_types, only    : nsp
      

      integer, intent(in)                 :: istep
      character(len=label_length+4), save :: fname,inpname
      character(len=label_length+4)       :: paste
      integer,                       save :: unitmdc
      logical,  save                      :: writmdc
      external          io_assign, io_close, paste
      CHARACTER(len=255) :: dir, workdir,comm,path_dir
      integer            :: i, j, ia, iua, na_ul, is
      integer            :: na_cx1, na_cx2, atm1, atm2
      real(dp), allocatable           :: xal(:,:)
      integer, allocatable            :: isal(:)
      logical                         :: ex 
      character(len=100)       :: tmp 
      integer, parameter       :: nk=1
      real(dp)                 :: factor

C**************************************************************************
      if (IOnode) then
      fname = 'STRUCT.fdf'
      inpname = 'RUN_S.fdf'
      
!     Create work directory for single point calculation
      call system('mkdir  MD-step-'// 
     &         adjustl(trim( dirname(istep) ) ) )   
      dir = '/MD-step-'//adjustl(trim(dirname(istep)))
      call getcwd(workdir)
      path_dir = trim(workdir)//trim(dir)


!     Write STRUCT.fdf of complex on workdir
      call io_assign(unitmdc)
      open (unitmdc, file=fname, form = 'formatted', 
     .      position='append', status='unknown',
     .      defaultfile=path_dir)
      
      atm1   = fdf_get('MDC.Adt1', na_cx1)
      atm2   = fdf_get('MDC.Adt2', na_cx2)
      na_ul  = 1 + (atm2 - atm1)  
      
      allocate (xal(3,na_ul))
      allocate (isal(na_ul))
      xal    = 0.0D0
      ia     = 0
      factor = 0.529177_dp  
      do iua = atm1,atm2
          ia = ia + 1
          do j  = 1,3
           xal(j,ia) = xa(j,iua)*factor
          enddo
      enddo
      isal(1:na_ul)  = isa(atm1:atm2)
      
      write(unitmdc,'(a,6x,i5)')'NumberOfAtoms:',na_ul
      write(unitmdc,'(a,4x,i5)')'NumberOfSpecies:',nsp 
      write(unitmdc,'(a)') ' '
      write(unitmdc,'(a,6x,a,6x,a)') 'LatticeConstant','1.0','Ang'
      write(unitmdc,'(a)') '%block LatticeVectors'
      write(unitmdc,'(3(3f12.6/))') 
     .            ((ucell(i,j)*factor,i=1,3),j=1,3)
      write(unitmdc,'(a)') '%endblock LatticeVectors'
      write(unitmdc,'(a)') ' '
      write(unitmdc,'(/a)') '%block ChemicalSpeciesLabel' 
        do is=1,nsp
           write(unitmdc,'(2(1x,i4),1x,a)')
     .          is,izofis(is),labelfis(is) 
        enddo
      write(unitmdc,'(a)') '%endblock ChemicalSpeciesLabel' 
      write(unitmdc,'(a)') ' '
      write(unitmdc,'(a)') 'AtomicCoordinatesFormat Ang'
      write(unitmdc,'(a)') '%block AtomicCoordinatesAndAtomicSpecies'
      do i = 1,na_ul
             write(unitmdc,'(3f10.6,4x,i5)') 
     .            (xal(j,i),j=1,3), isal(i)
      enddo
      write(unitmdc,'(a)') '%endblock AtomicCoordinatesAndAtomicSpecies'
      deallocate (xal)
      deallocate (isal)
      call io_close(unitmdc)

      
      call io_assign(unitmdc)
      open (unitmdc, file='tmp.txt' , form = 'formatted', 
     .           status='unknown')  
      write(unitmdc,'(a)') path_dir
      call io_close(unitmdc)

!     Execute bash script to copy files
      comm = './Files_setup.sh'
      call system(comm)
      endif

      END subroutine mdc_fdf
      END module m_mdc_fdf
      
      
