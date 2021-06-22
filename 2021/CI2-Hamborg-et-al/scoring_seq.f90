program tt
implicit none

integer, parameter :: NSEQ_MAX=100000
integer :: i=1,k,argcount,IARGC,NPOS,NSeq,output,identical,thr,start,fine
real,allocatable  :: J(:,:)
real,allocatable :: hf(:)
integer,allocatable :: seq(:)
character(len=200000) line
character(len=10000) line_a(NSEQ_MAX),line_r(NSEQ_MAX)
!real*8,allocatable  :: atom_pos(:,:,:,:)
real*8          :: r,Eh,Ej

character*50    :: wq_char,file_in1,file_in2,file_clean
!character*3     :: ResName(NRES_max,NCHAIN_max)

identical=0
output=0
thr=0
start=1
fine=10000

argcount = IARGC()

file_in1='XXX'
file_in2='XXX'

  if(argcount.lt.2)then
     write(6,*)'USAGE:'
     write(6,*)'scoring_seq.x -PAR potts_parameters_file -SEQ msa'
    stop 'ERROR: Inputs Missing'
  endif

  do i=1,argcount
     call getarg(i,wq_char)

     if (INDEX(wq_char,'-PAR').NE.0) then
        call getarg(i+1,wq_char)
        read(wq_char,*)file_in1
     endif
     if (INDEX(wq_char,'-SEQ').NE.0) then
        call getarg(i+1,wq_char)
        read(wq_char,*)file_in2
     endif

  enddo


write(file_clean, fmt='(a,a)') 'scores_',trim(file_in2)
write(*,*),'INPUT:  ',file_in2
write(*,*),'OUTPUT: ',file_clean


call read_file_in()


contains

subroutine read_file_in()

  implicit none
  integer :: ii, ll(NSEQ_MAX),ia,jj,id,h,a(10000)
  integer*8 tot,r(10000),thr2

  line_a(NSEQ_MAX)=''
  line_r(NSEQ_MAX)=''

  open(22,file=file_in1,status='old')
  open(unit=24,file=file_clean,status='unknown')
  thr2=0
  id=0
  ii=0
  do
     read(22,*,end=112) line

     ii=ii+1

  enddo
  112 continue

  rewind(22)
  NPOS=(ii-1)/21
  write(*,*),ii,NPOS
  allocate(J(ii-1,ii-1))
  allocate(hf(ii-1),seq(NPOS))

  jj=1
  do
     if(jj.lt.ii) then

        read(22,*,end=113) J(jj,1:ii-1)

        jj=jj+1

     else
        read(22,*,end=113) hf(1:ii-1)
     
     endif


  enddo
  113 continue
  open(23,file=file_in2,status='old')
     jj=0
     do
        read(23,*,end=114) seq(1:NPOS)
        jj=jj+1
        Eh=0.
        Ej=0.
        do i=1,NPOS
           Eh=Eh-hf((i-1)*21+seq(i))
           
           do k=i+1,NPOS
              Ej=Ej-J((i-1)*21+seq(i),(k-1)*21+seq(k))
           enddo
        enddo
        write(*,*),jj,Ej+Eh,Eh
     enddo
 
114  continue 
  close(22)
  close(23)
deallocate(J,hf,seq)

  end subroutine read_file_in

end program tt
