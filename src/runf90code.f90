
! Haplotype association sampler

! Fortran version, converted from R version

! John.Henshall@csiro.au
! Nov 2013

! Deliberate global variables to speed execution

subroutine runf90code()
!use variables

implicit none

            integer, allocatable :: seed(:)
            integer ::   j, n, un, istat, dt(8), pid, t(2), s
            integer(8) :: count, tms

integer :: nanis
integer :: nhap
integer :: nsamp

! AWG integer, dimension(:), allocatable :: id
character(len=10), dimension(:), allocatable :: id
integer, dimension(:), allocatable :: hapgeno,hapgenocand
integer, dimension(:), allocatable :: g1i,g2i,g,g1j,g2j

integer, dimension(:), allocatable :: h1,h2,ispair  ! haplotype alleles

character(len=10), dimension(:), allocatable :: pheno

integer, dimension(:), allocatable :: alleles
logical, dimension(:), allocatable :: wiped

real, dimension(:,:), allocatable :: likpen, Ppen, Hpen

real, dimension(:), allocatable :: randr

integer, dimension(:,:), allocatable :: hapgraph
integer, dimension(:), allocatable :: nhapgraph

! Logs
integer, dimension(:), allocatable :: depth_log
logical, dimension(:), allocatable :: accept_log
real, dimension(:), allocatable :: pi_i_log
real, dimension(:), allocatable :: pi_j_log
real, dimension(:), allocatable :: q_ij_log
real, dimension(:), allocatable :: q_ji_log
real, dimension(:), allocatable :: geno_hap_log
real, dimension(:), allocatable :: geno_ani_log


! from runit
integer :: i,depth,d,ali,l,jj,srep
integer :: old_geno, new_geno

real :: x,q_ij,q_ji,pi_i, pi_j, PP, PH, prob,lrat,xlim

logical :: accept

integer, dimension(:,:), allocatable :: local_hapgraph
integer, dimension(:), allocatable :: local_nhapgraph,nodes,noden,indi
logical, dimension(:), allocatable :: ssset
integer, dimension(:,:), allocatable :: tempi
real, dimension(:), allocatable :: likpen_vec
real, dimension(:), allocatable :: Ppen_vec,Hpen_vec
character(len=57) :: outfile




integer :: nhap1,nhap2


! First file parameters and run parameters 

open(unit=11,file="TmpData/forfortran.txt",action="read",status="unknown")
read (11,*) nsamp
read (11,*) nanis

! Allocate arrays
allocate(id(nanis))
allocate(pheno(nanis))
allocate(g1i(nanis),g2i(nanis),g1j(nanis),g2j(nanis),g(nanis))
allocate(h1(nanis),h2(nanis))
allocate(likpen(nanis,3),Ppen(nanis,3),Hpen(nanis,3))

do i = 1,nanis
  read (11,*) id(i)
  read (11,*) pheno(i)
  read (11,*) h1(i),h2(i)
  read (11,*) Ppen(i,:)
  read (11,*) Hpen(i,:)
  read (11,*) likpen(i,:)
enddo

close(11)

nhap1 = maxval(h1)
nhap2 = maxval(h2)
nhap = max(nhap1,nhap2)

allocate(hapgeno(nhap),hapgenocand(nhap),randr(nhap))

! Read in starting values

open(unit=11,file="TmpData/starting.txt",action="read",status="unknown")
read (11,*) nhap

do i = 1,nhap
  read (11,*) hapgeno(i)
enddo

close(11)

allocate(alleles(nhap),wiped(nhap),ispair(nanis))

allocate(hapgraph(nhap,nhap),nhapgraph(nhap))

! Logs
allocate(depth_log(nsamp))
allocate(accept_log(nsamp))
allocate(pi_i_log(nsamp))
allocate(pi_j_log(nsamp))
allocate(q_ij_log(nsamp))
allocate(q_ji_log(nsamp))
allocate(geno_hap_log(nhap))
allocate(geno_ani_log(nanis))

! print*,"All data read, ",nhap," haplotype alleles, ",nanis," animals"



!---------------------------
!    Prepare
!---------------------------



nhapgraph = 0
hapgraph = 0

! hapgraph is a full stored adjancy matrix, which is the hap if
! an allele pair occured in a diploid animal, zero otherwise.

hapgraph = 0

! First time through, for counting only
do j = 1,nanis
  hapgraph(h1(j),h2(j)) = 1
  hapgraph(h2(j),h1(j)) = 1
enddo

! Zero diagonals
do j = 1,nhap
  hapgraph(j,j) = 0
enddo

! counts

do i = 1,nhap
  nhapgraph(i) = sum(hapgraph(i,:))
enddo

! Redo, with col number this time
hapgraph = 0
do j = 1,nanis
  hapgraph(h1(j),h2(j)) = h2(j)
  hapgraph(h2(j),h1(j)) = h1(j)
enddo

! Zero diagonals
do j = 1,nhap
  hapgraph(j,j) = 0
enddo

j = min(nhap,20)

open(unit=13,file="Temp/hapgraph.txt",action="write",status="unknown")
do i = 1,j
  write(13,"(66I6)") nhapgraph(i),hapgraph(i,1:j)
enddo
close(13)


!-----------------------
!  init_random_seed
!-----------------------

          
            un = 44
            call random_seed(size = n)
            allocate(seed(n))
            ! First try if the OS provides a random number generator
            !open(newunit=un, file="/dev/urandom", access="stream", &
            !     form="unformatted", action="read", status="old", iostat=istat)
            !if (istat == 0) then
            !   read(un) seed
            !   close(un)
            !else
               ! Fallback to XOR:ing the current time and pid. The PID is
               ! useful in case one launches multiple instances of the same
               ! program in parallel.
               call system_clock(count)
               if (count /= 0) then
                  t = transfer(count, t)
               else
                  call date_and_time(values=dt)
                  tms = (dt(1) - 1970) * 365_8 * 24 * 60 * 60 * 1000 &
                       + dt(2) * 31_8 * 24 * 60 * 60 * 1000 &
                       + dt(3) * 24 * 60 * 60 * 60 * 1000 &
                       + dt(5) * 60 * 60 * 1000 &
                       + dt(6) * 60 * 1000 + dt(7) * 1000 &
                       + dt(8)
                  t = transfer(tms, t)
               end if
               s = ieor(t(1), t(2))
               pid = 1099279 ! Add a prime
               s = ieor(s, pid)
               if (n >= 3) then
                  seed(1) = t(1) + 36269
                  seed(2) = t(2) + 72551
                  seed(3) = pid
                  if (n > 3) then
                     seed(4:) = s + 37 * (/ (i, i = 0, n - 4) /)
                  end if
               else
                  seed = s + 37 * (/ (i, i = 0, n - 1 ) /)
               end if
            !end if
            call random_seed(put=seed)

!-------------------------------
!  runit
!-------------------------------


xlim = 15

allocate(local_hapgraph(nhap,nhap),local_nhapgraph(nhap),indi(nhap))
allocate(likpen_vec(nanis))
allocate(Ppen_vec(nanis))
allocate(Hpen_vec(nanis))

do i = 1,nhap
  indi(i) = i
enddo




  geno_ani_log = 0
  geno_hap_log = 0

  ! Assign genotypes
  g1i = hapgeno(h1)
  g2i = hapgeno(h2)
  g = g1i + g2i - 2

  ! Compute log-likelihood of initial sample

  ! likpen columns are 0, 1, 2 copies of H

  where (g == 0)
    likpen_vec = likpen(:,1)
  else where (g == 1)
    likpen_vec = likpen(:,2)
  else where (g == 2)
    likpen_vec = likpen(:,3)
  end where

  pi_i = sum(likpen_vec)


  do srep = 1,nsamp

    alleles = 0

    call random_number(x)
    i = 1
    j = nhap
    ali = int(x*(j+1-i))+i

      
    alleles(1) = ali

    ! Maybe not changing everything here, just some of the QTL genotypes?
    ! Sample a depth to wipe
    call random_number(x)
    i = 1
    j = nhap
    depth = int(x*(j+1-i))+i

    ! But force occasional shallow sample
    call random_number(x)
    if (x < 0.25) then
      call random_number(x)
      i = 1
      j = 2
      depth = int(x*(j+1-i))+i
    endif

    local_hapgraph = hapgraph
    local_nhapgraph = nhapgraph

    ! Zeroing will reduce counts for linkled haps
    where (local_hapgraph(:,ali) > 0)
      local_nhapgraph = local_nhapgraph - 1 
    end where

    ! Zero out column for first hap
    local_hapgraph(:,ali) = 0

    jj = min(nhap,20)

    ! Now get haps linked through heterozygous animals 
    do d = 1,depth-1

      ! Sample an already seen haplotype with linked nodes remaining
      noden = local_nhapgraph(alleles(1:d))

      if (sum(noden) == 0) then
        depth = d
        exit
      endif

      nodes = indi(alleles(1:d))
      nodes = pack(nodes,noden > 0)

      l = size(nodes)

      call random_number(x)
      i = 1
      j = l
      ali = int(x*(j+1-i))+i


      ! Sample a linked haplotype
      nodes = pack(indi,local_hapgraph(ali,:) > 0)

      l = size(nodes)

      call random_number(x)
      i = 1
      j = l
      ali = int(x*(j+1-i))+i


      ! Zeroing will reduce counts for linkled haps
      where (local_hapgraph(:,ali) > 0)
        local_nhapgraph = local_nhapgraph - 1 
      end where

      ! Zero out column 
      local_hapgraph(:,ali) = 0

      alleles(d+1) = ali

    enddo

    ! Wipe

    hapgenocand = hapgeno
    hapgenocand(alleles(1:depth)) = 0

    g1j = hapgenocand(h1)
    g2j = hapgenocand(h2)

    ! Resample, use same order as wiped.
    ! Accumulate qij and qji at the same time

    q_ij = 0
    q_ji = 0

    do d  = 1,depth

      ali = alleles(d)


      ! Ppen columns are other allele = P (1), H (2) or missing (0)
      ! Homozygs are already covered, they are in the missing column

      Ppen_vec = 0
      Hpen_vec = 0
      where (h1 == ali .and. g2j == 1)
        Ppen_vec = Ppen(:,1)
        Hpen_vec = Hpen(:,1)
      else where (h1 == ali .and. g2j == 2)
        Ppen_vec = Ppen(:,2)
        Hpen_vec = Hpen(:,2)
      else where (h1 == ali .and. g2j == 0)
        Ppen_vec = Ppen(:,3)
        Hpen_vec = Hpen(:,3)
      else where (h2 == ali .and. g1j == 1)
        Ppen_vec = Ppen(:,1)
        Hpen_vec = Hpen(:,1)
      else where (h2 == ali .and. g1j == 2)
        Ppen_vec = Ppen(:,2)
        Hpen_vec = Hpen(:,2)
      else where (h2 == ali .and. g1j == 0)
        Ppen_vec = Ppen(:,3)
        Hpen_vec = Hpen(:,3)
      end where

      PP = sum(Ppen_vec)
      PH = sum(Hpen_vec)

      ! Rescale
      x = min(PP,PH)
      PP = PP - x
      PH = PH - x
      
      if (PP > xlim) PP = xlim
      if (PH > xlim) PH = xlim

      PP = exp(PP)
      PH = exp(PH)

      prob = PP / (PP + PH) 

      ! Sample
      call random_number(x)
      if (x < prob) then        ! sample P
        hapgenocand(ali) = 1
        q_ij = q_ij + log(prob)
        where (h1 == ali) g1j = 1
        where (h2 == ali) g2j = 1
      else                      ! sample H
        hapgenocand(ali) = 2
        q_ij = q_ij + log(1-prob)
        where (h1 == ali) g1j = 2
        where (h2 == ali) g2j = 2
      endif

      ! Accumulate q_ji as well
      old_geno = hapgeno(ali)
      if (old_geno == 1) then
        q_ji = q_ji + log(prob)
      else
        q_ji = q_ji + log(1-prob)
      endif

    enddo


    ! Get likelihood
    g = g1j + g2j - 2

    where (g == 0)
      likpen_vec = likpen(:,1)
    else where (g == 1)
      likpen_vec = likpen(:,2)
    else where (g == 2)
      likpen_vec = likpen(:,3)
    end where

    pi_j = sum(likpen_vec)

    ! Metropolis Hastings step
    lrat = (pi_j + q_ji) - (pi_i + q_ij)

    accept = .FALSE. 
    if (lrat > 0) then
      accept = .TRUE.
    else 
      call random_number(x)
      if (exp(lrat) > x) then
        accept = .TRUE.
      endif
    endif

    depth_log(srep) = depth
    accept_log(srep) = accept
    pi_i_log(srep) = pi_i
    pi_j_log(srep) = pi_j
    q_ij_log(srep) = q_ij
    q_ji_log(srep) = q_ji

    if (accept) then
      hapgeno = hapgenocand
      pi_i = pi_j
      g1i = g1j
      g2i = g2j
    end if

    geno_hap_log = geno_hap_log + hapgeno
    geno_ani_log = geno_ani_log + g1i + g2i

  enddo

  geno_ani_log = geno_ani_log / nsamp
  geno_ani_log = geno_ani_log  - 2 ! copies of horned
  geno_ani_log = (2 - geno_ani_log)/2 ! Proportion polled

  geno_hap_log = geno_hap_log / nsamp
  geno_hap_log = geno_hap_log - 1 ! P horned
  geno_hap_log = 1 - geno_hap_log ! P polled

  ! Write out results

  !outfile = "Temp/samplog_" // char(rep+ ichar("0")) // ".txt"
  outfile = "Temp/samplog.txt"
  open(unit=11,file=outfile,action="write",status="replace")
  write(11,*) "sample depth pi_i pi_j q_ij q_ji accept"
  do srep = 1,nsamp 
    write(11,"(I10,I6,4F12.3,L3)"),srep,depth_log(srep), pi_i_log(srep), &
        pi_j_log(srep), q_ij_log(srep), q_ji_log(srep),accept_log(srep)
  enddo
  close(11)

  !outfile = "Temp/anilog_" // char(rep+ ichar("0")) // ".txt"
  outfile = "Temp/anilog.txt"
  open(unit=11,file=outfile,action="write",status="replace")
  write(11,*) "id pheno hap1 hap2 phap1 phap2 ppoll last1 last2"
  do srep = 1,nanis 
    write(11,"(A10,A1,A10,2I6,3F9.6,2I6)") id(srep)," ",pheno(srep), &
        h1(srep),h2(srep),geno_hap_log(h1(srep)), &
        geno_hap_log(h2(srep)),geno_ani_log(srep), &
        hapgeno(h1(srep)), hapgeno(h2(srep))
  enddo
  close(11)
  

  !outfile = "Temp/haplog_" // char(rep+ ichar("0")) // ".txt"
  outfile = "Temp/haplog.txt"
  open(unit=11,file=outfile,action="write",status="replace")
  write(11,*) "haplotype ppoll last"
  do srep = 1,nhap 
    write(11,"(I10,F9.6,I5)") srep,geno_hap_log(srep),hapgeno(srep)
  enddo
  close(11)
  



end subroutine runf90code

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!




