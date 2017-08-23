! Test the PopIII star module
! Sam Geen, August 2017

program testmodule

  use lookup_table_module
  use popiiistar_module
  implicit none
  !integer,parameter::dp=kind(1.0D0) ! real*8  
  !real(dp),parameter::Myrins = 3.1556926d13
  real(dp)::accin,massin
  real(dp)::radius,temperature
  character(len=128)::dir

  ! Test the basic table setup
  !debug_lkup=.true.
  dir = "../"
  call setup_popiii_tables(dir)
  accin = 2d21
  massin = 8d33
  write(*,*) "INPUT ACCRETION, MASS IN SOLAR UNITS", accin/solar_mass_per_year_cgs, massin/solar_mass_cgs
  call popiii_lookup(accin,massin,radius,temperature)
  write(*,*) "RADIUS, TEMPERATURE", radius, temperature
  write(*,*) "RADIUS IN LOG SOLAR UNITS", log10(radius/solar_radius_cgs)
  write(*,*) "TEMPERATURE IN LOG K", log10(temperature)

  ! Output a sample of tracks
!   do i=5,120,5
!      call output_track(real(i,dp))
!   end do
   call output_track(1d-5)
   call output_track(3d-5)
   call output_track(1d-4)
   call output_track(3d-4)
   call output_track(1d-3)
!   call output_track(62.5d0)
!   call output_track(1d0)
!   call output_track(200d0)

end program testmodule


SUBROUTINE Replace_Text (s,text,rep,out)
implicit none
CHARACTER(*):: text,rep
CHARACTER(1)::c
CHARACTER(20):: s, out
INTEGER:: i, nt, nr

out = s
DO i=1,LEN(s)
   c = s(i:i)
   if (c.eq.text) then
      out(i:i) = rep
   endif
END DO
END SUBROUTINE Replace_Text

SUBROUTINE output_track(accretion)

  use popiiistar_module
  implicit none
  
  real(dp),intent(in)::accretion
  real(dp)::mass,radius,temperature,mmax,acc
  integer::i
  integer,parameter::ilun=101
  integer,parameter::nt=1000
  character(len=20)::astr,tmp
  character(len=200)::filename
  real(dp),dimension(nt)::ms,rs,ts
  ! Set up file
  write(astr,'(F10.6)') accretion
  call replace_text(astr," ","0",tmp)
  call replace_text(tmp,".","p",astr)
  filename="../track"//TRIM(astr)//".dat"
  open(ilun,file=filename,form="unformatted")
  ! Read data
  mmax = 80d0*solar_mass_cgs
  write(*,*) accretion
  acc = accretion*solar_mass_per_year_cgs
  do i=1,nt
    mass = (i*mmax)/real(nt,dp)
    call popiii_lookup(acc,mass,radius,temperature)
    !write(*,*) "RADIUS, TEMPERATURE", radius, temperature
    !write(*,*) "RADIUS IN LOG SOLAR UNITS", log10(radius/solar_radius_cgs)
    !write(*,*) "TEMPERATURE IN LOG K", log10(temperature)
    ms(i) = mass/solar_mass_cgs
    rs(i) = log10(radius/solar_radius_cgs)
    ts(i) = log10(temperature)
  enddo
  ! Make cumulative values
  !do i=2,nt
  !   es(i) = es(i)+es(i-1)
  !   mls(i) = mls(i)+mls(i-1)
  !   np1s(i) = np1s(i)+np1s(i-1)
  !   np2s(i) = np2s(i)+np2s(i-1)
  !   np3s(i) = np3s(i)+np3s(i-1)
  !enddo
  ! Write
  write(ilun) accretion
  write(ilun) ms
  write(ilun) rs
  write(ilun) ts
  close(ilun)
  
END SUBROUTINE output_track
