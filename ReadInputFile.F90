!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                         ! 
!    FILE: ReadInputFile.F90                              !
!    CONTAINS: subroutine ReadInputFile                   !
!                                                         ! 
!    PURPOSE: Read parameters from bou.in file            !
!                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine ReadInputFile
      use param
      use mpih
      implicit none
      logical :: fexist
      character(len=4) :: dummy

       open(unit=15,file='bou.in',status='old')
        read(15,301) dummy
        read(15,*) nthm,nzm,nrm,nsst,idtv,nread
        read(15,301) dummy
        read(15,*) ntst,tout,tmax,ireset,walltimemax
        read(15,301) dummy
        read(15,*) alx3,rint,rext,istr,str,lamb
        read(15,301) dummy
        read(15,*) virotint,virotext,sizepert,omegapert
        read(15,301) dummy
        read(15,*) ren,romeg,resid,limitCFL
        read(15,301) dummy
        read(15,*) tsta,starea
        read(15,301) dummy       
        read(15,*) dtmin,dtmax,limitVel,cfllim  
        read(15,301) dummy       
        read(15,*) nslabr, dpthcut, dpstave, tframe2d
        read(15,301) dummy       
        read(15,*) dpfullf, tframe3d
        read(15,301) dummy       
        read(15,*) minstrt, minstrz, multmin
301     format(a4)                
      close(15)

      nth=nthm+1     
      nz=nzm+1                                                          
      nr=nrm+1

      if(idtv.eq.0) variabletstep = .false.
      if(nread.ne.0) readflow = .true.
      if(ireset.ne.0) resetlogstime = .true.

      if(nslabr.ne.0) then
       inquire(file='./rcuts.in', exist=fexist) 
       if(fexist) then
        dumpslabr = .true.
       else
        write(6,*) "rcuts.in not found turning off slab dump"
       end if
      endif

      if(dpthcut.ne.0) dumpslabt = .true.
      if(dpstave.ne.0) dumpstave = .true.
      if(dpfullf.ne.0) dumpfullf = .true.

      if(starea.ne.0) then 
       inquire(file='stats/stafield_master.h5', exist=fexist) 
       if(fexist) then
        readstats = .true.
        if (.not. readflow) write(6,*) 'Warning: Restarting', &
         ' flowfield with statistics read'
       else
        write(6,*) "stats data not found"
       end if
      endif

      return
      end
