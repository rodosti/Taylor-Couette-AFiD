      program papero
      use mpih
      use param
      use local_arrays, only: qr,qz,pr,qt
      use hdf5
      use decomp_2d
      use decomp_2d_fft
#ifdef _OPENMP
      use omp_lib
#endif
      implicit none
      integer :: n,ns,nthreads
      integer :: l, errorcode
      real    :: instCFL,dmax
      real    :: ti(2), tin(3)
      real :: ts,te,minwtdt

      call ReadInputFile

      call decomp_2d_init(nrm,nzm,nthm,0,0, &
       (/ .false.,.true.,.true. /))

      ts=MPI_WTIME()
      tin(1) = MPI_WTIME()

      call MpiBarrier

      call HdfStart

      if (nrank.eq.master) ismaster = .true.

      if (ismaster) write(6,*) 'MPI tasks=', nproc

#ifdef _OPENMP
      if (ismaster) then 
!$OMP PARALLEL
!$OMP MASTER
      nthreads = omp_get_num_threads()
!$OMP END MASTER
!$OMP END PARALLEL
          write(6,*) 'OMP threads=', nthreads
      end if
#endif

      if (ismaster) call OpenLogs

!RO   Define Pi=3.1412....

      pi=2.d0*dasin(1.d0)                          

      if(ismaster) then

      write(6,112) alx3/(rext-rint)
  112 format(//,20x,'TAYLOR-COUETTE',//,10x, &
       'with aspect-ratio:  D/H = ',f5.2)
      write(6,142) 
  142 format(//,8x,'Axial periodicity boundary condition')
      write(6,202) ren*(rext-rint),-romeg*(rext-rint)
  202 format(/,5x,'Parameters: ',' Re_s=',e10.3,' R_om= ',e10.3)
      if(idtv.eq.1) then
         write(6,204) limitCFL
  204 format(/,5x,'Variable dt and fixed cfl= ',e11.4,/ )            
      else 
         write(6,205) dtmax,limitCFL
  205 format(/,5x,'Fixed dt= ',e11.4,' and maximum cfl=',e11.4,/ )
      endif

      endif

      call InitTimeMarchScheme

      call InitVariables
      
      call CreateGrid

      call WriteGridInfo

      call HdfStart

      call InitStats

      if(dumpslabr) call InitDumpRCuts

      if(ismaster) then
      write(6,754)nth,nz,nr                  
  754 format(/,5x,'grid resolution: ',' nth= ',i5,' nz= ',i5, &
       ' nr= ',i5/)                       
      write(6,755) 1.d0/dxt,1.d0/dxz,1.d0/dxr,dt,ntst                  
  755 format(/,2x,' dxt=',e10.3,' dxz=',e10.3,' dxr=',e10.3,' dt=' &
       ,e10.3,' ntst=',i7,/)
      endif

      
      call InitPressureSolver

      vmax=0.


      if(readflow) then 

        if(ismaster) write(6,*) 'Reading initial condition from file'

        call ReadFlowField

      else

        if(ismaster) write(6,*) 'Creating initial condition'

        ntime=0                                                         
        time=0.d0
        instCFL=0.d0
        
        call CreateInitialConditions

      endif 

      call update_halo(qt,1,.TRUE.)
      call update_halo(qr,1,.TRUE.)
      call update_halo(qz,1,.TRUE.)
      call update_halo(pr,1,.TRUE.)

!EP   Check divergence. Should be reduced to machine precision after the first
!phcalc. Here it can still be high.

      call CheckDivergence(dmax)
      if(ismaster) write(6,*)' initial divg dmax  ',dmax

      if(ismaster) write(6,711) ntst,tout

711   format(3x,'check in cond :  ntst =',i8,2x,'tout =',f10.1//)


      if(variabletstep) then
       if(ismaster) write(6,*)ntime,time,vmax(1),vmax(2),vmax(3),dt,dmax
      else
       if(ismaster) write(6,*)ntime,time,vmax(1),vmax(2),vmax(3),instCFL,dmax
      endif

      if(ismaster) then
        tin(2) = MPI_WTIME()
        write(6,'(a,f6.2,a)') 'Initialization Time = ', tin(2) -tin(1), ' sec.'
      endif

     
!
!  ********* starts the time dependent calculation ***

      errorcode = 0 !EP set errocode to 0 (OK)
      minwtdt = huge(0.0d0) !EP initialize minimum time step walltime

!  
      do ntime=0,ntst                                           

       ti(1) = MPI_WTIME()
 
       call CalcMaxCFL(instCFL)
            
        if(variabletstep) then
          if(ntime.gt.1) then
            if(instCFL.lt.1.0d-8) then !EP prevent fp-overflow
              dt=dtmax
            else
              dt=limitCFL/instCFL
            endif
            if(dt.gt.dtmax) dt=dtmax
          else
            dt=dtmin
          endif
            if(dt.lt.dtmin) errorcode = 166
        else  
!RO    fixed time-step
          instCFL=instCFL*dt
          if(instCFL.gt.limitCFL) errorcode = 165
        endif

       call TimeMarcher

       time=time+dt

        if(ntime.eq.1.or.mod(time,tout).lt.dt) then
          call CalcVelocityStats
          if(vmax(1).gt.limitVel.and.vmax(2).gt.limitVel) errorcode = 266

          call CalcMaxCFL(instCFL)
          call CheckDivergence(dmax)
          call CalcTorque

           if(time.gt.tsta) then
            call CalcStats
            call CalcDissipation
            if(mod(time,tframe2d).lt.dt) then
              if (dumpslabr) call DumpRCuts
              if (dumpslabt) call DumpThCuts
              if (dumpstave) call DumpStreamwiseAvg
            endif
            if(mod(time,tframe3d).lt.dt) then
              if (dumpfullf) call DumpFlowField
            endif
           endif

            if(.not.variabletstep) instCFL=instCFL*dt

            if(abs(dmax).gt.resid) errorcode = 169
       end if

       if(time.gt.tmax) errorcode=333

        ti(2) = MPI_WTIME()
        minwtdt = min(minwtdt,ti(2) - ti(1))
        if(mod(time,tout).lt.dt) then
          if(ismaster) then
          write(6,*) 'Maximum divergence = ', dmax  
          write(6,*) 'ntime - time - vmax(1) - vmax(2) - vmax(3)' 
          write(6,*) ntime,time,vmax(1),vmax(2),vmax(3)
          write(6,'(a,f8.3,a)') 'Minimum Iteration Time = ', minwtdt, &
                  ' sec.'
          endif
          minwtdt = huge(0.0d0)
        endif

       if( (ti(2) - tin(1)) .gt. walltimemax) errorcode = 334

       if( ntime .eq. ntst ) errorcode = 555 


      call MpiBcastInt(errorcode)

!EP   Conditional exits
      if(errorcode.ne.0) then

!EP    dt too small
        if(errorcode.eq.166) call QuitRoutine(tin,.false.,errorcode)

!EP   cfl too high    
        if(errorcode.eq.165) call QuitRoutine(tin,.false.,errorcode)
      
!EP   velocities diverged
        if(errorcode.eq.266) call QuitRoutine(tin,.false.,errorcode)
          
!EP   mass not conserved
        if(errorcode.eq.169) call QuitRoutine(tin,.false.,errorcode)

!EP   Physical time exceeded tmax, no error; normal quit
        if(errorcode.eq.333) call QuitRoutine(tin,.true.,errorcode)

!EP   walltime exceeded walltimemax, no error; normal quit
        if(errorcode.eq.334) call QuitRoutine(tin,.true.,errorcode)

!RS   FFT input not correct
        if(errorcode.eq.444) call QuitRoutine(tin,.false.,errorcode)

!RS   maximum number of timesteps reached, no error; normal quit
        if(errorcode.eq.555) call QuitRoutine(tin,.true.,errorcode)

        errorcode = 100 !EP already finalized
      
        exit

      endif

      enddo !EP main loop

      call QuitRoutine(tin,.true.,errorcode)
      
      stop                                                              
      end                                                               
