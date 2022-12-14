!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                         ! 
!    FILE: interp.F90                                     !
!    CONTAINS: subroutine interp,gridnew,interptrilin     !
!                                                         ! 
!    PURPOSE: Trilinear interpolation to a new grid for   !
!     continuation files.                                 !
!                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine interp(arrold,arrnew,ntho,nzo,nro, &
       istro3,stro3,intvar,xs2o,xe2o,xs3o,xe3o)
      use param
      use mpih
      use decomp_2d, only: xstart,xend
      implicit none
      integer :: istro3
      real :: stro3, thmax
      integer,intent(in) :: intvar,nzo,nro,ntho
      integer,intent(in) :: xs2o,xe2o,xs3o,xe3o
 
      real,intent(out),dimension(1:nr,xstart(2):xend(2), &
       xstart(3):xend(3)) :: arrnew
      real,dimension(0:nro+1,xs2o-1:xe2o+1,xs3o-1:xe3o+1) :: arrold
      real,dimension(0:ntho+1) :: tcold,tmold
      real,dimension(1:nth) :: tcnew,tmnew
      real,dimension(0:nzo+1) :: zcold,zmold
      real,dimension(1:nz) :: zcnew,zmnew
      real,dimension(0:nro+1) :: rcold,rmold
      real,dimension(1:nr) :: rcnew,rmnew

      real, dimension(0:ntho+1) :: xold
      real, dimension(0:nzo+1) :: yold
      real, dimension(0:nro+1) :: zold
      real, dimension(1:nth) :: xnew
      real, dimension(1:nz) :: ynew
      real, dimension(1:nr) :: znew

      real :: bn(6),an(8)

      real :: bix1,bix2,biy1,biy2,biz1,biz2,bix,biy,biz
      integer :: j,k,i,l
      real :: sfcf,sfff,sccf,scff,sfcc,sffc,sccc,scfc
      integer :: bici,bifi,bicj,bifj,bick,bifk

!EP   Create old grid
      pi=2.d0*asin(1.d0)

      thmax=2.0*pi/lamb

      call gridnew(ntho,0.0d0,thmax,0,1.0,tcold(1:ntho),tmold(1:ntho))
      call gridnew(nzo,0.0d0,alx3,0,1.0,zcold(1:nzo),zmold(1:nzo))
      call gridnew(nro,rint,rext,istro3,stro3,rcold(1:nro),rmold(1:nro))

!EP   2nd order extrapolation of grid
      tcold(0) = 2*tcold(1)-tcold(2)
      tcold(ntho+1) = 2*tcold(ntho)-tcold(ntho-1)
      tmold(0) = 2*tmold(1)-tmold(2)
      tmold(ntho+1) = 2*tmold(ntho)-tmold(ntho-1)

      zcold(0) = 2*zcold(1)-zcold(2)
      zcold(nzo+1) = 2*zcold(nzo)-zcold(nzo-1)
      zmold(0) = 2*zmold(1)-zmold(2)
      zmold(nzo+1) = 2*zmold(nzo)-zmold(nzo-1)

      rcold(0) = 2*rcold(1)-rcold(2)
      rcold(nro+1) = 2*rcold(nro)-rcold(nro-1)
      rmold(0) = 2*rmold(1)-rmold(2)
      rmold(nro+1) = 2*rmold(nro)-rmold(nro-1)

!EP   Create new grid
      call gridnew(nth,0.0d0,thmax,0,1.0,tcnew,tmnew)
      call gridnew(nz,0.0d0,alx3,0,1.0,zcnew,zmnew)
      call gridnew(nr,rint,rext,istr,str,rcnew,rmnew)
      
      select case (intvar)
          case (1)
             xold=tcold
             yold=zmold
             zold=rmold
             xnew=tcnew
             ynew=zmnew
             znew=rmnew
          case (2)
             xold=tmold
             yold=zcold
             zold=rmold
             xnew=tmnew
             ynew=zcnew
             znew=rmnew
          case (3)
             xold=tmold
             yold=zmold
             zold=rcold
             xnew=tmnew
             ynew=zmnew
             znew=rcnew
          case (4)
             xold=tmold
             yold=zmold
             zold=rmold
             xnew=tmnew
             ynew=zmnew
             znew=rmnew
      end select

!EP   Extrapolate halo

      do i=xs3o,xe3o
        do j=xs2o,xe2o
          arrold(0    ,j,i)=2.0*arrold(1  ,j,i)-arrold(2,j,i)
          arrold(nro+1,j,i)=2.0*arrold(nro,j,i)-arrold(nro-1,j,i)
        enddo
      enddo

      do j=xs2o,xe2o
        do k=1,nro
          arrold(k,j,xs3o-1)=2.0*arrold(k,j,xs3o)-arrold(k,j,xs3o+1)
          arrold(k,j,xe3o+1)=2.0*arrold(k,j,xe3o)-arrold(k,j,xe3o-1)
        enddo
      enddo

      do i=xs3o,xe3o
        do k=1,nro
          arrold(k,xs2o-1,i)=2.0*arrold(k,xs2o,i)-arrold(k,xs2o+1,i)
          arrold(k,xe2o+1,i)=2.0*arrold(k,xe2o,i)-arrold(k,xe2o-1,i)
        enddo
      enddo


!EP   INTERP
      do i=xstart(3),xend(3)
       do j=xstart(2),xend(2)
        do k=1,nr
!    Find nearest grid value
      bix=xnew(i)
      biy=ynew(j)
      biz=znew(k)

      bifi=ntho-1
      bici=ntho
      do l=1,ntho
      if(xold(l).ge.bix) then
      bifi=l-1
      bici=l
      goto 10
      endif
      enddo
10    continue

      bifj=nzo-1
      bicj=nzo
      do l=1,nzo
      if(yold(l).ge.biy) then
      bifj=l-1
      bicj=l
      goto 20
      endif
      enddo
20    continue

      bifk=nro-1
      bick=nro
      do l=1,nro
      if(zold(l).ge.biz) then
      bifk=l-1
      bick=l
      goto 30
      endif
      enddo
30    continue

!EP   Define
      bix1 = xold(bifi) 
      bix2 = xold(bici)
      biy1 = yold(bifj) 
      biy2 = yold(bicj)
      biz1 = zold(bifk) 
      biz2 = zold(bick)

!EP   Send data

       sfcf = arrold(bifk,bicj,bifi)
       sfff = arrold(bifk,bifj,bifi)
       sccf = arrold(bifk,bicj,bici)
       scff = arrold(bifk,bifj,bici)

       sfcc = arrold(bifk,bicj,bifi)
       sffc = arrold(bifk,bifj,bifi)
       sccc = arrold(bifk,bicj,bici)
       scfc = arrold(bifk,bifj,bici)

       an(1) = sfcf
       an(2) = sfff
       an(3) = sccf
       an(4) = scff
       an(5) = sfcc
       an(6) = sffc
       an(7) = sccc
       an(8) = scfc
       bn(1) = bix1
       bn(2) = bix2
       bn(3) = biy1
       bn(4) = biy2
       bn(5) = biz1
       bn(6) = biz2
       call interptrilin(bix,biy,biz,an,bn,arrnew(k,j,i))

        enddo
       enddo
      enddo

      end

      subroutine interptrilin(bix,biy,biz,an,bn,ans)
      implicit none
      
      real, intent(in) :: bix,biy,biz
      real,intent(in) :: bn(6),an(8)
      real, intent(out) :: ans

      real :: bix1,bix2,biy1,biy2,biz1,biz2
      real :: afifjck,acifjck,aficjck,acicjck
      real :: afifjfk,acifjfk,aficjfk,acicjfk,dxdydz
      

      aficjfk = an(1)
      afifjfk = an(2)
      acicjfk = an(3)
      acifjfk = an(4)
      aficjck = an(5)
      afifjck = an(6)
      acicjck = an(7)
      acifjck = an(8)
      
      bix1 = bn(1)
      bix2 = bn(2)
      biy1 = bn(3)
      biy2 = bn(4)
      biz1 = bn(5)
      biz2 = bn(6)

      dxdydz = 1.0/((bix2-bix1)*(biy2-biy1)*(biz2-biz1))

      ans = (aficjfk*(bix2-bix)*(biy-biy1)*(biz2-biz) &
            +afifjfk*(bix2-bix)*(biy2-biy)*(biz2-biz) &
            +acicjfk*(bix-bix1)*(biy-biy1)*(biz2-biz) &
            +acifjfk*(bix-bix1)*(biy2-biy)*(biz2-biz) &
            +aficjck*(bix2-bix)*(biy-biy1)*(biz-biz1) &
            +afifjck*(bix2-bix)*(biy2-biy)*(biz-biz1) &
            +acicjck*(bix-bix1)*(biy-biy1)*(biz-biz1) &
            +acifjck*(bix-bix1)*(biy2-biy)*(biz-biz1) &
            )*dxdydz

      end

      subroutine gridnew(n,rint,rext,istr,str,rc,rm)
      implicit none
      double precision,intent(in) :: rint,rext,str
      double precision,dimension(1:n),intent(out) :: rc,rm
      double precision,dimension(1:n) :: etaz
      double precision,dimension(1:n+400) :: etazm
      double precision :: x3,etain,delet,pi
      integer,intent(in) :: n,istr
      integer :: i,nrmo,nclip

      if (istr.eq.0) then
        do i=1,n
          x3=real(i-1)/real(n-1)
          rc(i)=rext*x3+rint
        enddo
       endif

      if(istr.eq.6) then
      pi=2.0d0*asin(1.0d0)
      nclip = int(str)
      nrmo = n+nclip+nclip
      do i=1,nrmo
        etazm(i)=+cos(pi*(float(i)-0.5)/float(nrmo))
      end do
      do i=1,n
        etaz(i)=etazm(i+nclip)
      end do
      delet = etaz(1)-etaz(n)
      etain = etaz(1)
      do i=1,n
        etaz(i)=etaz(i)/(0.5*delet)
      end do
      rc(1) = rint
      do i=2,n-1
        rc(i) = (rext-rint)*(1.-etaz(i))*0.5+rint
      end do
      rc(n) = rext
      endif


      do i=1,n-1
        rm(i)=(rc(i)+rc(i+1))*0.5d0
      enddo
      rm(n) = 2*rc(n)-rm(n-1)

      return
      end
