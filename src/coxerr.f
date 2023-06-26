      subroutine coxerr(t,dlt,wuz,size,npred,mcode,derr,
     &     bt,va,succind,s,ew,db,nb,ef,deri,s1,s2)
      integer size,npred,dlt(size),mcode,succind,s(size)
      double precision t(size),wuz(0:npred,size),derr,
     &     bt(npred),va(npred,npred),ew(size),db(npred),
     &     nb(npred),ef(npred),deri(npred,npred),
     &     s1(0:npred),s2(npred,npred)

      external prop1,prop2
      logical succ,init
      integer i,j,p,q,tmp,ifr
      double precision ob

      do i=1,size
         s(i)=i
      enddo
      do i=2,size
         j=i
         do while(j .gt. 1 .and.
     &        (t(s(j-1)) .gt. t(s(j)) .or.
     &        (t(s(j-1)) .eq. t(s(j)) .and.
     &        (dlt(s(j-1)) .lt. dlt(s(j))))))
            tmp=s(j)
            s(j)=s(j-1)
            s(j-1)=tmp
            j=j-1
         enddo
      enddo

      if(mcode .eq. 1)
     &     call solver(prop1,derr,
     &     dlt,wuz,s,size,npred,bt,succ,ew,db,nb,ef,deri,s1,s2)
      if(mcode .eq. 2)
     &     call solver(prop2,derr,
     &     dlt,wuz,s,size,npred,bt,succ,ew,db,nb,ef,deri,s1,s2)

      succind=0
      if(succ) succind=1

      call inverse(deri,npred,ifr,va,s2)
      init=.false.
      if(mcode .eq. 1)
     &     call prop1(1,dlt,wuz,s,size,npred,bt,ob,ef,va,
     &     init,ew,s1,s2)
      if(mcode .eq. 2)
     &     call prop2(1,dlt,wuz,s,size,npred,bt,ob,ef,va,
     &     init,ew,s1,s2)

      do p=1,npred
         do q=1,npred
            s2(p,q)=0.0d0
            do tmp=1,npred
               s2(p,q)=s2(p,q)+deri(p,tmp)*va(tmp,q)
            enddo
         enddo
      enddo
      do p=1,npred
         do q=1,p
            va(p,q)=0.0d0
            do tmp=1,npred
               va(p,q)=va(p,q)+s2(p,tmp)*deri(q,tmp)
            enddo
         enddo
      enddo
      do p=1,npred
         do q=p+1,npred
            va(p,q)=va(q,p)
         enddo
      enddo

      end

      subroutine solver(method,derr,
     &     dlt,wuz,s,size,npred,bt,succ,ew,db,nb,ef,deri,s1,s2)
      external method
      logical succ
      integer size,npred,dlt(size),s(size)
      double precision derr,wuz(0:npred,size),
     &     bt(npred),ew(size),db(npred),nb(npred),ef(npred),
     &     deri(npred,npred),s1(0:npred),s2(npred,npred)

      logical init,cont
      integer p,ifr,cnt
      double precision ob,dtmp

      succ=.true.
      init=.true.
      call method(0,dlt,wuz,s,size,npred,bt,ob,ef,deri,
     &     init,ew,s1,s2)

      cont=.true.
      if(ob .le. derr) cont=.false.
      do while(cont)
         dtmp=ob
         do p=1,npred
            db(p)=ef(p)
         enddo
         call axb(deri,npred,ifr,db)

         cnt=0
         do while(cnt .eq. 0 .or.
     &        (ob .ge. dtmp .and. cnt .le. 20))
            do p=1,npred
               nb(p)=bt(p)-db(p)/2.0d0**dble(cnt)
            enddo
            call method(0,dlt,wuz,s,size,npred,nb,ob,ef,deri,
     &           init,ew,s1,s2)
            cnt=cnt+1
            if(isnan(ob)) cnt=21
         enddo

         if(cnt .eq. 21) then
            cont=.false.
            succ=.false.
            call method(0,dlt,wuz,s,size,npred,bt,ob,ef,deri,
     &           init,ew,s1,s2)
         else
            do p=1,npred
               bt(p)=nb(p)
            enddo
            if(ob .le. derr) cont=.false.
         endif
      enddo

      end

      subroutine prop1(ord,dlt,wuz,s,size,npred,bt,ob,ef,deri,
     &     init,ew,s1,s2)
      logical init
      integer ord,size,npred,dlt(size),s(size)
      double precision wuz(0:npred,size),bt(npred),
     &     ob,ef(npred),deri(npred,npred),ew(size),
     &     s1(0:npred),s2(npred,npred)

      integer i,p,q
      double precision s0,dtmp

      s0=0.0d0
      do p=1,npred
         s1(p)=0.0d0
      enddo
      if(ord .eq. 0) then
         do p=1,npred
            ef(p)=0.0d0
         enddo
         s1(0)=0.0d0
         do p=1,npred
            do q=p,npred
               deri(p,q)=0.0d0
               s2(p,q)=0.0d0
            enddo
         enddo
         do p=2,npred
            deri(p,1)=0.0d0
            s2(p,1)=0.0d0
         enddo
      endif
      do i=size,1,-1
         dtmp=bt(1)*wuz(0,s(i))/2.0d0
         do p=2,npred
            dtmp=dtmp+bt(p)*wuz(p,s(i))
         enddo
         dtmp=dexp(dtmp)

         s0=s0+dtmp
         do p=1,npred
            s1(p)=s1(p)+dtmp*wuz(p,s(i))
         enddo
         if(ord .eq. 0) then
            s1(0)=s1(0)+dtmp*wuz(0,s(i))/2.0d0
            do p=1,npred
               do q=p,npred
                  if(q .ne. 1) s2(p,q)=s2(p,q)
     &                 +dtmp*wuz(p,s(i))*wuz(q,s(i))
               enddo
            enddo
            do p=1,npred
               s2(p,1)=s2(p,1)+dtmp*wuz(p,s(i))
     &              *wuz(0,s(i))/2.0d0
            enddo
            if(dlt(s(i)) .eq. 1) then
               dtmp=dexp(-bt(1)*wuz(0,s(i))/2.0d0)
               do p=1,npred
                  ef(p)=ef(p)+(wuz(p,s(i))-s1(p)/s0)*dtmp
               enddo
               do p=1,npred
                  do q=p,npred
                     if(q .ne. 1) deri(p,q)=deri(p,q)
     &                    -(s2(p,q)/s0-s1(p)*s1(q)/s0**2.0d0)
     &                    *dtmp
                  enddo
               enddo
               do p=1,npred
                  deri(p,1)=deri(p,1)-(s2(p,1)/s0
     &                 -s1(p)*s1(0)/s0**2.0d0)*dtmp
     &                 -(wuz(p,s(i))-s1(p)/s0)*dtmp
     &                 *wuz(0,s(i))/2.0d0
               enddo
            endif
         endif
      enddo

      if(ord .eq. 0) then
         do p=1,npred
            ef(p)=ef(p)/dble(size)
         enddo
         do p=1,npred
            do q=p,npred
               deri(p,q)=deri(p,q)/dble(size)
            enddo
         enddo
         do p=2,npred
            deri(p,1)=deri(p,1)/dble(size)
         enddo 
         do p=3,npred
            do q=2,p-1
               deri(p,q)=deri(q,p)
            enddo
         enddo

         ob=0.0d0
         do p=1,npred
            ob=ob+ef(p)**2.0d0
         enddo
         ob=dsqrt(ob)
      else
         c0=0.0d0
         do p=1,npred
            s2(p,1)=0.0d0
            do q=p,npred
               deri(p,q)=0.0d0
            enddo
         enddo
         do i=1,size
            if(dlt(s(i)) .eq. 1) then
               dtmp=dexp(-bt(1)*wuz(0,s(i))/2.0d0)
               c0=c0+dtmp/s0
               do p=1,npred
                  s2(p,1)=s2(p,1)+s1(p)/s0**2.0d0*dtmp
               enddo
               do p=1,npred
                  ef(p)=(wuz(p,s(i))-s1(p)/s0)*dtmp
               enddo
            else
               do p=1,npred
                  ef(p)=0.0d0
               enddo
            endif

            dtmp=bt(1)*wuz(0,s(i))/2.0d0
            do p=2,npred
               dtmp=dtmp+bt(p)*wuz(p,s(i))
            enddo
            dtmp=dexp(dtmp)
            do p=1,npred
               ef(p)=ef(p)-(c0*wuz(p,s(i))-s2(p,1))*dtmp
            enddo

            do p=1,npred
               do q=p,npred
                  deri(p,q)=deri(p,q)+ef(p)*ef(q)
               enddo
            enddo

            s0=s0-dtmp
            do p=1,npred
               s1(p)=s1(p)-dtmp*wuz(p,s(i))
            enddo
         enddo

         do p=1,npred
            do q=p,npred
               deri(p,q)=deri(p,q)/dble(size)**2.0d0
            enddo
         enddo
         do p=1,npred
            do q=1,p-1
               deri(p,q)=deri(q,p)
            enddo
         enddo
      endif

      end

      subroutine prop2(ord,dlt,wuz,s,size,npred,bt,ob,ef,deri,
     &     init,ew,s1,s2)
      logical init
      integer ord,size,npred,dlt(size),s(size)
      double precision wuz(0:npred,size),bt(npred),
     &     ob,ef(npred),deri(npred,npred),ew(size),
     &     s1(0:npred),s2(npred,npred)

      integer i,p,q,ifr
      double precision s0,dtmp

      if(init) then
         s1(0)=0.0d0
         do p=1,npred
            ef(p)=0.0d0
            s1(p)=0.0d0
            do q=p,npred
               deri(p,q)=0.0d0
            enddo
         enddo
         do i=1,size
            s1(0)=s1(0)+wuz(0,i)
            do p=1,npred
               ef(p)=ef(p)+wuz(0,i)*wuz(p,i)
               s1(p)=s1(p)+wuz(p,i)
               do q=p,npred
                  deri(p,q)=deri(p,q)+wuz(p,i)*wuz(q,i)
               enddo
            enddo
         enddo
         do p=0,npred
            s1(p)=s1(p)/dble(size)
         enddo
         do p=1,npred
            ef(p)=ef(p)/dble(size)-s1(0)*s1(p)
            do q=p,npred
               deri(p,q)=deri(p,q)/dble(size)-s1(p)*s1(q)
            enddo
         enddo
         do p=1,npred
            do q=1,p-1
               deri(p,q)=deri(q,p)
            enddo
         enddo
         call axb(deri,npred,ifr,ef)
         do i=1,size
            ew(i)=0.0d0
            do p=1,npred
               ew(i)=ew(i)+ef(p)*wuz(p,i)
            enddo
         enddo

         init=.false.
      endif

      s0=0.0d0
      do p=1,npred
         s1(p)=0.0d0
      enddo
      if(ord .eq. 0) then
         do p=1,npred
            ef(p)=0.0d0
         enddo
         s1(0)=0.0d0
         do p=1,npred
            do q=p,npred
               deri(p,q)=0.0d0
               s2(p,q)=0.0d0
            enddo
         enddo
         do p=2,npred
            deri(p,1)=0.0d0
            s2(p,1)=0.0d0
         enddo
      endif
      do i=size,1,-1
         dtmp=bt(1)*(wuz(0,s(i))+ew(s(i)))/2.0d0
         do p=2,npred
            dtmp=dtmp+bt(p)*wuz(p,s(i))
         enddo
         dtmp=dexp(dtmp)

         s0=s0+dtmp
         do p=1,npred
            s1(p)=s1(p)+dtmp*wuz(p,s(i))
         enddo
         if(ord .eq. 0) then
            s1(0)=s1(0)+dtmp*(wuz(0,s(i))+ew(s(i)))/2.0d0
            do p=1,npred
               do q=p,npred
                  if(q .ne. 1) s2(p,q)=s2(p,q)
     &                 +dtmp*wuz(p,s(i))*wuz(q,s(i))
               enddo
            enddo
            do p=1,npred
               s2(p,1)=s2(p,1)+dtmp*wuz(p,s(i))
     &              *(wuz(0,s(i))+ew(s(i)))/2.0d0
            enddo
            if(dlt(s(i)) .eq. 1) then
               dtmp=dexp(-bt(1)*(wuz(0,s(i))-ew(s(i)))/2.0d0)
               do p=1,npred
                  ef(p)=ef(p)+(wuz(p,s(i))-s1(p)/s0)*dtmp
               enddo
               do p=1,npred
                  do q=p,npred
                     if(q .ne. 1) deri(p,q)=deri(p,q)
     &                    -(s2(p,q)/s0-s1(p)*s1(q)/s0**2.0d0)
     &                    *dtmp
                  enddo
               enddo
               do p=1,npred
                  deri(p,1)=deri(p,1)-(s2(p,1)/s0
     &                 -s1(p)*s1(0)/s0**2.0d0)*dtmp
     &                 -(wuz(p,s(i))-s1(p)/s0)*dtmp
     &                 *(wuz(0,s(i))-ew(s(i)))/2.0d0
               enddo
            endif
         endif
      enddo

      if(ord .eq. 0) then
         do p=1,npred
            ef(p)=ef(p)/dble(size)
         enddo
         do p=1,npred
            do q=p,npred
               deri(p,q)=deri(p,q)/dble(size)
            enddo
         enddo
         do p=2,npred
            deri(p,1)=deri(p,1)/dble(size)
         enddo 
         do p=3,npred
            do q=2,p-1
               deri(p,q)=deri(q,p)
            enddo
         enddo

         ob=0.0d0
         do p=1,npred
            ob=ob+ef(p)**2.0d0
         enddo
         ob=dsqrt(ob)
      else
         c0=0.0d0
         do p=1,npred
            s2(p,1)=0.0d0
            do q=p,npred
               deri(p,q)=0.0d0
            enddo
         enddo
         do i=1,size
            if(dlt(s(i)) .eq. 1) then
               dtmp=dexp(-bt(1)*(wuz(0,s(i))-ew(s(i)))/2.0d0)
               c0=c0+dtmp/s0
               do p=1,npred
                  s2(p,1)=s2(p,1)+s1(p)/s0**2.0d0*dtmp
               enddo
               do p=1,npred
                  ef(p)=(wuz(p,s(i))-s1(p)/s0)*dtmp
               enddo
            else
               do p=1,npred
                  ef(p)=0.0d0
               enddo
            endif

            dtmp=bt(1)*(wuz(0,s(i))+ew(s(i)))/2.0d0
            do p=2,npred
               dtmp=dtmp+bt(p)*wuz(p,s(i))
            enddo
            dtmp=dexp(dtmp)
            do p=1,npred
               ef(p)=ef(p)-(c0*wuz(p,s(i))-s2(p,1))*dtmp
            enddo

            do p=1,npred
               do q=p,npred
                  deri(p,q)=deri(p,q)+ef(p)*ef(q)
               enddo
            enddo

            s0=s0-dtmp
            do p=1,npred
               s1(p)=s1(p)-dtmp*wuz(p,s(i))
            enddo
         enddo

         do p=1,npred
            do q=p,npred
               deri(p,q)=deri(p,q)/dble(size)**2.0d0
            enddo
         enddo
         do p=1,npred
            do q=1,p-1
               deri(p,q)=deri(q,p)
            enddo
         enddo
      endif

      end

      subroutine inverse(a,dim,ifr,b,c)
      integer dim,ifr
      double precision a(dim,dim),b(dim,dim),
     &     c(dim,dim)

      integer p,q,n

      ifr=1
      do p=1,dim
         do q=1,p-1
            b(q,p)=0.0d0
            do n=1,dim
               b(q,p)=b(q,p)+a(n,p)*a(n,q)
            enddo
            do n=1,dim
               a(n,p)=a(n,p)-b(q,p)*a(n,q)
            enddo
         enddo
         b(p,p)=0.0d0
         do n=1,dim
            b(p,p)=b(p,p)+a(n,p)**2
         enddo
         b(p,p)=dsqrt(b(p,p))
         if(b(p,p) .gt. 1.0d-10) then
            do n=1,dim
               a(n,p)=a(n,p)/b(p,p)
            enddo
         else
            ifr=0
            b(p,p)=1.0d0
         endif
      enddo

      do p=dim,1,-1
         do q=1,p-1
            c(p,q)=0.0d0
         enddo
         c(p,p)=1.0d0/b(p,p)
         do q=p+1,dim
            c(p,q)=0.0d0
            do n=p+1,q
               c(p,q)=c(p,q)-b(p,n)*c(n,q)
            enddo
            c(p,q)=c(p,q)/b(p,p)
         enddo
      enddo

      do p=1,dim
         do q=1,dim
            b(p,q)=0.0d0
            do n=1,dim
               b(p,q)=b(p,q)+c(p,n)*a(q,n)
            enddo
         enddo
      enddo

      do p=1,dim
         do q=1,dim
            a(p,q)=b(p,q)
         enddo
      enddo

      end

      subroutine axb(a,dim,ifr,b)
      integer dim,ifr
      double precision a(dim,dim),b(dim)

      integer p,q,i
      double precision small,dtmp,dcmp

      small=1.0d-10
      ifr=1

      do p=1,dim-1
         i=p
         dtmp=dabs(a(p,p))
         if(dtmp .le. small) then
            do q=p+1,dim
               if(dabs(a(q,p)) .gt. dtmp) then
                  dtmp=dabs(a(q,p))
                  i=q
               endif
            enddo
            if(dtmp .gt. small) then
               do q=p,dim
                  dcmp=a(p,q)
                  a(p,q)=a(i,q)
                  a(i,q)=dcmp
               enddo
               dcmp=b(p)
               b(p)=b(i)
               b(i)=dcmp
            endif
         endif

         if(dtmp .gt. small) then
            do q=p+1,dim
               dcmp=-a(q,p)/a(p,p)
               do i=p+1,dim
                  a(q,i)=a(q,i)+a(p,i)*dcmp
               enddo
               b(q)=b(q)+b(p)*dcmp
            enddo
         endif
      enddo

      do p=dim,1,-1
         do q=p+1,dim
            b(p)=b(p)-b(q)*a(p,q)
         enddo
         if(dabs(a(p,p)) .le. small) then
            ifr=0
            b(p)=0.0d0
         else
            b(p)=b(p)/a(p,p)
         endif
      enddo

      end
