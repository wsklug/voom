c**********************************************************************c

      subroutine deltri(numpts,n,x,list,bin,v,e,numtri)
      implicit double precision (a-h,o-z)
c
      integer n,i,list(*),v(3,*),e(3,*),numtri,bin(*),p,numpts
c
      real*8 xmin,xmax,ymin,ymax,dmax,c00001,fact,x(2,*)
c
      parameter (c00001=1.0)
c
c
c     Compute min and max coords for x and y
c     Compute max overall dimension
c
      xmin = x(1,list(1))
      xmax = xmin
      ymin = x(2,list(1))
      ymax = ymin
c
      do 5 i = 2, n
         p    = list(i)
         xmin = min(xmin,x(1,p))
         xmax = max(xmax,x(1,p))
         ymin = min(ymin,x(2,p))
         ymax = max(ymax,x(2,p))
    5 continue
      dmax = max(xmax-xmin,ymax-ymin)
c
c     Normalize (x,y) coordinates of points
c
      fact = c00001/dmax
      do 10 i = 1, n
         p    = list(i)
         x(1,p) = (x(1,p)-xmin) * fact
         x(2,p) = (x(2,p)-ymin) * fact
   10 continue
c
c     Sort points into bins
c     This call is optional
c
      call bsort (n,x,xmin,xmax,ymin,ymax,dmax,bin,list)
c
c     Compute Delaunay triangulation
c
      call delaun (numpts,n,x,list,bin,v,e,numtri)
c
c     Reset (x,y) coordinates to original values
c
      do 20 i = 1, n
         p    = list(i)
         x(1,p) = x(1,p) * dmax + xmin
         x(2,p) = x(2,p) * dmax + ymin
   20 continue
c
      return
      end

c**********************************************************************c

      subroutine delaun(numpts,n,x,list,stack,v,e,numtri)
      implicit double precision (a-h,o-z)
c
      integer v(3,*),n,i,t,list(*),numtri,p,e(3,*),maxstk,topstk,
     +        v1,v2,v3,l,r,pop,a,b,c,erl,era,erb,edg,triloc,numpts,
     +        tstrt,tstop,stack(*)
c
      real*8  x(2,*),xp,yp,c00000,c00100
c
      logical swap
c
      parameter (c00000=0.0,
     +           c00100=100.0)
c
c     Define vertex and adjacent lists for supertriangle
c
      v1 = numpts + 1
      v2 = numpts + 2
      v3 = numpts + 3
      v(1,1) = v1
      v(2,1) = v2
      v(3,1) = v3
      e(1,1) = 0
      e(2,1) = 0
      e(3,1) = 0
c
c     Set coordinates of supertriangle
c
      x(1,v1) = -c00100
      x(1,v2) =  c00100
      x(1,v3) =  c00000
      x(2,v1) = -c00100
      x(2,v2) = -c00100
      x(2,v3) =  c00100
c
c     Loop over each point
c
      numtri = 1
      topstk = 0
      maxstk = numpts
      do 100 i = 1, n
         p  = list(i)
         xp = x(1,p)
         yp = x(2,p)
c
c     Locate triangle in which point lies
c
      t = triloc(xp,yp,x,v,e,numtri)
c
c     Create new vertex and adjacency lists for triangle t
c
      a = e(1,t)
      b = e(2,t)
      c = e(3,t)
      v1 = v(1,t)
      v2 = v(2,t)
      v3 = v(3,t)
      v(1,t) = p
      v(2,t) = v1
      v(3,t) = v2
      e(1,t) = numtri+2
      e(2,t) = a
      e(3,t) = numtri+1
c
c     Create new triangles
c
      numtri = numtri + 1
      v(1,numtri) = p
      v(2,numtri) = v2
      v(3,numtri) = v3
      e(1,numtri) = t
      e(2,numtri) = b
      e(3,numtri) = numtri+1
      numtri = numtri + 1
      v(1,numtri) = p
      v(2,numtri) = v3
      v(3,numtri) = v1
      e(1,numtri) = numtri - 1
      e(2,numtri) = c
      e(3,numtri) = t
c
c     Put each edge of triangle t on stack.
c     Store triangles on left side of each edge.
c     Update adjacent lists for adjacent triangles.
c     Adjacency list for element a does not need to be updated.
c
      if (a.ne.0) then
         call pushdel(t,maxstk,topstk,stack)
      endif
      if (b.ne.0) then
         e(edg(b,t,e),b) = numtri - 1
         call pushdel (numtri-1,maxstk,topstk,stack)
      endif
      if (c.ne.0) then
         e(edg(c,t,e),c) = numtri
         call pushdel (numtri,maxstk,topstk,stack)
      endif
c
c     Loop while stack is not empty
c
   50 if (topstk.gt.0) then
         l = pop(topstk,stack)
         r = e(2,l)
c
c     Check if new point is in circumcircle for triangle r
c
         erl = edg(r,l,e)
         era = mod(erl,3) + 1
         erb = mod(era,3) + 1
         v1  = v(erl,r)
         v2  = v(era,r)
         v3  = v(erb,r)
         if 
     1   (swap(x(1,v1),x(2,v1),x(1,v2),x(2,v2),x(1,v3),x(2,v3),xp,yp)) 
     2   then
c
c     New point is inside circumcircle for triangle r
c     swap diagonal for convex quad formed by p-v2-v3-v1
c
            a = e(era,r)
            b = e(erb,r)
            c = e(3,l)
c
c     Update vertex and adjacency list for triangle l
c
            v(3,l) = v3
            e(2,l) = a
            e(3,l) = r
c
c     Update vertex and adjacency list for triangle r
c    
            v(1,r) = p
            v(2,r) = v3
            v(3,r) = v1
            e(1,r) = l
            e(2,r) = b
            e(3,r) = c
c
c     Put edges l-a and r-b on stack
c     Update adjacent lists for triangles a and c
c
            if (a.ne.0) then
               e(edg(a,r,e),a) = l
               call pushdel(l,maxstk,topstk,stack)
            endif
            if (b.ne.0) then
               call pushdel(r,maxstk,topstk,stack)
            endif
            if (c.ne.0) then
               e(edg(c,l,e),c) = r
            endif
         endif
         go to 50
      endif
  100 continue
c
c     Check consistency of triangulation
c
      if (numtri.ne.2*n+1) then
         write(6,1000)
         write(6,2000)
         stop
      endif
c
c     Remove all triangles containing supertriangle vertices
c     Find first triangle to be deleted (triangle t)
c     Update adjacency lists for triangles adjacent to t
c
      do 120 t = 1, numtri
         if ((v(1,t).gt.numpts) .or.
     +       (v(2,t).gt.numpts) .or.
     +       (v(3,t).gt.numpts)) then
            do 110 i = 1, 3
               a = e(i,t)
               if (a.ne.0) then
                  e(edg(a,t,e),a) = 0
               endif
  110       continue
            go to 125
         endif
  120 continue
  125 tstrt  = t + 1
      tstop  = numtri
      numtri = t - 1
c
c     Remove triangles
c
      do 200 t = tstrt, tstop
         if ((v(1,t).gt.numpts) .or.
     +       (v(2,t).gt.numpts) .or.
     +       (v(3,t).gt.numpts)) then
c
c     Triangle t is to be deleted
c     Updated adjaceny list for triangles adjacent to t
c
      do 130 i = 1, 3
         a = e(i,t)
         if (a.ne.0) then
            e(edg(a,t,e),a) = 0
         endif
  130 continue
      else
c
c     Triangle t is not to be deleted
c     Put triangle t in place of triangle numtri
c     Update adjaceny lists for triangles adjacent to t
c
         numtri = numtri + 1
         do 140 i = 1, 3
            a = e(i,t)
            e(i,numtri) = a
            v(i,numtri) = v(i,t)
            if (a.ne.0) then
               e(edg(a,t,e),a) = numtri
            endif
  140    continue
      endif
  200 continue
c
 1000 format(' **** Error in subroutine delaun ****')
 2000 format(' **** Incorrect number of triangles formed ****')
c
      return
      end

c**********************************************************************c

      function edg (l,k,e)
      implicit double precision (a-h,o-z)
c
      integer l,k,i,e(3,*),edg
c
      do 10 i = 1, 3
         if (e(i,l).eq.k) then
            edg = i
            return
      endif
   10 continue
c
      write(6,1000)
      write(6,2000)
c
 1000 format(' **** Error in function edg ****')
 2000 format(' **** Elements not adjacent ****')
      stop
c
      end

c**********************************************************************c

      function pop (topstk,stack)
      implicit double precision (a-h,o-z)
c
      integer pop,topstk,stack(*)
c
      if (topstk.gt.0) then
         pop = stack(topstk)
         topstk = topstk - 1
      else
         write(6,1000)
         write(6,2000)
         stop
      endif
c
 1000 format(' **** Error in function pop')
 2000 format(' **** Stack underflow')
c
      return
      end

c**********************************************************************c

      subroutine pushdel(item,maxstk,topstk,stack)
      implicit double precision (a-h,o-z)
c
      integer topstk,maxstk,stack(*),item
c
      topstk = topstk + 1
      if (topstk.gt.maxstk) then
         write(6,1000)
         write(6,2000)
         stop
      else
         stack(topstk) = item
      endif
c
 1000 format(' **** Error in subroutine pushdel ****')
 2000 format(' **** Stack overflow ****')
c
      return
      end

c**********************************************************************c

      subroutine qsorti(n,list,key)
      implicit double precision (a-h,o-z)
c
      integer list(*),key(*),n,ll,lr,lm,nl,nr,ltemp,stktop,maxstk,guess
c
      parameter (maxstk=32)
c
      integer lstack(maxstk), rstack(maxstk)
c
      ll = 1
      lr = n
      stktop = 0
   10 if (ll.lt.lr) then
         nl = ll
         nr = lr
         lm = (ll+lr) / 2
         guess = key(list(lm))
c
c     Find keys for exchange
c
   20    if (key(list(nl)).lt.guess) then
            nl = nl + 1
            go to 20
         endif
   30    if (guess.lt.key(list(nr))) then
            nr = nr - 1
            go to 30
         endif
         if (nl.lt.(nr-1)) then
            ltemp    = list(nl)
            list(nl) = list(nr)
            list(nr) = ltemp
            nl = nl + 1
            nr = nr - 1
            go to 20
         endif
c
c     Deal with crossing of pointers
c
         if (nl.le.nr) then
            if (nl.lt.nr) then
               ltemp    = list(nl)
               list(nl) = list(nr)
               list(nr) = ltemp
            endif
            nl = nl + 1
            nr = nr - 1
         endif
c
c     Select sub-list to be processed next
c
         stktop = stktop + 1
         if (nr.lt.lm) then
            lstack(stktop) = nl
            rstack(stktop) = lr
            lr = nr
         else
            lstack(stktop) = ll
            rstack(stktop) = nr
            ll = nl
         endif
         go to 10
      endif 
c
c     Process any stacked sub-lists
c
      if (stktop.ne.0) then
         ll = lstack(stktop)
         lr = rstack(stktop)
         stktop = stktop - 1
         go to 10
      endif
c
      return
      end

c**********************************************************************c

      function swap (x1,y1,x2,y2,x3,y3,xp,yp)
      implicit double precision (a-h,o-z)
c
      real*8 x1,y1,x2,y2,x3,y3,xp,yp,x13,y13,x23,y23,x1p,y1p,
     1 x2p,y2p,cosa,cosb,sina,sinb,c00000
c
      logical swap
c
      parameter (c00000=0.0)
c
      x13 = x1 - x3
      y13 = y1 - y3
      x23 = x2 - x3
      y23 = y2 - y3
      x1p = x1 - xp
      y1p = y1 - yp
      x2p = x2 - xp
      y2p = y2 - yp
c
      cosa = x13 * x23 + y13 * y23
      cosb = x2p * x1p + y1p * y2p
c
      if ((cosa.ge.c00000) .and. (cosb.ge.c00000)) then
         swap = .false.
      elseif ((cosa.lt.c00000) .and. (cosb.lt.c00000)) then
         swap = .true.
      else
         sina = x13 * y23 - x23 * y13
         sinb = x2p * y1p - x1p * y2p
         if ((sina*cosb+sinb*cosa).lt.c00000) then
            swap = .true.
         else
            swap = .false.
         endif
      endif
c
      end

c**********************************************************************c

      function triloc (xp,yp,x,v,e,numtri)
      implicit double precision (a-h,o-z)
c
      integer v(3,*),e(3,*),numtri,v1,v2,i,t,triloc
c
      real*8  x(2,*),xp,yp
c
      t = numtri
   10 continue
      do 20 i = 1, 3
         v1 = v(i,t)
         v2 = v(mod(i,3)+1,t)
         if 
     1   ((x(2,v1)-yp)*(x(1,v2)-xp).gt.(x(1,v1)-xp)*(x(2,v2)-yp)) 
     2   then
            t = e(i,t)
            go to 10
         endif
   20 continue
c
c     Triangle has been found
c
      triloc = t
c
      end

c**********************************************************************c
c**********************************************************************c
c**********************************************************************c
      logical function inside(x,xtri)
      implicit double precision (a-h,o-z)
c
      dimension x(2),xtri(2,1)
      data tol/1.d-7/
c
      det = xtri(1,2)*xtri(2,3) - xtri(1,3)*xtri(2,2)
     1    + xtri(1,3)*xtri(2,1) - xtri(1,1)*xtri(2,3)
     2    + xtri(1,1)*xtri(2,2) - xtri(1,2)*xtri(2,1)
c
      p1  =(xtri(1,2)*xtri(2,3) - xtri(1,3)*xtri(2,2)
     1    + xtri(1,3)*x(2)      - x(1)*xtri(2,3)
     2    + x(1)*xtri(2,2)      - xtri(1,2)*x(2))/det
c
      p2  =(x(1)*xtri(2,3)      - xtri(1,3)*x(2)
     1    + xtri(1,3)*xtri(2,1) - xtri(1,1)*xtri(2,3)
     2    + xtri(1,1)*x(2)      - x(1)*xtri(2,1))/det
c
      p3  =(xtri(1,2)*x(2)      - x(1)*xtri(2,2)
     1    + x(1)*xtri(2,1)      - xtri(1,1)*x(2)
     2    + xtri(1,1)*xtri(2,2) - xtri(1,2)*xtri(2,1))/det
c
      pmin = - tol
      pmax = 1.d0 + tol
      if ((p1.ge.pmin).and.(p1.le.pmax).and.
     1    (p2.ge.pmin).and.(p2.le.pmax).and.
     2    (p3.ge.pmin).and.(p3.le.pmax)) then
      inside = .true.
      else
      inside = .false.
      end if
      return
      end
c**********************************************************************c
      subroutine bsort(n,x,xmin,xmax,ymin,ymax,dmax,bin,list)
      implicit double precision (a-h,o-z)
c
      integer list(*),bin(*),n,i,j,k,p,ndiv
c
      real*8 x(2,*),factx,facty,xmin,xmax,ymin,ymax,dmax
c
c     Compute number of bins in (x,y) coordinates
c     Compute inverse of bin size in (x,y) directions
c
      ndiv  = nint(real(n)**0.25)
      factx = real(ndiv) / ((xmax-xmin)*1.01/dmax)
      facty = real(ndiv) / ((ymax-ymin)*1.01/dmax)
c
c     Assign bin numbers to each point
c
      do 10 k = 1, n
         p = list(k)
         i = int(x(2,p)*facty)
         j = int(x(1,p)*factx)
         if (mod(i,2).eq.0) then
            bin(p) = i*ndiv+j+1
         else
            bin(p) = (i+1)*ndiv-j
         endif
   10 continue
c
c     Sort points in ascending sequence of bin number
c
      call qsorti (n,list,bin)
c
      end

