! Copyright 2019
!
! Zdenek Masin with contributions from others (see the UK-AMOR website)                               
!
! This file is part of GBTOlib.
!
!     GBTOlib is free software: you can redistribute it and/or modify
!     it under the terms of the GNU General Public License as published by
!     the Free Software Foundation, either version 3 of the License, or
!     (at your option) any later version.
!
!     GBTOlib is distributed in the hope that it will be useful,
!     but WITHOUT ANY WARRANTY; without even the implied warranty of
!     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!     GNU General Public License for more details.
!
!     You should have received a copy of the GNU General Public License
!     along with  GBTOlib (in trunk/COPYING). Alternatively, you can also visit
!     <https://www.gnu.org/licenses/>.
!
module sort_gbl
  use precisn_gbl
  use utils_gbl, only: xermsg

  implicit none

  private

  !> \todo interface the rest of the routines in this module
  interface cfp_sort_float_int_1d
     module procedure wp_sort_float_int_1d, ep1_sort_float_int_1d
  end interface

  interface sort_int_1d
     module procedure sort_int_1d_int32, sort_int_1d_int64
  end interface

  public multiway_merge_float_int

  !> todo OpenMP parallelize the routines below
  public sort_float, sort_int_1d, sort_int_float, sort_float_int, sort_int_int, cfp_sort_float_int_1d, heap_sort_int_float

contains

   !> Sorts the first n elements of the floating point array arr(:,d2).
   subroutine sort_float(n,d2,arr)
     implicit none
   
     integer :: n, d2
     real(kind=cfp), intent(inout) :: arr(:,:)
     integer, parameter :: m=7,nstack=50
     integer :: i,ir,j,jstack,k,l,istack(nstack)
     real(kind=cfp) :: a,temp
        jstack=0
        l=1
        ir=n
   1      if(ir-l.lt.m)then
          do j=l+1,ir
            a=arr(j,d2)
            do i=j-1,1,-1
              if(arr(i,d2).le.a)goto 2
              arr(i+1,d2)=arr(i,d2)
            end do
            i=0
   2        arr(i+1,d2)=a
          end do
          if(jstack.eq.0)return
          ir=istack(jstack)
          l=istack(jstack-1)
          jstack=jstack-2
        else
          k=(l+ir)/2
          temp=arr(k,d2)
          arr(k,d2)=arr(l+1,d2)
          arr(l+1,d2)=temp
          if(arr(l+1,d2).gt.arr(ir,d2))then
            temp=arr(l+1,d2)
            arr(l+1,d2)=arr(ir,d2)
            arr(ir,d2)=temp
          endif
          if(arr(l,d2).gt.arr(ir,d2))then
            temp=arr(l,d2)
            arr(l,d2)=arr(ir,d2)
            arr(ir,d2)=temp
          endif
          if(arr(l+1,d2).gt.arr(l,d2))then
            temp=arr(l+1,d2)
            arr(l+1,d2)=arr(l,d2)
            arr(l,d2)=temp
          endif
          i=l+1
          j=ir
          a=arr(l,d2)
   3      continue
            i=i+1
          if(arr(i,d2).lt.a)goto 3
   4      continue
            j=j-1
          if(arr(j,d2).gt.a)goto 4
          if(j.lt.i)goto 5
          temp=arr(i,d2)
          arr(i,d2)=arr(j,d2)
          arr(j,d2)=temp
          goto 3
   5      arr(l,d2)=arr(j,d2)
          arr(j,d2)=a
          jstack=jstack+2
          if(jstack.gt.nstack) call xermsg('sort','sort_float','nstack parameter too small.',1,1)
          if(ir-i+1.ge.j-l)then
            istack(jstack)=ir
            istack(jstack-1)=i
            ir=j-1
          else
            istack(jstack)=j-1
            istack(jstack-1)=l
            l=i
          endif
        endif
        goto 1
   end subroutine sort_float

   !> Sorts the first n elements of the 1D integer array arr(:).
   subroutine sort_int_1d_int32(n,arr)
     implicit none
   
     integer(shortint) :: n
     integer(shortint), intent(inout) :: arr(:)
     integer, parameter :: m=7,nstack=50
     integer :: i,ir,j,jstack,k,l,istack(nstack)
     integer :: a,temp
        jstack=0
        l=1
        ir=n
   1      if(ir-l.lt.m)then
          do j=l+1,ir
            a=arr(j)
            do i=j-1,1,-1
              if(arr(i).le.a)goto 2
              arr(i+1)=arr(i)
            end do
            i=0
   2        arr(i+1)=a
          end do
          if(jstack.eq.0)return
          ir=istack(jstack)
          l=istack(jstack-1)
          jstack=jstack-2
        else
          k=(l+ir)/2
          temp=arr(k)
          arr(k)=arr(l+1)
          arr(l+1)=temp
          if(arr(l+1).gt.arr(ir))then
            temp=arr(l+1)
            arr(l+1)=arr(ir)
            arr(ir)=temp
          endif
          if(arr(l).gt.arr(ir))then
            temp=arr(l)
            arr(l)=arr(ir)
            arr(ir)=temp
          endif
          if(arr(l+1).gt.arr(l))then
            temp=arr(l+1)
            arr(l+1)=arr(l)
            arr(l)=temp
          endif
          i=l+1
          j=ir
          a=arr(l)
   3      continue
            i=i+1
          if(arr(i).lt.a)goto 3
   4      continue
            j=j-1
          if(arr(j).gt.a)goto 4
          if(j.lt.i)goto 5
          temp=arr(i)
          arr(i)=arr(j)
          arr(j)=temp
          goto 3
   5      arr(l)=arr(j)
          arr(j)=a
          jstack=jstack+2
          if(jstack.gt.nstack) call xermsg('sort','sort_int_1d_shortint','nstack parameter too small.',1,1)
          if(ir-i+1.ge.j-l)then
            istack(jstack)=ir
            istack(jstack-1)=i
            ir=j-1
          else
            istack(jstack)=j-1
            istack(jstack-1)=l
            l=i
          endif
        endif
        goto 1
   end subroutine sort_int_1d_int32

   !> Sorts the first n elements of the 1D integer array arr(:).
   subroutine sort_int_1d_int64(n,arr)
     implicit none
   
     integer(longint) :: n
     integer(longint), intent(inout) :: arr(:)
     integer, parameter :: m=7,nstack=50
     integer :: i,ir,j,jstack,k,l,istack(nstack)
     integer :: a,temp
        jstack=0
        l=1
        ir=n
   1      if(ir-l.lt.m)then
          do j=l+1,ir
            a=arr(j)
            do i=j-1,1,-1
              if(arr(i).le.a)goto 2
              arr(i+1)=arr(i)
            end do
            i=0
   2        arr(i+1)=a
          end do
          if(jstack.eq.0)return
          ir=istack(jstack)
          l=istack(jstack-1)
          jstack=jstack-2
        else
          k=(l+ir)/2
          temp=arr(k)
          arr(k)=arr(l+1)
          arr(l+1)=temp
          if(arr(l+1).gt.arr(ir))then
            temp=arr(l+1)
            arr(l+1)=arr(ir)
            arr(ir)=temp
          endif
          if(arr(l).gt.arr(ir))then
            temp=arr(l)
            arr(l)=arr(ir)
            arr(ir)=temp
          endif
          if(arr(l+1).gt.arr(l))then
            temp=arr(l+1)
            arr(l+1)=arr(l)
            arr(l)=temp
          endif
          i=l+1
          j=ir
          a=arr(l)
   3      continue
            i=i+1
          if(arr(i).lt.a)goto 3
   4      continue
            j=j-1
          if(arr(j).gt.a)goto 4
          if(j.lt.i)goto 5
          temp=arr(i)
          arr(i)=arr(j)
          arr(j)=temp
          goto 3
   5      arr(l)=arr(j)
          arr(j)=a
          jstack=jstack+2
          if(jstack.gt.nstack) call xermsg('sort','sort_int_1d_longint','nstack parameter too small.',1,1)
          if(ir-i+1.ge.j-l)then
            istack(jstack)=ir
            istack(jstack-1)=i
            ir=j-1
          else
            istack(jstack)=j-1
            istack(jstack-1)=l
            l=i
          endif
        endif
        goto 1
   end subroutine sort_int_1d_int64

   !> Heap-sorts the first n elements of the array arr and simultaneously the array brr. Array arr(:) is the list of integers and the array brr(:) is the list of floats.
   subroutine heap_sort_int_float(n,arr,brr)
     implicit none
     integer :: n
     integer, intent(inout) :: arr(:)
     real(kind=cfp), intent(inout) :: brr(:)

     integer :: i,ir,j,l, arr_tmp
     real(kind=cfp) :: brr_tmp
      
        if (n.lt.2) return
        
        l=n/2+1
        ir=n
        
  10    continue
        
        if(l.gt.1)then
        
           l=l-1
           
           arr_tmp=arr(l)
           brr_tmp=brr(l)
        
        else
        
           arr_tmp=arr(ir)
           brr_tmp=brr(ir)
           
           arr(ir)=arr(1)
           brr(ir)=brr(1)
           
           ir=ir-1
           
           if(ir.eq.1)then
           
              arr(1)=arr_tmp
              brr(1)=brr_tmp
              
              return
           
           endif
        
        endif
        
        i=l
        
        j=l+l
        
  20    if(j.le.ir)then
        
           if(j.lt.ir)then
           
              if(arr(j).lt.arr(j+1))j=j+1
           
           endif
        
           if(arr_tmp.lt.arr(j))then
           
              arr(i)=arr(j)
              brr(i)=brr(j)
              
              i=j
              
              j=j+j
           
           else
           
              j=ir+1
           
           endif
        
           goto 20
        
        endif
        
        arr(i)=arr_tmp
        brr(i)=brr_tmp
        
        goto 10

   end subroutine heap_sort_int_float

   !> Sorts the first n elements of the array arr and simultaneously the array brr. Array arr(:,d2) is the list of integers and the array brr(:,d2) is the list of floats.
   subroutine sort_int_float(n,d2,arr,brr)
     implicit none
   
     integer :: n, d2
     integer, intent(inout) :: arr(:,:)
     real(kind=cfp), intent(inout) :: brr(:,:)
     integer, parameter :: m=7,nstack=50
     integer :: i,ir,j,jstack,k,l,istack(nstack)
     real(kind=cfp) :: a,b,temp
        jstack=0
        l=1
        ir=n
   1      if(ir-l.lt.m)then
          do j=l+1,ir
            a=arr(j,d2)
            b=brr(j,d2)
            do i=j-1,1,-1
              if(arr(i,d2).le.a)goto 2
              arr(i+1,d2)=arr(i,d2)
              brr(i+1,d2)=brr(i,d2)
            end do
            i=0
   2        arr(i+1,d2)=a
            brr(i+1,d2)=b
          end do
          if(jstack.eq.0)return
          ir=istack(jstack)
          l=istack(jstack-1)
          jstack=jstack-2
        else
          k=(l+ir)/2
          temp=arr(k,d2)
          arr(k,d2)=arr(l+1,d2)
          arr(l+1,d2)=temp
          temp=brr(k,d2)
          brr(k,d2)=brr(l+1,d2)
          brr(l+1,d2)=temp
          if(arr(l+1,d2).gt.arr(ir,d2))then
            temp=arr(l+1,d2)
            arr(l+1,d2)=arr(ir,d2)
            arr(ir,d2)=temp
            temp=brr(l+1,d2)
            brr(l+1,d2)=brr(ir,d2)
            brr(ir,d2)=temp
          endif
          if(arr(l,d2).gt.arr(ir,d2))then
            temp=arr(l,d2)
            arr(l,d2)=arr(ir,d2)
            arr(ir,d2)=temp
            temp=brr(l,d2)
            brr(l,d2)=brr(ir,d2)
            brr(ir,d2)=temp
          endif
          if(arr(l+1,d2).gt.arr(l,d2))then
            temp=arr(l+1,d2)
            arr(l+1,d2)=arr(l,d2)
            arr(l,d2)=temp
            temp=brr(l+1,d2)
            brr(l+1,d2)=brr(l,d2)
            brr(l,d2)=temp
          endif
          i=l+1
          j=ir
          a=arr(l,d2)
          b=brr(l,d2)
   3      continue
            i=i+1
          if(arr(i,d2).lt.a)goto 3
   4      continue
            j=j-1
          if(arr(j,d2).gt.a)goto 4
          if(j.lt.i)goto 5
          temp=arr(i,d2)
          arr(i,d2)=arr(j,d2)
          arr(j,d2)=temp
          temp=brr(i,d2)
          brr(i,d2)=brr(j,d2)
          brr(j,d2)=temp
          goto 3
   5      arr(l,d2)=arr(j,d2)
          arr(j,d2)=a
          brr(l,d2)=brr(j,d2)
          brr(j,d2)=b
          jstack=jstack+2
          if(jstack.gt.nstack) call xermsg('sort','sort_int_float','nstack parameter too small.',1,1)
          if(ir-i+1.ge.j-l)then
            istack(jstack)=ir
            istack(jstack-1)=i
            ir=j-1
          else
            istack(jstack)=j-1
            istack(jstack-1)=l
            l=i
          endif
        endif
        goto 1
   end subroutine sort_int_float

   !> Sorts the first n elements of the array arr and simultaneously the array brr. Array arr(:,d2) is the list of floats and the array brr(:,d2) is the list of integers.
   subroutine sort_float_int(n,d2,arr,brr)
     implicit none
   
     integer :: n, d2
     integer, intent(inout) :: brr(:,:)
     real(kind=cfp), intent(inout) :: arr(:,:)
     integer, parameter :: m=7,nstack=50
     integer :: i,ir,j,jstack,k,l,istack(nstack)
     real(kind=cfp) :: a,b,temp
        jstack=0
        l=1
        ir=n
   1      if(ir-l.lt.m)then
          do j=l+1,ir
            a=arr(j,d2)
            b=brr(j,d2)
            do i=j-1,1,-1
              if(arr(i,d2).le.a)goto 2
              arr(i+1,d2)=arr(i,d2)
              brr(i+1,d2)=brr(i,d2)
            end do
            i=0
   2        arr(i+1,d2)=a
            brr(i+1,d2)=b
          end do
          if(jstack.eq.0)return
          ir=istack(jstack)
          l=istack(jstack-1)
          jstack=jstack-2
        else
          k=(l+ir)/2
          temp=arr(k,d2)
          arr(k,d2)=arr(l+1,d2)
          arr(l+1,d2)=temp
          temp=brr(k,d2)
          brr(k,d2)=brr(l+1,d2)
          brr(l+1,d2)=temp
          if(arr(l+1,d2).gt.arr(ir,d2))then
            temp=arr(l+1,d2)
            arr(l+1,d2)=arr(ir,d2)
            arr(ir,d2)=temp
            temp=brr(l+1,d2)
            brr(l+1,d2)=brr(ir,d2)
            brr(ir,d2)=temp
          endif
          if(arr(l,d2).gt.arr(ir,d2))then
            temp=arr(l,d2)
            arr(l,d2)=arr(ir,d2)
            arr(ir,d2)=temp
            temp=brr(l,d2)
            brr(l,d2)=brr(ir,d2)
            brr(ir,d2)=temp
          endif
          if(arr(l+1,d2).gt.arr(l,d2))then
            temp=arr(l+1,d2)
            arr(l+1,d2)=arr(l,d2)
            arr(l,d2)=temp
            temp=brr(l+1,d2)
            brr(l+1,d2)=brr(l,d2)
            brr(l,d2)=temp
          endif
          i=l+1
          j=ir
          a=arr(l,d2)
          b=brr(l,d2)
   3      continue
            i=i+1
          if(arr(i,d2).lt.a)goto 3
   4      continue
            j=j-1
          if(arr(j,d2).gt.a)goto 4
          if(j.lt.i)goto 5
          temp=arr(i,d2)
          arr(i,d2)=arr(j,d2)
          arr(j,d2)=temp
          temp=brr(i,d2)
          brr(i,d2)=brr(j,d2)
          brr(j,d2)=temp
          goto 3
   5      arr(l,d2)=arr(j,d2)
          arr(j,d2)=a
          brr(l,d2)=brr(j,d2)
          brr(j,d2)=b
          jstack=jstack+2
          if(jstack.gt.nstack) call xermsg('sort','sort_int_float','nstack parameter too small.',1,1)
          if(ir-i+1.ge.j-l)then
            istack(jstack)=ir
            istack(jstack-1)=i
            ir=j-1
          else
            istack(jstack)=j-1
            istack(jstack-1)=l
            l=i
          endif
        endif
        goto 1
   end subroutine sort_float_int

   !> Sorts the first n elements of the array arr and simultaneously the array brr. Array arr(:,d2) is the list of integers and the array brr(:,d2) is the list of integers.
   subroutine sort_int_int(n,d2,arr,brr)
     implicit none
   
     integer :: n, d2
     integer, intent(inout) :: arr(:,:)
     integer, intent(inout) :: brr(:,:)
     integer, parameter :: m=7,nstack=50
     integer :: i,ir,j,jstack,k,l,istack(nstack)
     integer :: a,b,temp
        jstack=0
        l=1
        ir=n
   1      if(ir-l.lt.m)then
          do j=l+1,ir
            a=arr(j,d2)
            b=brr(j,d2)
            do i=j-1,1,-1
              if(arr(i,d2).le.a)goto 2
              arr(i+1,d2)=arr(i,d2)
              brr(i+1,d2)=brr(i,d2)
            end do
            i=0
   2        arr(i+1,d2)=a
            brr(i+1,d2)=b
          end do
          if(jstack.eq.0)return
          ir=istack(jstack)
          l=istack(jstack-1)
          jstack=jstack-2
        else
          k=(l+ir)/2
          temp=arr(k,d2)
          arr(k,d2)=arr(l+1,d2)
          arr(l+1,d2)=temp
          temp=brr(k,d2)
          brr(k,d2)=brr(l+1,d2)
          brr(l+1,d2)=temp
          if(arr(l+1,d2).gt.arr(ir,d2))then
            temp=arr(l+1,d2)
            arr(l+1,d2)=arr(ir,d2)
            arr(ir,d2)=temp
            temp=brr(l+1,d2)
            brr(l+1,d2)=brr(ir,d2)
            brr(ir,d2)=temp
          endif
          if(arr(l,d2).gt.arr(ir,d2))then
            temp=arr(l,d2)
            arr(l,d2)=arr(ir,d2)
            arr(ir,d2)=temp
            temp=brr(l,d2)
            brr(l,d2)=brr(ir,d2)
            brr(ir,d2)=temp
          endif
          if(arr(l+1,d2).gt.arr(l,d2))then
            temp=arr(l+1,d2)
            arr(l+1,d2)=arr(l,d2)
            arr(l,d2)=temp
            temp=brr(l+1,d2)
            brr(l+1,d2)=brr(l,d2)
            brr(l,d2)=temp
          endif
          i=l+1
          j=ir
          a=arr(l,d2)
          b=brr(l,d2)
   3      continue
            i=i+1
          if(arr(i,d2).lt.a)goto 3
   4      continue
            j=j-1
          if(arr(j,d2).gt.a)goto 4
          if(j.lt.i)goto 5
          temp=arr(i,d2)
          arr(i,d2)=arr(j,d2)
          arr(j,d2)=temp
          temp=brr(i,d2)
          brr(i,d2)=brr(j,d2)
          brr(j,d2)=temp
          goto 3
   5      arr(l,d2)=arr(j,d2)
          arr(j,d2)=a
          brr(l,d2)=brr(j,d2)
          brr(j,d2)=b
          jstack=jstack+2
          if(jstack.gt.nstack) call xermsg('sort','sort_int_float','nstack parameter too small.',1,1)
          if(ir-i+1.ge.j-l)then
            istack(jstack)=ir
            istack(jstack-1)=i
            ir=j-1
          else
            istack(jstack)=j-1
            istack(jstack-1)=l
            l=i
          endif
        endif
        goto 1
   end subroutine sort_int_int

   !> Sorts the first n elements of the array arr and simultaneously the array brr. Array arr(:) is the list of floats and the array brr(:) is the list of integers.
   subroutine wp_sort_float_int_1d(n,arr,brr)
     implicit none
   
     integer :: n
     integer, intent(inout) :: brr(:)
     real(kind=wp), intent(inout) :: arr(:)
     integer, parameter :: m=7,nstack=50
     integer :: i,ir,j,jstack,k,l,istack(nstack)
     real(kind=wp) :: a,b,temp
        jstack=0
        l=1
        ir=n
   1      if(ir-l.lt.m)then
          do j=l+1,ir
            a=arr(j)
            b=brr(j)
            do i=j-1,1,-1
              if(arr(i).le.a)goto 2
              arr(i+1)=arr(i)
              brr(i+1)=brr(i)
            end do
            i=0
   2        arr(i+1)=a
            brr(i+1)=b
          end do
          if(jstack.eq.0)return
          ir=istack(jstack)
          l=istack(jstack-1)
          jstack=jstack-2
        else
          k=(l+ir)/2
          temp=arr(k)
          arr(k)=arr(l+1)
          arr(l+1)=temp
          temp=brr(k)
          brr(k)=brr(l+1)
          brr(l+1)=temp
          if(arr(l+1).gt.arr(ir))then
            temp=arr(l+1)
            arr(l+1)=arr(ir)
            arr(ir)=temp
            temp=brr(l+1)
            brr(l+1)=brr(ir)
            brr(ir)=temp
          endif
          if(arr(l).gt.arr(ir))then
            temp=arr(l)
            arr(l)=arr(ir)
            arr(ir)=temp
            temp=brr(l)
            brr(l)=brr(ir)
            brr(ir)=temp
          endif
          if(arr(l+1).gt.arr(l))then
            temp=arr(l+1)
            arr(l+1)=arr(l)
            arr(l)=temp
            temp=brr(l+1)
            brr(l+1)=brr(l)
            brr(l)=temp
          endif
          i=l+1
          j=ir
          a=arr(l)
          b=brr(l)
   3      continue
            i=i+1
          if(arr(i).lt.a)goto 3
   4      continue
            j=j-1
          if(arr(j).gt.a)goto 4
          if(j.lt.i)goto 5
          temp=arr(i)
          arr(i)=arr(j)
          arr(j)=temp
          temp=brr(i)
          brr(i)=brr(j)
          brr(j)=temp
          goto 3
   5      arr(l)=arr(j)
          arr(j)=a
          brr(l)=brr(j)
          brr(j)=b
          jstack=jstack+2
          if(jstack.gt.nstack) call xermsg('sort','sort_int_float','nstack parameter too small.',1,1)
          if(ir-i+1.ge.j-l)then
            istack(jstack)=ir
            istack(jstack-1)=i
            ir=j-1
          else
            istack(jstack)=j-1
            istack(jstack-1)=l
            l=i
          endif
        endif
        goto 1
   end subroutine wp_sort_float_int_1d

   !> Sorts the first n elements of the array arr and simultaneously the array brr. Array arr(:) is the list of floats and the array brr(:) is the list of integers.
   subroutine ep1_sort_float_int_1d(n,arr,brr)
     implicit none
   
     integer :: n
     integer, intent(inout) :: brr(:)
     real(kind=ep1), intent(inout) :: arr(:)
     integer, parameter :: m=7,nstack=50
     integer :: i,ir,j,jstack,k,l,istack(nstack)
     real(kind=ep1) :: a,b,temp
        jstack=0
        l=1
        ir=n
   1      if(ir-l.lt.m)then
          do j=l+1,ir
            a=arr(j)
            b=brr(j)
            do i=j-1,1,-1
              if(arr(i).le.a)goto 2
              arr(i+1)=arr(i)
              brr(i+1)=brr(i)
            end do
            i=0
   2        arr(i+1)=a
            brr(i+1)=b
          end do
          if(jstack.eq.0)return
          ir=istack(jstack)
          l=istack(jstack-1)
          jstack=jstack-2
        else
          k=(l+ir)/2
          temp=arr(k)
          arr(k)=arr(l+1)
          arr(l+1)=temp
          temp=brr(k)
          brr(k)=brr(l+1)
          brr(l+1)=temp
          if(arr(l+1).gt.arr(ir))then
            temp=arr(l+1)
            arr(l+1)=arr(ir)
            arr(ir)=temp
            temp=brr(l+1)
            brr(l+1)=brr(ir)
            brr(ir)=temp
          endif
          if(arr(l).gt.arr(ir))then
            temp=arr(l)
            arr(l)=arr(ir)
            arr(ir)=temp
            temp=brr(l)
            brr(l)=brr(ir)
            brr(ir)=temp
          endif
          if(arr(l+1).gt.arr(l))then
            temp=arr(l+1)
            arr(l+1)=arr(l)
            arr(l)=temp
            temp=brr(l+1)
            brr(l+1)=brr(l)
            brr(l)=temp
          endif
          i=l+1
          j=ir
          a=arr(l)
          b=brr(l)
   3      continue
            i=i+1
          if(arr(i).lt.a)goto 3
   4      continue
            j=j-1
          if(arr(j).gt.a)goto 4
          if(j.lt.i)goto 5
          temp=arr(i)
          arr(i)=arr(j)
          arr(j)=temp
          temp=brr(i)
          brr(i)=brr(j)
          brr(j)=temp
          goto 3
   5      arr(l)=arr(j)
          arr(j)=a
          brr(l)=brr(j)
          brr(j)=b
          jstack=jstack+2
          if(jstack.gt.nstack) call xermsg('sort','sort_int_float','nstack parameter too small.',1,1)
          if(ir-i+1.ge.j-l)then
            istack(jstack)=ir
            istack(jstack-1)=i
            ir=j-1
          else
            istack(jstack)=j-1
            istack(jstack-1)=l
            l=i
          endif
        endif
        goto 1
   end subroutine ep1_sort_float_int_1d


   !> \brief   Merge sorted arrays
   !> \authors J Benda
   !> \date    2025
   !>
   !> On input, this subroutine expect real and integer arrays `reals` and `ints` of length `n` that consist of
   !> consecutive chunks presorted (in acending way) on `ints`. The subroutine will make a copy of both these arrays
   !> and then perform multiway merge of the presorted chunks into the original storage.
   !>
   !> If non-empty `scratch_prefix_opt` is provided on input, the copy of the arrays will be mapped to disk prefixed
   !> with the given location. Anonymous mapping (i.e., plain allocation) is used otherwise.
   !>
   !> \param[in]    nelem   Length of the arrays `floats` and `ints`.
   !> \param[inout] reals   Array of reals.
   !> \param[inout] ints    Array if integers.
   !> \param[in]    sizes   Lengths of presorted chunks (according to `ints`).
   !> \param[in]    scratch_prefix_opt   Optional prefix of temporary scratch files.
   !>
   subroutine multiway_merge_float_int(nelem, reals, ints, sizes, scratch_prefix_opt)

      use file_mapping_gbl, only: file_mapping, int64
      use iso_c_binding,    only: c_f_pointer
      use utils_gbl,        only: xermsg

      integer,   intent(in)    :: nelem, sizes(:)
      integer,   intent(inout) :: ints(:)
      real(cfp), intent(inout) :: reals(:)

      character(len=*), optional, intent(in) :: scratch_prefix_opt

      type(file_mapping), target    :: reals_mapping, ints_mapping
      character(len=:), allocatable :: scratch

      integer, allocatable :: positions(:), ends(:)
      integer, pointer     :: ints0(:)
      real(cfp), pointer   :: reals0(:)
      integer              :: i, j, chunk, nchunks
      integer(int64)       :: offset, length
      logical              :: use_scratch

      ! early exit when there is nothing to merge
      if (size(sizes) < 2) return

      ! sanity check
      if (sum(sizes) /= nelem) then
         call xermsg('sort_gbl', 'multiway_merge_float_int', 'Incompatible input parameters!', 1, 1)
      end if

      ! find out if scratch use is requested
      if (present(scratch_prefix_opt)) then
         scratch = scratch_prefix_opt
         use_scratch = .true.
      else
         scratch = ''
         use_scratch = .false.
      end if

      ! set up mapping parameters (extend length by one for convenience of access in the merging below)
      offset = 0
      length = nelem

      ! allocate memory (on or off disk, depending on the input)
      call reals_mapping % init(scratch, offset, length*storage_size(reals)/8, .true., use_scratch, use_scratch)
      call ints_mapping % init(scratch, offset, length*storage_size(ints)/8, .true., use_scratch, use_scratch)
      call c_f_pointer(reals_mapping%ptr, reals0, [nelem])
      call c_f_pointer(ints_mapping%ptr, ints0, [nelem])

      ! make a copy of the input arrays into the mapping
      !$omp parallel do
      do i = 1, nelem
         reals0(i) = reals(i)
         ints0(i) = ints(i)
      end do

      ! prepare chunk information arrays
      nchunks = size(sizes)
      allocate (positions(nchunks), ends(nchunks))

      ! compute chunk information
      ends = sizes
      do i = 2, nchunks
         ends(i) = ends(i) + ends(i - 1)
      end do
      positions = ends - sizes + 1

      ! merge elements from all chunks
      do i = 1, nelem
         ! find non-depleted chunk that starts with the smallest integer
         chunk = 0
         do j = 1, nchunks
            if (positions(j) <= ends(j)) then
               if (chunk > 0) then
                  if (ints0(positions(chunk)) <= ints0(positions(j))) cycle
               end if
               chunk = j
            end if
         end do
         ! copy the data from the beginning of that chunk
         reals(i) = reals0(positions(chunk))
         ints(i) = ints0(positions(chunk))
         ! advance chunk pointer to the next element
         positions(chunk) = positions(chunk) + 1
      end do

   end subroutine multiway_merge_float_int

end module sort_gbl
