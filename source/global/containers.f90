! Copyright 2019
!
! For a comprehensive list of the developers that contributed to these codes
! see the UK-AMOR website.
!
! This file is part of UKRmol-in (UKRmol+ suite).
!
!     UKRmol-in is free software: you can redistribute it and/or modify
!     it under the terms of the GNU General Public License as published by
!     the Free Software Foundation, either version 3 of the License, or
!     (at your option) any later version.
!
!     UKRmol-in is distributed in the hope that it will be useful,
!     but WITHOUT ANY WARRANTY; without even the implied warranty of
!     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!     GNU General Public License for more details.
!
!     You should have received a copy of the GNU General Public License
!     along with  UKRmol-in (in source/COPYING). Alternatively, you can also visit
!     <https://www.gnu.org/licenses/>.

!> \brief   Standard containers
!> \authors J Benda
!> \date    2018 - 2019
!>
!> This module contains a special data structure, \ref bstree, used (e.g.) by CONGEN to carry out quick look-ups in the
!> storage of determinants.
!>
module containers

    use iso_fortran_env, only: output_unit

    implicit none

    ! internal subroutines used by bstree
    private bstree_deallocate
    private bstree_destroy
    private bstree_insert
    private bstree_locate

    ! internal subroutines used by bstnode
    private bstnode_insert
    private bstnode_insert_perform
    private bstnode_insert_repair
    private bstnode_locate
    private bstnode_rotate_left
    private bstnode_rotate_right
    private bstnode_uncle

    !> \brief   Node of the binary search tree
    !> \authors J Benda
    !> \date    2018
    !>
    !> Every node carries pointer to the parent node (or null if a node is the root node),
    !> pointers to two child nodes, color (black/red) and the interesting payload that defines
    !> value and order of the nodes. In this implementation it is the index of a specific
    !> data item in the storage of data items.
    !>
    type bstnode
        type(bstnode), pointer :: parent => null()  !< parent node (or null if root)
        type(bstnode), pointer :: left   => null()  !< left child nodes; when not allocated, corresponds to a black leaf
        type(bstnode), pointer :: right  => null()  !< right child nodes; when not allocated, corresponds to a black leaf
        integer :: id                               !< payload, often index to an array with actual data
        integer :: color                            !< black(0)/red(1) indicator for this node
    contains
        procedure, pass :: insert  => bstnode_insert
        procedure, pass :: locate  => bstnode_locate
        procedure, pass :: output  => bstnode_output
        procedure, pass :: uncle   => bstnode_uncle
        final           :: bstnode_deallocate
    end type bstnode

    !> \brief Binary search tree wrapper used for fast search of data items
    !> \authors J Benda
    !> \date    2018
    !>
    !> This is a simple implementation of the standard red-black self-balancing binary search tree.
    !> It is an order-2 recursive data structure, each node of which is associated with one data item.
    !> Data items are represented by integers, typically pointers to some linear data storage.
    !> The crutial property of binary search trees is that their nodes are ordered in such a way that
    !> \verbatim
    !>    left < parent < right
    !> \endverbatim
    !> where \c left and \c right stand for the data items associated with the two child nodes L and R
    !> of any node P, \c parent stands for the data item associated with the node P, and the *less than*
    !> relation is defined when needed by an auxiliary comparator. This data structure allows
    !> addition and, most notably, look-up of nodes (i.e. data items) in O(log(N)) time, where N is the
    !> number of nodes (data items) already in the graph.
    !>
    !> The tree data structure relies on use of pointers and dynamic allocation of child nodes.
    !>
    type bstree
        type(bstnode), pointer :: root => null()    !< root node of the tree
    contains
        procedure, pass :: insert  => bstree_insert
        procedure, pass :: locate  => bstree_locate
        procedure, pass :: output  => bstree_output
        procedure, pass :: dtwrite => bstree_dtwrite
        procedure, pass :: compare => bstree_compare
        procedure, pass :: destroy => bstree_destroy
        final           :: bstree_deallocate
    end type bstree

contains

    !> \brief   Destructor for the binary search tree
    !> \authors J Benda
    !> \date    2018
    !>
    subroutine bstree_deallocate (this)

        type(bstree) :: this

        call this % destroy

    end subroutine bstree_deallocate


    !> \brief   Deallocate whole tree
    !> \authors J Benda
    !> \date    2018
    !>
    !> Deallocates the root node wrapped in the \c bstree object.
    !>
    subroutine bstree_destroy (this)

        class(bstree) :: this

        ! deallocation of the root node recursively deallocates all nodes
        ! via their own destructors
        if (associated(this % root)) then
            deallocate(this % root)
        end if

    end subroutine bstree_destroy


    !> \brief   Wrapper around bstnode::insert
    !> \authors J Benda
    !> \date    2018 - 2019
    !>
    !> Insert new element (data item) to the tree. This subroutine forwards the
    !> call to bstnode::insert of the root node. Insertion of element, followed by
    !> tree re-balancing may result in change of the root node. If this happens,
    !> this subroutine adequately adjusts the root pointer bstree::root.
    !>
    !> \param this  Tree where the element is to be inserted
    !> \param idx   Element to add
    !> \param data  Optional pointer to data to pass to the comparator
    !>
    subroutine bstree_insert (this, idx, data)

        use iso_c_binding, only: c_ptr

        class(bstree), intent(inout)        :: this
        integer,       intent(in)           :: idx
        type(c_ptr),   intent(in), optional :: data

        ! if the root node is not set, just use the new element as root and return
        if (.not. associated(this % root)) then
            allocate(this % root)
            this % root % parent => null()
            this % root % color = 0
            this % root % id = idx
        else
            this % root => this % root % insert(this, idx, data)
        end if

    end subroutine bstree_insert


    !> \brief   Wrapper around bstnode::locate
    !> \authors J Benda
    !> \date    2018 - 2019
    !>
    !> \param this  Tree where the data item is to be searched for
    !> \param val   Value to find
    !> \param data  Optional pointer to data to pass to the comparator
    !>
    !> \return Index of the found value or -1.
    !>
    integer function bstree_locate (this, val, data) result (idx)

        use iso_c_binding, only: c_ptr

        class(bstree), intent(in)           :: this
        integer,       intent(in)           :: val
        type(c_ptr),   intent(in), optional :: data

        if (associated(this % root)) then
            idx = this % root % locate(this, val, data)
        else
            idx = -1
        end if

    end function bstree_locate


    !> \brief   Write the tree (for debugging purposes).
    !> \authors J Benda
    !> \date    2018
    !>
    !> \param this    Binary search tree
    !> \param indent  White-space indentation of lines
    !>
    subroutine bstree_output (this, indent)

        class(bstree), intent(in) :: this
        integer,       intent(in) :: indent

        integer :: i

        if (associated(this % root)) then
            write(*,'(a)',advance='no') repeat(' ', indent)
            call this % root % output(this, indent)
        end if

    end subroutine bstree_output


    !> \brief   Write a single-line representation of the data item
    !> \authors J Benda
    !> \date    2019
    !>
    !> This is used in debugging output of the tree. Every data item will be given an opportunity
    !> to print itself to a (if possible) short single line. The default implementation of this
    !> subroutine writes an empty line, but derived types can override it according to their desire.
    !>
    !> \param this  Binary search tree
    !> \param lu    Unit for writing
    !> \param id    Tree node id
    !>
    subroutine bstree_dtwrite (this, lu, id)

        class(bstree), intent(in) :: this
        integer,       intent(in) :: lu, id

        write (lu, *)

    end subroutine bstree_dtwrite


    !> \brief   Binary predicate
    !> \authors J Benda
    !> \date    2019
    !>
    !> This is the relation (predicate) that defines ordering for the data items represented by node IDs.
    !> In the basic implementation this just compares node IDs using spaceship operator (i.e. -1 when i < j,
    !> 0 when i == j, and +1 when i > j). In derived types, this should be overriden to correctly describe
    !> ordering of non-trivial data items.
    !>
    !> \param this  Binary search tree
    !> \param i     First node ID
    !> \param j     Second node ID
    !> \param data  Optional pointer to data to pass to the comparator
    !>
    integer function bstree_compare (this, i, j, data) result (verdict)

        use iso_c_binding, only: c_ptr

        class(bstree), intent(in)           :: this
        integer,       intent(in)           :: i, j
        type(c_ptr),   intent(in), optional :: data

        verdict = max(-1, min(+1, i - j))

    end function bstree_compare


    !> \brief   Destructor for the binary search tree node
    !> \authors J Benda
    !> \date    2018
    !>
    !> Deallocates the leaves of the node.
    !>
    recursive subroutine bstnode_deallocate (this)

        type(bstnode) :: this

        if (associated(this % left))  deallocate(this % left)
        if (associated(this % right)) deallocate(this % right)

    end subroutine bstnode_deallocate


    !> \brief   Insert new node to the binary search tree
    !> \authors J Benda
    !> \date    2018 - 2019
    !>
    !> Inserts new node into the tree and re-balances the tree to satisfy the
    !> red-black tree constraints.
    !>
    !> \param this  Tree node where the element is to be inserted
    !> \param tree  Parent tree object
    !> \param idx   Index of the data item to add
    !> \param data  Optional pointer to data to pass to the comparator
    !>
    !> \return Pointer to the newly added node.
    !>
    function bstnode_insert (this, tree, idx, data) result (cur)

        use iso_c_binding, only: c_ptr

        class(bstnode)            :: this
        class(bstree), intent(in) :: tree
        integer,       intent(in) :: idx
        type(c_ptr),   intent(in), optional :: data

        type(bstnode), pointer :: cur

        ! insert node into already existing tree
        cur => bstnode_insert_perform(this, tree, idx, data)

        ! repair the tree
        call bstnode_insert_repair(cur)

        ! find and return the new root node
        do while (associated(cur % parent))
            cur => cur % parent
        end do

    end function bstnode_insert


    !> \brief   Insert new element to the search tree
    !> \authors J Benda
    !> \date    2018
    !>
    !> This searches the binary search tree for a proper place for the new element and
    !> adds it as a new node. If the element already exists somewhere in the tree, the
    !> tree is not modified by this function.
    !>
    !> The insertion is implemented in a non-recursive way (as a loop with moving pointer).
    !>
    !> \param bst   Tree node where the insertion will begin
    !> \param tree  Parent tree object
    !> \param idx   Index of the data item to add
    !> \param data  Additional optional data pointer to pass to comparator
    !>
    !> \return Pointer to the newly added node.
    !>
    function bstnode_insert_perform (bst, tree, idx, data) result (cur)

        use iso_c_binding, only: c_ptr

        class(bstnode), target, intent(inout) :: bst
        class(bstree),          intent(in)    :: tree
        integer,                intent(in)    :: idx
        type(c_ptr), optional,  intent(in)    :: data

        type(bstnode), pointer :: cur

        cur => bst                                      ! start the addition at the given node

        do while (associated(cur))

            select case (tree % compare(idx, cur % id, data))

                case (-1)                               ! supplied data item is smaller than this node's data item

                    if (associated(cur % left)) then    ! insert new node into the left subtree (if available) ...
                        cur => cur % left
                    else                                ! ...or as a new (red) left leaf
                        allocate (cur % left)
                        cur % left % parent => cur
                        cur % left % color  =  1
                        cur % left % id     =  idx
                        cur => cur % left
                        return
                    end if

                case (0)                                ! supplied data item is equal to this node's data item

                    return                              ! node already exists -> exit without doing enything

                case (1)                                ! supplied data item is larger than this node's data item

                    if (associated(cur % right)) then   ! insert new node into the right subtree (if available)...
                        cur => cur % right
                    else                                ! ...or as a new (red) right leaf
                        allocate (cur % right)
                        cur % right % parent => cur
                        cur % right % color  =  1
                        cur % right % id     =  idx
                        cur => cur % right
                        return
                    end if

            end select

        end do

    end function bstnode_insert_perform


    !> \brief   Repair the tree after insertion
    !> \authors J Benda
    !> \date    2018
    !>
    !> Insertion of new (= red) nodes into the tree destroys the defining properties of the tree. Not completely,
    !> the tree is still a valid binary search tree, so all look-ups will work correctly, but some other
    !> properties related to the red-black colouring may be violated, which then leads to performance
    !> degradation (i.e. an unbalanced tree).
    !>
    !> This subroutine will straighten it up again, to preserve the asymptotic O(log(N)) property of all operations.
    !>
    !> The rules that need to be fixed are the following:
    !>   1. The root node has to be black.
    !>   2. If a node is red, then both its children are black.
    !>   3. Every path from a given node to any of its descendant leaves contains the same number of black nodes.
    !>
    !> \param bst   Tree node where the repairs will begin
    !>
    subroutine bstnode_insert_repair (bst)

        class(bstnode), target :: bst

        type(bstnode), pointer :: n, p, g, u, glr, grl

        n => bst

        do while (associated(n))

            ! If given a (red) root node, paint it black and return.

            if (.not. associated(n % parent)) then
                n % color = 0
                return
            end if

            ! If given a (red) non-root node with a black parent, return.

            if (n % parent % color == 0) then
                return
            end if

            ! If given a (red) non-root node with red parent and red uncle,
            ! paint both parent and uncle black, grandparent red, and continue to
            ! fix the grandparent.

            p => n % parent
            u => n % uncle()
            g => p % parent

            if (associated(u)) then
                if (u % color == 1) then
                    p % color = 0
                    u % color = 0
                    g % color = 1
                    n => g
                    cycle
                end if
            end if

            ! Being here means that the current (red) node is non-root, its parent
            ! is red (this is what bothers us) and its uncle is black. This is more
            ! involved to get right, because the nodes need to be reconnected in a
            ! very special way, called "tree rotation". See the literature on red-black
            ! trees for details.

            ! First, if the current node is (horizontally) between its parent and uncle,
            ! rotate the parent-node edge  so that the node gets to the other side of its parent
            ! (left or right) than the uncle is (right or left).

            glr => null() ; if (associated(g % left))  glr => g % left % right
            grl => null() ; if (associated(g % right)) grl => g % right % left

            if (associated(n, glr)) then
                call bstnode_rotate_left (p)
                n => n % left
            else if (associated(n, grl)) then
                call bstnode_rotate_right (p)
                n => n % right
            end if

            ! Second, promote parent to grandparent, making the former grandparent
            ! a sibling of the current node, and making the former uncle (if any) a nephew of
            ! the current node.

            p => n % parent
            g => p % parent

            if (associated(n, p % left)) then
                call bstnode_rotate_right (g)
            else
                call bstnode_rotate_left (g)
            end if

            ! Finally, switch the colors of parent and grandparent.
            p % color = 0
            g % color = 1

            return

        end do

    end subroutine bstnode_insert_repair


    !> \brief   Find the supplied data item in the storage using the binary tree
    !> \authors J Benda
    !> \date    2018
    !>
    !> Search the binary tree for a data item that evaluates the predicate to zero.
    !> The search starts at the root node. At every node, the given data item
    !> is compared to the node's data item. If equal, the node's index is returned.
    !> If less, the search continues in the left sub-tree. Otherwise the search continues
    !> in the right sub-tree.
    !>
    !> \param this  Tree node where the data item is to be searched for
    !> \param tree  Parent tree object.
    !> \param val   Value to search for; has no meaning when predicate is present
    !> \param data  Optional pointer to data to pass to the comparator
    !>
    !> \return Index of the data item in the storage or -1 if not found.
    !>
    integer function bstnode_locate (this, tree, val, data) result (idx)

        use iso_c_binding, only: c_ptr

        class(bstnode), target, intent(in) :: this
        class(bstree),          intent(in) :: tree
        integer,                intent(in) :: val
        type(c_ptr),  optional, intent(in) :: data

        type(bstnode), pointer :: cur

        cur => this

        do while (associated(cur))

            select case (tree % compare(val, cur % id, data))
                case (-1) ; cur => cur % left
                case ( 0) ; idx = cur % id ; return
                case ( 1) ; cur => cur % right
            end select

        end do

        idx = -1  ! not found

    end function bstnode_locate


    !> \brief   Write the sub-tree (for debugging purposes).
    !> \authors J Benda
    !> \date    2018
    !>
    !> \param this    Node of the binary search tree
    !> \param tree    Parent tree object
    !> \param indent  White-space indentation of lines
    !>
    recursive subroutine bstnode_output (this, tree, indent)

        class(bstnode), intent(in) :: this
        class(bstree),  intent(in) :: tree
        integer,        intent(in) :: indent

        integer :: i, pid

        if (associated(this % parent)) then
            pid = this % parent % id
        else
            pid = 0
        end if

        if (this % color == 0) then
            write(*,'(I0,"(black, child of ",I0,")")', advance = 'no') this % id, pid
        else
            write(*,'(I0,"(red, child of ",I0,")")', advance = 'no') this % id, pid
        end if

        call tree % dtwrite(output_unit, this % id)

        if (associated(this % left)) then
            write(*,'(2a)',advance='no') repeat(' ', indent), 'L '
            call this % left % output(tree, indent + 2)
        end if

        if (associated(this % right)) then
            write(*,'(2a)',advance='no') repeat(' ', indent), 'R '
            call this % right % output(tree, indent + 2)
        end if

    end subroutine bstnode_output


    !> \brief Left tree rotation
    !>
    !> Left tree rotation is used in re-balancing of the binary search tree.
    !> On entry, the given node has a right child. On entry, the relation
    !> between the node and its child will be changed, making the node a left
    !> child of its former child (it becoming the parent).
    !>
    !> The operation is schematically shown below, which illustrates that it can
    !> be carried out without destroying the horizontal order (here "1,P,3,N,5,C,7")
    !> that is relied upon during searches and insertions. The parent node P
    !> in this diagram is optinal (N can be originally the root node). The nubered
    !> nodes can be either leaves or arbitrary subtrees.
    !>
    !> \verbatim
    !>     (P)                     (P)
    !>    /   \                   /   \
    !>   1     N                 1     C
    !>        / \         =>          / \
    !>       3   C                   N   7
    !>          / \                 / \
    !>         5   7               3   5
    !> \endverbatim
    !>
    !> This is called "left" rotation, because the edge N-C rotates left.
    !>
    subroutine bstnode_rotate_left (N)

        type(bstnode), pointer :: N, P, C

        C => N % right
        P => N % parent

        ! move left leaf of the child to the vacated leaf of this node
        N % right => C % left
        if (associated(N % right)) N % right % parent => N

        ! move the node to the left leaf of the child (it becoming parent)
        C % left => N

        ! fix the family relations
        if (associated(P)) then
            if (associated(P % left,  N)) P % left  => C
            if (associated(P % right, N)) P % right => C
        end if
        C % parent => P
        N % parent => C

    end subroutine bstnode_rotate_left


    !> \brief Right tree rotation
    !>
    !> Right tree rotation is used in re-balancing of the binary search tree.
    !> On entry, the given node has a left child. On entry, the relation
    !> between the node and its child will be changed, making the node a right
    !> child of its former child (it becoming the parent).
    !>
    !> The operation is schematically shown below, which illustrates that it can
    !> be carried out without destroying the horizontal order (here "1,C,3,N,5,p,7")
    !> that is relied upon during searches and insertions. The parent node P
    !> in this diagram is optinal (N can be originally the root node). The nubered
    !> nodes can be either leaves or arbitrary subtrees.
    !>
    !> \verbatim
    !>         (P)                     (P)
    !>        /   \                   /   \
    !>       N     7                 C     7
    !>      / \           =>        / \
    !>     C   5                   1   N
    !>    / \                         / \
    !>   1   3                       3   5
    !> \endverbatim
    !>
    !> This is called "right" rotation, because the edge C-N rotates right.
    !>
    subroutine bstnode_rotate_right (N)

        type(bstnode), pointer :: N, P, C

        C => N % left
        P => N % parent

        ! move right leaf of the child to the vacated leaf of this node
        N % left => C % right
        if (associated(N % left)) N % left % parent => N

        ! move the node to the right leaf of the child (it becoming parent)
        C % right => N

        ! fix the family relations
        if (associated(P)) then
            if (associated(P % left,  N)) P % left  => C
            if (associated(P % right, N)) P % right => C
        end if
        C % parent => P
        N % parent => C

    end subroutine bstnode_rotate_right


    !> \brief Get uncle node
    !>
    !> Returns pointer to the uncle node (sibling of parent) of the supplied node.
    !> If such node does not exist, it halts the program.
    !>
    function bstnode_uncle (this) result (cur)

        class(bstnode), intent(in), target :: this
        type(bstnode), pointer :: cur

        cur => this

        if (.not. associated(cur % parent)) then
            ! no parent -> no uncle
            stop 1
        end if

        if (.not. associated(cur % parent % parent)) then
            ! no granparent -> no uncle
            stop 1
        end if

        if (associated(cur % parent, cur % parent % parent % left)) then
            cur => cur % parent % parent % right
            return
        end if

        if (associated(cur % parent, cur % parent % parent % right)) then
            cur => cur % parent % parent % left
            return
        end if

        ! broken tree
        stop 1

    end function bstnode_uncle

end module containers
