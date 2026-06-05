!
! Copyright (c) 2000-2008, Roland Schmehl. All rights reserved.
!
! This software is distributable under the BSD license. See the terms of the
! BSD license in the documentation provided with this software.
!
module fparser2
   !------- -------- --------- --------- --------- --------- --------- --------- -------
   ! Fortran 90 function parser v1.1
   !------- -------- --------- --------- --------- --------- --------- --------- -------
   !
   ! This function parser module is intended for applications where a set of mathematical
   ! fortran-style expressions is specified at runtime and is then evaluated for a large
   ! number of variable values. This is done by compiling the set of function strings
   ! into byte code, which is interpreted efficiently for the various variable values.
   !
   ! The source code is available from http://fparser.sourceforge.net
   !
   ! Please send comments, corrections or questions to the author:
   ! Roland Schmehl <roland.schmehl@alumni.uni-karlsruhe.de>
   !
   !------- -------- --------- --------- --------- --------- --------- --------- -------
   ! The function parser concept is based on a C++ class library written by  Juha
   ! Nieminen <warp@iki.fi> available from http://warp.povusers.org/FunctionParser/
   !------- -------- --------- --------- --------- --------- --------- --------- -------
   use SIZES ! Import KIND parameters
   implicit none
   !------- -------- --------- --------- --------- --------- --------- --------- -------
   public                     :: initf2, & ! Initialize function parser for n functions
                                 parsef2, & ! Parse single function string
                                 evalf2, & ! Evaluate single function
                                 EvalErrMsg2 ! Error message (Use only when EvalErrType>0)
   integer, public            :: EvalErrType2 ! =0: no error occured, >0: evaluation error
   !------- -------- --------- --------- --------- --------- --------- --------- -------
   private
   save
   integer(sz), parameter :: cImmed = 1, &
                             cNeg = 2, &
                             cAdd = 3, &
                             cSub = 4, &
                             cMul = 5, &
                             cDiv = 6, &
                             cPow = 7, &
                             cAbs = 8, &
                             cExp = 9, &
                             cLog10 = 10, &
                             cLog = 11, &
                             cSqrt = 12, &
                             cSinh = 13, &
                             cCosh = 14, &
                             cTanh = 15, &
                             cSin = 16, &
                             cCos = 17, &
                             cTan = 18, &
                             cAsin = 19, &
                             cAcos = 20, &
                             cAtan = 21, &
                             VarBegin = 22
   character(LEN=1), dimension(cAdd:cPow), parameter :: Ops = ['+', &
                                                               '-', &
                                                               '*', &
                                                               '/', &
                                                               '^']
   character(LEN=5), dimension(cAbs:cAtan), parameter :: Funcs = ['abs  ', &
                                                                  'exp  ', &
                                                                  'log10', &
                                                                  'log  ', &
                                                                  'sqrt ', &
                                                                  'sinh ', &
                                                                  'cosh ', &
                                                                  'tanh ', &
                                                                  'sin  ', &
                                                                  'cos  ', &
                                                                  'tan  ', &
                                                                  'asin ', &
                                                                  'acos ', &
                                                                  'atan ']
   type tComp
      integer(sz), dimension(:), pointer :: ByteCode
      integer                            :: ByteCodeSize
      real(sz), dimension(:), pointer :: Immed
      integer                            :: ImmedSize
      real(sz), dimension(:), pointer :: Stack
      integer                            :: StackSize, &
                                            StackPtr
   end type tComp
   type(tComp), dimension(:), pointer :: Comp ! Bytecode
   integer, dimension(:), allocatable :: ipos ! Associates function strings
   !
contains
   !
   subroutine initf2(n)
      !----- -------- --------- --------- --------- --------- --------- --------- -------
      ! Initialize function parser for n functions
      !----- -------- --------- --------- --------- --------- --------- --------- -------
      implicit none
      integer, intent(in) :: n ! Number of functions
      integer             :: i
      !----- -------- --------- --------- --------- --------- --------- --------- -------
      allocate (Comp(n))
      do i = 1, n
         nullify (Comp(i)%ByteCode, Comp(i)%Immed, Comp(i)%Stack)
      end do
   end subroutine initf2
   !
   subroutine parsef2(i, FuncStr, Var)
      !----- -------- --------- --------- --------- --------- --------- --------- -------
      ! Parse ith function string FuncStr and compile it into bytecode
      !----- -------- --------- --------- --------- --------- --------- --------- -------
      implicit none
      integer, intent(in) :: i ! Function identifier
      character(LEN=*), intent(in) :: FuncStr ! Function string
      character(LEN=*), dimension(:), intent(in) :: Var ! Array with variable names
      character(LEN=len(FuncStr))                :: Func ! Function string, local use
      !----- -------- --------- --------- --------- --------- --------- --------- -------
      if (i < 1 .or. i > size(Comp)) then
         write (*, *) '*** Parser error: Function number ', i, ' out of range'
         stop
      end if
      allocate (ipos(len_trim(FuncStr))) ! Char. positions in orig. string
      Func = FuncStr ! Local copy of function string
      call Replace('**', '^ ', Func) ! Exponent into 1-Char. format
      call RemoveSpaces(Func) ! Condense function string
      !print*,Func,FuncStr,Var
      call CheckSyntax(Func, FuncStr, Var)
      deallocate (ipos)
      call Compile(i, Func, Var) ! Compile into bytecode
   end subroutine parsef2
   !
   function evalf2(i, Val) result(res)
      !----- -------- --------- --------- --------- --------- --------- --------- -------
      ! Evaluate bytecode of ith function for the values passed in array Val(:)
      !----- -------- --------- --------- --------- --------- --------- --------- -------
      implicit none
      integer, intent(in) :: i ! Function identifier
      real(sz), dimension(:), intent(in) :: Val ! Variable values
      real(sz)                           :: res ! Result
      integer, target                     :: IP, & ! Instruction pointer
                                             DP, & ! Data pointer
                                             SP ! Stack pointer
      real(sz), parameter :: zero = 0._sz
      !----- -------- --------- --------- --------- --------- --------- --------- -------
      DP = 1
      SP = 0
      do IP = 1, Comp(i)%ByteCodeSize
         select case (Comp(i)%ByteCode(IP))

         case (cImmed); SP = SP + 1; Comp(i)%Stack(SP) = Comp(i)%Immed(DP); DP = DP + 1
         case (cNeg); Comp(i)%Stack(SP) = -Comp(i)%Stack(SP)
         case (cAdd); Comp(i)%Stack(SP - 1) = Comp(i)%Stack(SP - 1) + Comp(i)%Stack(SP); SP = SP - 1
         case (cSub); Comp(i)%Stack(SP - 1) = Comp(i)%Stack(SP - 1) - Comp(i)%Stack(SP); SP = SP - 1
         case (cMul); Comp(i)%Stack(SP - 1) = Comp(i)%Stack(SP - 1)*Comp(i)%Stack(SP); SP = SP - 1
         case (cDiv); if (Comp(i)%Stack(SP) == 0._sz) then; EvalErrType2 = 1; res = zero; return; end if
            Comp(i)%Stack(SP - 1) = Comp(i)%Stack(SP - 1)/Comp(i)%Stack(SP); SP = SP - 1
         case (cPow); Comp(i)%Stack(SP - 1) = Comp(i)%Stack(SP - 1)**Comp(i)%Stack(SP); SP = SP - 1
         case (cAbs); Comp(i)%Stack(SP) = abs(Comp(i)%Stack(SP))
         case (cExp); Comp(i)%Stack(SP) = exp(Comp(i)%Stack(SP))
         case (cLog10); if (Comp(i)%Stack(SP) <= 0._sz) then; EvalErrType2 = 3; res = zero; return; end if
            Comp(i)%Stack(SP) = log10(Comp(i)%Stack(SP))
         case (cLog); if (Comp(i)%Stack(SP) <= 0._sz) then; EvalErrType2 = 3; res = zero; return; end if
            Comp(i)%Stack(SP) = log(Comp(i)%Stack(SP))
         case (cSqrt); if (Comp(i)%Stack(SP) < 0._sz) then; EvalErrType2 = 3; res = zero; return; end if
            Comp(i)%Stack(SP) = sqrt(Comp(i)%Stack(SP))
         case (cSinh); Comp(i)%Stack(SP) = sinh(Comp(i)%Stack(SP))
         case (cCosh); Comp(i)%Stack(SP) = cosh(Comp(i)%Stack(SP))
         case (cTanh); Comp(i)%Stack(SP) = tanh(Comp(i)%Stack(SP))
         case (cSin); Comp(i)%Stack(SP) = sin(Comp(i)%Stack(SP))
         case (cCos); Comp(i)%Stack(SP) = cos(Comp(i)%Stack(SP))
         case (cTan); Comp(i)%Stack(SP) = tan(Comp(i)%Stack(SP))
         case (cAsin); if ((Comp(i)%Stack(SP) < -1._sz) .or. (Comp(i)%Stack(SP) > 1._sz)) then
               EvalErrType2 = 4; res = zero; return; end if
            Comp(i)%Stack(SP) = asin(Comp(i)%Stack(SP))
         case (cAcos); if ((Comp(i)%Stack(SP) < -1._sz) .or. (Comp(i)%Stack(SP) > 1._sz)) then
               EvalErrType2 = 4; res = zero; return; end if
            Comp(i)%Stack(SP) = acos(Comp(i)%Stack(SP))
         case (cAtan); Comp(i)%Stack(SP) = atan(Comp(i)%Stack(SP))
         case DEFAULT; SP = SP + 1; Comp(i)%Stack(SP) = Val(Comp(i)%ByteCode(IP) - VarBegin + 1)
         end select
      end do
      EvalErrType2 = 0
      res = Comp(i)%Stack(1)
   end function evalf2
   !
   subroutine CheckSyntax(Func, FuncStr, Var)
      !----- -------- --------- --------- --------- --------- --------- --------- -------
      ! Check syntax of function string,  returns 0 if syntax sz ok
      !----- -------- --------- --------- --------- --------- --------- --------- -------
      implicit none
      character(LEN=*), intent(in) :: Func ! Function string without spaces
      character(LEN=*), intent(in) :: FuncStr ! Original function string
      character(LEN=*), dimension(:), intent(in) :: Var ! Array with variable names
      integer(sz)                                 :: n
      character(LEN=1)                           :: c
      real(sz)                                    :: r
      logical                                     :: err
      integer                                     :: ParCnt, & ! Parenthesis counter
                                                     j, ib, in, lFunc
      !----- -------- --------- --------- --------- --------- --------- --------- -------
      j = 1
      ParCnt = 0
      lFunc = len_trim(Func)
      !print*,lFunc
      step: do
         if (j > lFunc) call ParseErrMsg(j, FuncStr)
         c = Func(j:j)
         !-- -------- --------- --------- --------- --------- --------- --------- -------
         ! Check for valid operand (must appear)
         !-- -------- --------- --------- --------- --------- --------- --------- -------
         if (c == '-' .or. c == '+') then ! Check for leading - or +
            j = j + 1
            if (j > lFunc) call ParseErrMsg(j, FuncStr, 'Missing operand')
            c = Func(j:j)
            if (any(c == Ops)) call ParseErrMsg(j, FuncStr, 'Multiple operators')
         end if
         n = MathFunctionIndex(Func(j:))
         if (n > 0) then ! Check for math function
            j = j + len_trim(Funcs(n))
            if (j > lFunc) call ParseErrMsg(j, FuncStr, 'Missing function argument')
            c = Func(j:j)
            if (c /= '(') call ParseErrMsg(j, FuncStr, 'Missing opening parenthesis')
         end if
         if (c == '(') then ! Check for opening parenthesis
            ParCnt = ParCnt + 1
            j = j + 1
            cycle step
         end if
         if (scan(c, '0123456789.') > 0) then ! Check for number
            r = RealNum(Func(j:), ib, in, err)
            if (err) call ParseErrMsg(j, FuncStr, 'Invalid number format:  '//Func(j + ib - 1:j + in - 2))
            j = j + in - 1
            if (j > lFunc) exit
            c = Func(j:j)
         else ! Check for variable
            n = VariableIndex(Func(j:), Var, ib, in)
            !print*,Func(J:),Var,ib,in,'wtf?',n
            if (n == 0) call ParseErrMsg(j, FuncStr, 'Invalid element: '//Func(j + ib - 1:j + in - 2))
            j = j + in - 1
            if (j > lFunc) exit
            c = Func(j:j)
         end if
         do while (c == ')') ! Check for closing parenthesis
            ParCnt = ParCnt - 1
            if (ParCnt < 0) call ParseErrMsg(j, FuncStr, 'Mismatched parenthesis')
            if (Func(j - 1:j - 1) == '(') call ParseErrMsg(j - 1, FuncStr, 'Empty parentheses')
            j = j + 1
            if (j > lFunc) exit
            c = Func(j:j)
         end do
         !-- -------- --------- --------- --------- --------- --------- --------- -------
         ! Now, we have a legal operand: A legal operator or end of string must follow
         !-- -------- --------- --------- --------- --------- --------- --------- -------
         if (j > lFunc) exit
         if (any(c == Ops)) then ! Check for multiple operators
            if (j + 1 > lFunc) call ParseErrMsg(j, FuncStr)
            if (any(Func(j + 1:j + 1) == Ops)) call ParseErrMsg(j + 1, FuncStr, 'Multiple operators')
         else ! Check for next operand
            call ParseErrMsg(j, FuncStr, 'Missing operator')
         end if
         !-- -------- --------- --------- --------- --------- --------- --------- -------
         ! Now, we have an operand and an operator: the next loop will check for another
         ! operand (must appear)
         !-- -------- --------- --------- --------- --------- --------- --------- -------
         j = j + 1
      end do step
      if (ParCnt > 0) call ParseErrMsg(j, FuncStr, 'Missing )')
   end subroutine CheckSyntax
   !
   function EvalErrMsg2() result(msg)
      !----- -------- --------- --------- --------- --------- --------- --------- -------
      ! Return error message
      !----- -------- --------- --------- --------- --------- --------- --------- -------
      implicit none
      character(LEN=*), dimension(4), parameter :: m = ['Division by zero                ', &
                                                        'Argument of SQRT negative       ', &
                                                        'Argument of LOG negative        ', &
                                                        'Argument of ASIN or ACOS illegal']
      character(LEN=len(m))                     :: msg
      !----- -------- --------- --------- --------- --------- --------- --------- -------
      if (EvalErrType2 < 1 .or. EvalErrType2 > size(m)) then
         msg = ''
      else
         msg = m(EvalErrType2)
      end if
   end function EvalErrMsg2
   !
   subroutine ParseErrMsg(j, FuncStr, Msg)
      !----- -------- --------- --------- --------- --------- --------- --------- -------
      ! Print error message and terminate program
      !----- -------- --------- --------- --------- --------- --------- --------- -------
      implicit none
      integer, intent(in) :: j
      character(LEN=*), intent(in) :: FuncStr ! Original function string
      character(LEN=*), optional, intent(in) :: Msg
      integer                                 :: k
      !----- -------- --------- --------- --------- --------- --------- --------- -------
      if (present(Msg)) then
         write (*, *) '*** Error in syntax of function string: '//Msg
      else
         write (*, *) '*** Error in syntax of function string:'
      end if
      write (*, *)
      write (*, '(A)') ' '//FuncStr
      do k = 1, ipos(j)
         write (*, '(A)', ADVANCE='NO') ' ' ! Advance to the jth position
      end do
      write (*, '(A)') '?'
      stop
   end subroutine ParseErrMsg
   !
   function OperatorIndex(c) result(n)
      !----- -------- --------- --------- --------- --------- --------- --------- -------
      ! Return operator index
      !----- -------- --------- --------- --------- --------- --------- --------- -------
      implicit none
      character(LEN=1), intent(in) :: c
      integer(sz)                   :: n, j
      !----- -------- --------- --------- --------- --------- --------- --------- -------
      n = 0
      do j = cAdd, cPow
         if (c == Ops(j)) then
            n = j
            exit
         end if
      end do
   end function OperatorIndex
   !
   function MathFunctionIndex(str) result(n)
      !----- -------- --------- --------- --------- --------- --------- --------- -------
      ! Return index of math function beginnig at 1st position of string str
      !----- -------- --------- --------- --------- --------- --------- --------- -------
      implicit none
      character(LEN=*), intent(in) :: str
      integer(sz)                   :: n, j
      integer                       :: k
      character(LEN=len(Funcs))    :: fun
      !----- -------- --------- --------- --------- --------- --------- --------- -------
      n = 0
      do j = cAbs, cAtan ! Check all math functions
         k = min(len_trim(Funcs(j)), len(str))
         call LowCase(str(1:k), fun)
         if (fun == Funcs(j)) then ! Compare lower case letters
            n = j ! Found a matching function
            exit
         end if
      end do
   end function MathFunctionIndex
   !
   function VariableIndex(str, Var, ibegin, inext) result(n)
      !----- -------- --------- --------- --------- --------- --------- --------- -------
      ! Return index of variable at begin of string str (returns 0 if no variable found)
      !----- -------- --------- --------- --------- --------- --------- --------- -------
      implicit none
      character(LEN=*), intent(in) :: str ! String
      character(LEN=*), dimension(:), intent(in) :: Var ! Array with variable names
      integer(sz)                                 :: n ! Index of variable
      integer, optional, intent(out) :: ibegin, & ! Start position of variable name
                                        inext ! Position of character after name
      integer                                     :: j, ib, in, lstr
      !----- -------- --------- --------- --------- --------- --------- --------- -------
      n = 0
      lstr = len_trim(str)
      if (lstr > 0) then
         do ib = 1, lstr ! Search for first character in str
            if (str(ib:ib) /= ' ') exit ! When lstr>0 at least 1 char in str
         end do
         do in = ib, lstr ! Search for name terminators
            if (scan(str(in:in), '+-*/^) ') > 0) exit
         end do
         do j = 1, size(Var)
            if (str(ib:in - 1) == Var(j)) then
               n = j ! Variable name found
               exit
            end if
         end do
      end if
      if (present(ibegin)) ibegin = ib
      if (present(inext)) inext = in
   end function VariableIndex
   !
   subroutine RemoveSpaces(str)
      !----- -------- --------- --------- --------- --------- --------- --------- -------
      ! Remove Spaces from string, remember positions of characters in old string
      !----- -------- --------- --------- --------- --------- --------- --------- -------
      implicit none
      character(LEN=*), intent(inout) :: str
      integer                          :: k, lstr
      !----- -------- --------- --------- --------- --------- --------- --------- -------
      lstr = len_trim(str)
      ipos = [(k, k=1, lstr)]
      k = 1
      do while (str(k:lstr) /= ' ')
         if (str(k:k) == ' ') then
            str(k:lstr) = str(k + 1:lstr)//' ' ! Move 1 character to left
            ipos(k:lstr) = [ipos(k + 1:lstr), 0] ! Move 1 element to left
            k = k - 1
         end if
         k = k + 1
      end do
   end subroutine RemoveSpaces
   !
   subroutine Replace(ca, cb, str)
      !----- -------- --------- --------- --------- --------- --------- --------- -------
      ! Replace ALL appearances of character set ca in string str by character set cb
      !----- -------- --------- --------- --------- --------- --------- --------- -------
      implicit none
      character(LEN=*), intent(in) :: ca
      character(LEN=len(ca)), intent(in) :: cb ! LEN(ca) must be LEN(cb)
      character(LEN=*), intent(inout) :: str
      integer                             :: j, lca
      !----- -------- --------- --------- --------- --------- --------- --------- -------
      lca = len(ca)
      do j = 1, len_trim(str) - lca + 1
         if (str(j:j + lca - 1) == ca) str(j:j + lca - 1) = cb
      end do
   end subroutine Replace
   !
   subroutine Compile(i, F, Var)
      !----- -------- --------- --------- --------- --------- --------- --------- -------
      ! Compile i-th function string F into bytecode
      !----- -------- --------- --------- --------- --------- --------- --------- -------
      implicit none
      integer, intent(in) :: i ! Function identifier
      character(LEN=*), intent(in) :: F ! Function string
      character(LEN=*), dimension(:), intent(in) :: Var ! Array with variable names
      integer                                     :: istat
      !----- -------- --------- --------- --------- --------- --------- --------- -------
      if (associated(Comp(i)%ByteCode)) deallocate (Comp(i)%ByteCode, &
                                                    Comp(i)%Immed, &
                                                    Comp(i)%Stack)
      Comp(i)%ByteCodeSize = 0
      Comp(i)%ImmedSize = 0
      Comp(i)%StackSize = 0
      Comp(i)%StackPtr = 0
      call CompileSubstr(i, F, 1, len_trim(F), Var) ! Compile string to determine size
      allocate (Comp(i)%ByteCode(Comp(i)%ByteCodeSize), &
                Comp(i)%Immed(Comp(i)%ImmedSize), &
                Comp(i)%Stack(Comp(i)%StackSize), &
                STAT=istat)
      if (istat /= 0) then
         write (*, *) '*** Parser error: Memmory allocation for byte code failed'
         stop
      else
         Comp(i)%ByteCodeSize = 0
         Comp(i)%ImmedSize = 0
         Comp(i)%StackSize = 0
         Comp(i)%StackPtr = 0
         call CompileSubstr(i, F, 1, len_trim(F), Var) ! Compile string into bytecode
      end if
      !
   end subroutine Compile
   !
   subroutine AddCompiledByte(i, b)
      !----- -------- --------- --------- --------- --------- --------- --------- -------
      ! Add compiled byte to bytecode
      !----- -------- --------- --------- --------- --------- --------- --------- -------
      implicit none
      integer, intent(in) :: i ! Function identifier
      integer(sz), intent(in) :: b ! Value of byte to be added
      !----- -------- --------- --------- --------- --------- --------- --------- -------
      Comp(i)%ByteCodeSize = Comp(i)%ByteCodeSize + 1
      if (associated(Comp(i)%ByteCode)) Comp(i)%ByteCode(Comp(i)%ByteCodeSize) = b
   end subroutine AddCompiledByte
   !
   function MathItemIndex(i, F, Var) result(n)
      !----- -------- --------- --------- --------- --------- --------- --------- -------
      ! Return math item index, if item is real number, enter it into Comp-structure
      !----- -------- --------- --------- --------- --------- --------- --------- -------
      implicit none
      integer, intent(in) :: i ! Function identifier
      character(LEN=*), intent(in) :: F ! Function substring
      character(LEN=*), dimension(:), intent(in) :: Var ! Array with variable names
      integer(sz)                                 :: n ! Byte value of math item
      !----- -------- --------- --------- --------- --------- --------- --------- -------
      n = 0
      if (scan(F(1:1), '0123456789.') > 0) then ! Check for begin of a number
         Comp(i)%ImmedSize = Comp(i)%ImmedSize + 1
         if (associated(Comp(i)%Immed)) Comp(i)%Immed(Comp(i)%ImmedSize) = RealNum(F)
         n = cImmed
      else ! Check for a variable
         n = VariableIndex(F, Var)
         if (n > 0) n = VarBegin + n - 1
      end if
   end function MathItemIndex
   !
   function CompletelyEnclosed(F, b, e) result(res)
      !----- -------- --------- --------- --------- --------- --------- --------- -------
      ! Check if function substring F(b:e) is completely enclosed by a pair of parenthesis
      !----- -------- --------- --------- --------- --------- --------- --------- -------
      implicit none
      character(LEN=*), intent(in) :: F ! Function substring
      integer, intent(in) :: b, e ! First and last pos. of substring
      logical                       :: res
      integer                       :: j, k
      !----- -------- --------- --------- --------- --------- --------- --------- -------
      res = .false.
      if (F(b:b) == '(' .and. F(e:e) == ')') then
         k = 0
         do j = b + 1, e - 1
            if (F(j:j) == '(') then
               k = k + 1
            elseif (F(j:j) == ')') then
               k = k - 1
            end if
            if (k < 0) exit
         end do
         if (k == 0) res = .true. ! All opened parenthesis closed
      end if
   end function CompletelyEnclosed
   !
   recursive subroutine CompileSubstr(i, F, b, e, Var)
      !----- -------- --------- --------- --------- --------- --------- --------- -------
      ! Compile i-th function string F into bytecode
      !----- -------- --------- --------- --------- --------- --------- --------- -------
      implicit none
      integer, intent(in) :: i ! Function identifier
      character(LEN=*), intent(in) :: F ! Function substring
      integer, intent(in) :: b, e ! Begin and end position substring
      character(LEN=*), dimension(:), intent(in) :: Var ! Array with variable names
      integer(sz)                                 :: n
      integer                                     :: b2, j, k, io
      character(LEN=*), parameter :: calpha = 'abcdefghijklmnopqrstuvwxyz'// &
                                     'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
      !----- -------- --------- --------- --------- --------- --------- --------- -------
      ! Check for special cases of substring
      !----- -------- --------- --------- --------- --------- --------- --------- -------
      if (F(b:b) == '+') then ! Case 1: F(b:e) = '+...'
!      WRITE(*,*)'1. F(b:e) = "+..."'
         call CompileSubstr(i, F, b + 1, e, Var)
         return
      elseif (CompletelyEnclosed(F, b, e)) then ! Case 2: F(b:e) = '(...)'
!      WRITE(*,*)'2. F(b:e) = "(...)"'
         call CompileSubstr(i, F, b + 1, e - 1, Var)
         return
      elseif (scan(F(b:b), calpha) > 0) then
         n = MathFunctionIndex(F(b:e))
         if (n > 0) then
            b2 = b + index(F(b:e), '(') - 1
            if (CompletelyEnclosed(F, b2, e)) then ! Case 3: F(b:e) = 'fcn(...)'
!            WRITE(*,*)'3. F(b:e) = "fcn(...)"'
               call CompileSubstr(i, F, b2 + 1, e - 1, Var)
               call AddCompiledByte(i, n)
               return
            end if
         end if
      elseif (F(b:b) == '-') then
         if (CompletelyEnclosed(F, b + 1, e)) then ! Case 4: F(b:e) = '-(...)'
!         WRITE(*,*)'4. F(b:e) = "-(...)"'
            call CompileSubstr(i, F, b + 2, e - 1, Var)
            call AddCompiledByte(i, cNeg)
            return
         elseif (scan(F(b + 1:b + 1), calpha) > 0) then
            n = MathFunctionIndex(F(b + 1:e))
            if (n > 0) then
               b2 = b + index(F(b + 1:e), '(')
               if (CompletelyEnclosed(F, b2, e)) then ! Case 5: F(b:e) = '-fcn(...)'
!               WRITE(*,*)'5. F(b:e) = "-fcn(...)"'
                  call CompileSubstr(i, F, b2 + 1, e - 1, Var)
                  call AddCompiledByte(i, n)
                  call AddCompiledByte(i, cNeg)
                  return
               end if
            end if
         end if
      end if
      !----- -------- --------- --------- --------- --------- --------- --------- -------
      ! Check for operator in substring: check only base level (k=0), exclude expr. in ()
      !----- -------- --------- --------- --------- --------- --------- --------- -------
      do io = cAdd, cPow ! Increasing priority +-*/^
         k = 0
         do j = e, b, -1
            if (F(j:j) == ')') then
               k = k + 1
            elseif (F(j:j) == '(') then
               k = k - 1
            end if
            if (k == 0 .and. F(j:j) == Ops(io) .and. IsBinaryOp(j, F)) then
               if (any(F(j:j) == Ops(cMul:cPow)) .and. F(b:b) == '-') then ! Case 6: F(b:e) = '-...Op...' with Op > -
!               WRITE(*,*)'6. F(b:e) = "-...Op..." with Op > -'
                  call CompileSubstr(i, F, b + 1, e, Var)
                  call AddCompiledByte(i, cNeg)
                  return
               else ! Case 7: F(b:e) = '...BinOp...'
!               WRITE(*,*)'7. Binary operator',F(j:j)
                  call CompileSubstr(i, F, b, j - 1, Var)
                  call CompileSubstr(i, F, j + 1, e, Var)
                  call AddCompiledByte(i, OperatorIndex(Ops(io)))
                  Comp(i)%StackPtr = Comp(i)%StackPtr - 1
                  return
               end if
            end if
         end do
      end do
      !----- -------- --------- --------- --------- --------- --------- --------- -------
      ! Check for remaining items, i.e. variables or explicit numbers
      !----- -------- --------- --------- --------- --------- --------- --------- -------
      b2 = b
      if (F(b:b) == '-') b2 = b2 + 1
      n = MathItemIndex(i, F(b2:e), Var)
!   WRITE(*,*)'8. AddCompiledByte ',n
      call AddCompiledByte(i, n)
      Comp(i)%StackPtr = Comp(i)%StackPtr + 1
      if (Comp(i)%StackPtr > Comp(i)%StackSize) Comp(i)%StackSize = Comp(i)%StackSize + 1
      if (b2 > b) call AddCompiledByte(i, cNeg)
   end subroutine CompileSubstr
   !
   function IsBinaryOp(j, F) result(res)
      !----- -------- --------- --------- --------- --------- --------- --------- -------
      ! Check if operator F(j:j) in string F is binary operator
      ! Special cases already covered elsewhere:              (that is corrected in v1.1)
      ! - operator character F(j:j) is first character of string (j=1)
      !----- -------- --------- --------- --------- --------- --------- --------- -------
      implicit none
      integer, intent(in) :: j ! Position of Operator
      character(LEN=*), intent(in) :: F ! String
      logical                       :: res ! Result
      integer                       :: k
      logical                       :: Dflag, Pflag
      !----- -------- --------- --------- --------- --------- --------- --------- -------
      res = .true.
      if (F(j:j) == '+' .or. F(j:j) == '-') then ! Plus or minus sign:
         if (j == 1) then ! - leading unary operator ?
            res = .false.
         elseif (scan(F(j - 1:j - 1), '+-*/^(') > 0) then ! - other unary operator ?
            res = .false.
         elseif (scan(F(j + 1:j + 1), '0123456789') > 0 .and. & ! - in exponent of real number ?
                 scan(F(j - 1:j - 1), 'eEdD') > 0) then
            Dflag = .false.; Pflag = .false.
            k = j - 1
            do while (k > 1) !   step to the left in mantissa
               k = k - 1
               if (scan(F(k:k), '0123456789') > 0) then
                  Dflag = .true.
               elseif (F(k:k) == '.') then
                  if (Pflag) then
                     exit !   * EXIT: 2nd appearance of '.'
                  else
                     Pflag = .true. !   * mark 1st appearance of '.'
                  end if
               else
                  exit !   * all other characters
               end if
            end do
            if (Dflag .and. (k == 1 .or. scan(F(k:k), '+-*/^(') > 0)) res = .false.
         end if
      end if
   end function IsBinaryOp
   !
   function RealNum(str, ibegin, inext, error) result(res)
      !----- -------- --------- --------- --------- --------- --------- --------- -------
      ! Get real number from string - Format: [blanks][+|-][nnn][.nnn][e|E|d|D[+|-]nnn]
      !----- -------- --------- --------- --------- --------- --------- --------- -------
      implicit none
      character(LEN=*), intent(in) :: str ! String
      real(sz)                       :: res ! Real number
      integer, optional, intent(out) :: ibegin, & ! Start position of real number
                                        inext ! 1st character after real number
      logical, optional, intent(out) :: error ! Error flag
      integer                        :: ib, in, istat
      logical                        :: Bflag, & ! .T. at begin of number in str
                                        InMan, & ! .T. in mantissa of number
                                        Pflag, & ! .T. after 1st '.' encountered
                                        Eflag, & ! .T. at exponent identifier 'eEdD'
                                        InExp, & ! .T. in exponent of number
                                        DInMan, & ! .T. if at least 1 digit in mant.
                                        DInExp, & ! .T. if at least 1 digit in exp.
                                        err ! Local error flag
      !----- -------- --------- --------- --------- --------- --------- --------- -------
      Bflag = .true.; InMan = .false.; Pflag = .false.; Eflag = .false.; InExp = .false.
      DInMan = .false.; DInExp = .false.
      ib = 1
      in = 1
      do while (in <= len_trim(str))
         select case (str(in:in))
         case (' ') ! Only leading blanks permitted
            ib = ib + 1
            if (InMan .or. Eflag .or. InExp) exit
         case ('+', '-') ! Permitted only
            if (Bflag) then
               InMan = .true.; Bflag = .false. ! - at beginning of mantissa
            elseif (Eflag) then
               InExp = .true.; Eflag = .false. ! - at beginning of exponent
            else
               exit ! - otherwise STOP
            end if
         case ('0':'9') ! Mark
            if (Bflag) then
               InMan = .true.; Bflag = .false. ! - beginning of mantissa
            elseif (Eflag) then
               InExp = .true.; Eflag = .false. ! - beginning of exponent
            end if
            if (InMan) DInMan = .true. ! Mantissa contains digit
            if (InExp) DInExp = .true. ! Exponent contains digit
         case ('.')
            if (Bflag) then
               Pflag = .true. ! - mark 1st appearance of '.'
               InMan = .true.; Bflag = .false. !   mark beginning of mantissa
            elseif (InMan .and. .not. Pflag) then
               Pflag = .true. ! - mark 1st appearance of '.'
            else
               exit ! - otherwise STOP
            end if
         case ('e', 'E', 'd', 'D') ! Permitted only
            if (InMan) then
               Eflag = .true.; InMan = .false. ! - following mantissa
            else
               exit ! - otherwise STOP
            end if
         case DEFAULT
            exit ! STOP at all other characters
         end select
         in = in + 1
      end do
      err = (ib > in - 1) .or. (.not. DInMan) .or. ((Eflag .or. InExp) .and. .not. DInExp)
      if (err) then
         res = 0.0_sz
      else
         read (str(ib:in - 1), *, IOSTAT=istat) res
         err = istat /= 0
      end if
      if (present(ibegin)) ibegin = ib
      if (present(inext)) inext = in
      if (present(error)) error = err
   end function RealNum
   !
   subroutine LowCase(str1, str2)
      !----- -------- --------- --------- --------- --------- --------- --------- -------
      ! Transform upper case letters in str1 into lower case letters, result is str2
      !----- -------- --------- --------- --------- --------- --------- --------- -------
      implicit none
      character(LEN=*), intent(in) :: str1
      character(LEN=*), intent(out) :: str2
      integer                        :: j, k
      character(LEN=*), parameter :: lc = 'abcdefghijklmnopqrstuvwxyz'
      character(LEN=*), parameter :: uc = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
      !----- -------- --------- --------- --------- --------- --------- --------- -------
      str2 = str1
      do j = 1, len_trim(str1)
         k = index(uc, str1(j:j))
         if (k > 0) str2(j:j) = lc(k:k)
      end do
   end subroutine LowCase
   !
end module fparser2
