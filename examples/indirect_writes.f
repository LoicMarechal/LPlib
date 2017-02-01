
      program indirect_writes

      include 'lplib3.ins'
      include 'libmeshb7.ins'

      integer*4 MaxVer,MaxTet,NmbItr
      parameter (MaxVer=10000000)
      parameter (MaxTet=60000000)
      parameter (NmbItr=100)

      integer*4 i,j,k,NmbCpu,ver,dim,ref,TetTyp,VerTyp,ret
      integer*4 VerDeg(MaxVer),TetVer(4,MaxTet)
      integer*4 NmbPth,DepTab(4)
      integer*8 NmbVer,NmbTet,ParIdx,InpMsh
      real*8 VerTem(MaxVer),TetTem(MaxTet),t
      real*4 sta(2),acc

      external vertemp,tettemp

      InpMsh = gmfopenmesh("tet.meshb", GmfRead,ver,dim)
      if(InpMsh.eq.0) STOP ' cannot open file tet.meshb'

      NmbVer = gmfstatkwd(InpMsh, GmfVertices)
      if(NmbVer.gt.MaxVer) STOP ' too many vertices'

      NmbTet = gmfstatkwd(InpMsh, GmfTetrahedra)
      if(NmbTet.gt.MaxTet) STOP ' too many tetrahedra'

      print*, 'version  = ',ver,   'dimension  = ',dim
      print*, 'vertices = ',NmbVer,'tetrahedra = ',NmbTet

      ret = gmfgetblock(InpMsh,GmfTetrahedra,1,NmbVer,0_8,
     +                  GmfInt, TetVer(1,1), TetVer(1,NmbVer),
     +                  GmfInt, TetVer(2,1), TetVer(2,NmbVer),
     +                  GmfInt, TetVer(3,1), TetVer(3,NmbVer),
     +                  GmfInt, TetVer(4,1), TetVer(4,NmbVer),
     +                  GmfInt, ref, ref)

      ret = gmfclosemesh(InpMsh)

      ParIdx = initparallel(0)
      if(ParIdx.eq.0) STOP ' cannot init LPLIB3'
      NmbPth = getnumberofcores()
      print*, 'lib = ', ParIdx, 'threads = ', NmbPth

      VerTyp = newtype(ParIdx, NmbVer)
      if(VerTyp.eq.0) STOP ' cannot define vertex type'

      TetTyp = newtype(ParIdx, NmbTet)
      if(TetTyp.eq.0) STOP ' cannot define tetrahedron type'

      print*, 'VerTyp = ', VerTyp, 'TetTyp = ', TetTyp

      do i=1,NmbVer
          VerDeg(i) = 0
          VerTem(i) = i
      end do

      ret = begindependency(ParIdx, TetTyp, VerTyp)
      do i=1,NmbTet
          do j=1,4
              ret = adddependency(ParIdx, i, TetVer(j,i))
              VerDeg(TetVer(j,i)) = VerDeg(TetVer(j,i)) + 1
          end do
      end do
      ret = enddependency(ParIdx, sta)
      print*, 'mean = ', sta(1), 'max = ', sta(2)

      t = getwallclock()

      do i=1,NmbItr
          acc = acc + launchparallel(ParIdx,TetTyp,VerTyp,tettemp,3,
     +                               TetVer,TetTem,VerTem)
           acc = acc + launchparallel(ParIdx,TetTyp,VerTyp,vertemp,4,
     +                               TetTem,VerTem,VerDeg,TetVer)
      end do

      print*, 'wall time = ', getwallclock() - t, 'accelerate = ',
     +         acc / (NmbItr*2)

      call stopparallel(ParIdx)

      end


      subroutine tettemp(BegIdx,EndIdx,PthIdx,TetVer,TetTem,VerTem)

      integer i,j,BegIdx,EndIdx,PthIdx,TetVer(4,*)
      real*8 sum,TetTem(*),VerTem(*)

      do i = BegIdx,EndIdx
          sum = 0
          do j=1,4
              sum = sum + VerTem(TetVer(j,i))
          end do
          TetTem(i) = sum
      end do

      return
      end


      subroutine vertemp(BegIdx,EndIdx,PthIdx,TetTem,
     +                   VerTem,VerDeg,TetVer)

      integer i,j,k,BegIdx,EndIdx,PthIdx,TetVer(4,*),VerDeg(*)
      real*8 TetTem(*),VerTem(*)

      do i = BegIdx,EndIdx
          do j=1,4
              k = TetVer(j,i)
              VerTem(k) = VerTem(k) + TetTem(i) / VerDeg(k)
          end do
      end do

      return
      end
