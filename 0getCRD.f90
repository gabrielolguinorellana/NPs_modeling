  Program REPLICAT
  Implicit real(a-h,o-z)
  Integer CRD(99999), nn(99999)
  Dimension At(99999,3)
  Character*3 el(9)
  Character*40 DAT
  
  open(1,file="Tar")
        read(1,*)DAT
        read(1,*)dmin
        read(1,*)m
        read(1,*)el(1:m)
        el(m+1)="Cu2"
  
  
  open(11,file=DAT)
  open(12,file="atoms.dat")
  open(16,file="atoms.xyz")
  

  
! ---------------------------------- Read DAT
    read(11,*)
    read(11,*)Nat
    read(11,*)
    read(11,*)
    read(11,*)xi,xf
    read(11,*)yi,yf
    read(11,*)zi,zf
    read(11,*)
    read(11,*)
    read(11,*)
                        
                        write(12,*)
                        write(12,*)Nat, "atoms"
                        write(12,*)m+1, "atom types"
                        write(12,*)
                        write(12,*)xi,xf, "xlo xhi"
                        write(12,*)yi,yf, "ylo yhi"
                        write(12,*)zi,zf, "zlo zhi"
                        write(12,*)
                        write(12,*)"Atoms"
                        write(12,*)
                        
    do j=1,Nat
        read(11,*)i, nn(i), At(i,:)
    enddo
    
!-------------------------------------- Get coordination
    CRD=0
    do j=1,Nat-1
    if(nn(j).eq.1)then
            do i=j+1,Nat
            if(nn(i).eq.1)then
                    rx=At(i,1)-At(j,1)
                    ry=At(i,2)-At(j,2)
                    rz=At(i,3)-At(j,3)
                    
                    dist=sqrt(rx**2+ry**2+rz**2)
                    
                    if(dist.le.dmin)then
                        CRD(i)=CRD(i)+1
                        CRD(j)=CRD(j)+1
                    endif
            endif
            enddo
    endif
    enddo
    
    
!-------------------------------------- Print DAT & xyz with coordination
    nlc=0
    write(16,*)Nat; write(16,*)
    
    do i=1,Nat
        k=nn(i)
        if(k.ne.1)then
            write(12,*)i, k,  At(i,:)
            write(16,*)el(k), At(i,:)
        else
                if(CRD(i).gt.8)then
                    write(12,*)i, k,  At(i,:)
                    write(16,*)el(k), At(i,:)
                else
                    write(12,*)i, m+1, At(i,:)
                    write(16,*)"Ku"  , At(i,:)
                    nlc=nlc+1
                endif
        endif
    enddo
    
    print*
    print*,"Number of atoms with low coordination:",nlc
    print*,"Done!"
    print*

            
    END
