DIR="$(cd "$(dirname "$0")" >/dev/null 2>&1 ; pwd -P)"

# INPUT : $DAT lammps file
# OUTPUT: atoms.dat
# 
# Takes dat lammps file and change Cu element number if its coordination is below 9
# 
# Run like:
#            ./getCRD.sh Cu586n14.dat

DAT=$1
dmin=2.8
m=4             # Number of elements in original file
elm="Cu C H Ar"    # The order is important


echo $DAT    > Tar
echo $dmin  >> Tar
echo $m     >> Tar
echo $elm   >> Tar

gfortran -o $DIR/0getCRD $DIR/0getCRD.f90

$DIR/0getCRD
