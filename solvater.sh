#!/bin/bash
export LC_NUMERIC="en_US.UTF-8" # sirve para encontrar mayor con sort
DIR="$(cd "$(dirname "$0")" >/dev/null 2>&1 ; pwd -P)"

printUsage() {
        echo -e "Usage:
 	i) input file
        a) atmosphere element
	e) distance to the edge (A)
	d) density (kg/m^3)
	m) material (Cu or CH or CuCH)
	"
}

while getopts "i:a:e:d:m:" flag; do
    case $flag in
        i) in="$OPTARG";;
        a) atmos="$OPTARG";;
        e) edge="$OPTARG";;
	d) dens="$OPTARG";;
	m) material="$OPTARG";;
        *) printUsage
                exit 1;;
    esac
done

if [ -z $in ] || [ -z $atmos ] || [ -z $edge ] || [ -z $dens ] || [ -z $material ]; then
	echo -e "Malformed arguments."
        printUsage
        exit 1
fi

#---- Constants definition ----------- 

pi=3.1415926535897932384
av=6.02*10^23

if [ $atmos = "Ar" ]; then
	densg=$(echo "scale=30;$dens/1000" | bc) #g/cm3
	mAtmos=39.948 # g/mol
fi

#if [ $atmos = "C2H2" ]; then

#fi

#------------------------------------

nohup vmd -dispdev text -e $DIR/moveModel.tcl -args ${in%.*} 0 0 0 &> modelFlakes.log

mv ${in%.*}_c.xyz $in

awk '{print sqrt($2*$2), sqrt($3*$3), sqrt($4*$4)}' $in > abs

x=$(sort -k1 -n abs | awk '{print $1}' | tail -n 1)
y=$(sort -k2 -n abs | awk '{print $2}' | tail -n 1)
z=$(sort -k3 -n abs | awk '{print $3}' | tail -n 1)

rm abs modelFlakes.log

if (( $(echo "$x < 1" | bc -l) )); then x=0.001; fi
if (( $(echo "$y < 1" | bc -l) )); then y=0.001; fi
if (( $(echo "$z < 1" | bc -l) )); then z=0.001; fi

#if (( $(echo "$x >= $y" | bc -l) )) && (( $(echo "$x >= $z" | bc -l) )); then
#	npRA=$x
#elif (( $(echo "$y > $x" | bc -l) )) && (( $(echo "$y >= $z" | bc -l) )); then
#	npRA=$y
#else
#	npRA=$z
#fi


#if [ $shape=="s" ]; then
#	# NP volume calculation
#	npVA=$(echo "scale=30;4/3*$pi*$npRA^3" | bc)
#	npRcm=$(echo "scale=30;$npRA*10^(-8)" | bc)
#	npVcm=$(echo "scale=30;4/3*$pi*$npRcm^3" | bc)
#	# box volume (margin $edge A)
#	boxRA=$(echo "scale=30;$npRA+$edge" | bc)
#fi
#if [ $shape=="c" ]; then
	# NP volume calculation	
	npVA=$(echo "scale=30;($x*2)*($y*2)*($z*2)" | bc)
	npLA=$(echo "scale=30;$(echo $npVA | awk '{ print $1^(1/3) }')" | bc)
	npLcm=$(echo "scale=30;$npLA*10^(-8)" | bc)
	npVcm=$(echo "scale=30;$npLcm^3" | bc)
	# box volume (margin $edge A)
	boxRA=$(echo "scale=30;(($x+$y+$z)/3)+$edge" | bc)
#fi
boxLA=$(echo "scale=30;$boxRA*2" | bc)
boxVA=$(echo "scale=30;$boxLA^3" | bc)
boxLcm=$(echo "scale=30;$boxLA*10^(-8)" | bc)
boxVcm=$(echo "scale=30;$boxLcm^3" | bc)

#efective volume (boxVol-npVol) in cm3
effBoxVA=$(echo "scale=30;$boxVA-$npVA" | bc)
effBoxV=$(echo "scale=30;$boxVcm-$npVcm" | bc)

#Se calcula n° de moleculas para la densisdad 298 K y 1 atm (STP)

#gramos de sustancia 
# GivenDens (g) ---> 1 (cm3)
#         x (g) ---> BoxV (cm3)
xDensEmp=$(echo "scale=30;$densg*$boxVcm" | bc)
xDensNP=$(echo "scale=30;$densg*$effBoxV" | bc)

#echo "gramos de sustancia en el volumen de la caja: $xDensEmp"
#echo "gramos de sustancia en el volumen efectivo: $xDensNP"

#Moles
#    mAtmos  (g)  --->    1 (mol)
#    xDens  (g)  --->    x (mol)
xMolEmp=$(echo "scale=30;$xDensEmp/$mAtmos" | bc)
xMolNP=$(echo "scale=30;$xDensNP/$mAtmos" | bc)

#echo "moles de sustancia en el volumen de la caja: $xMolEmp"
#echo "moles de sustancia en el volumen efectivo: $xMolNP"

#Moléculas
#       1 (mol)            --->   6.02E23 (moléculas)
#   moles en la caja (mol) --->   x (moléculas)
xMolecEmp=$(echo "scale=30;$xMolEmp*$av" | bc)
xMolecNP=$(echo "scale=30;$xMolNP*$av" | bc)

#echo "moléculas de sustancia en el volumen de la caja: $xMolecEmp"
#echo "moléculas  de sustancia en el volumen efectivo: $xMolecNP"

npVA=$(printf "%.3f" $npVA)
boxRA=$(printf "%.3f" $boxRA)
boxVA=$(printf "%.3f" $boxVA)
effBoxVA=$(printf "%.3f" $effBoxVA)
xMolecEmp=$(printf "%.0f" $xMolecEmp)
xMolecNP=$(printf "%.0f" $xMolecNP)

# Corroboración
# (at/A3)(mol/at)(g/mol)(A3/mL)
corro=$(echo -e "scale=30;($xMolecNP/$effBoxVA)*(1/($av))*($mAtmos/1)*(1/10^-24)*1000" | bc)
corro=$(printf "%.3f" $corro)

echo -e "
	
COMMAND: $0 $@
	
### $atmos in a empty box

L/2(box)= $boxRA A
v(box)=  $boxVA A^3
molecules inside v(box) = $xMolecEmp

### input model in a $atmos box

v(model)= $npVA  A 
v(effective)= $effBoxVA A^3
molecules inside v(eff) = $xMolecNP

### Reverse calculation (corroboration):
Theoretical Ro for a box of $boxVA A^3 whith $xMolecNP molecules= $corro kg/m^3

Successfully finished 
	" > solvater_i${in%.*}_d${dens}_m${edge}.log

	cat solvater_i${in%.*}_d${dens}_m${edge}.log
	
	out=$(echo $in | cut -f1 -d'.')
	out=$(echo $out$atmos)

#if [ shape ==  "s" ]; 
#then	
#	echo -e "
## like-spherical nanoparticle in an $atmos atmosphere

#tolerance 2.0
#filetype xyz
#output $out.xyz
	
#structure $DIR/simulationFiles/models/$atmos.xyz
#  number $xMolecNP
#  inside box -$boxRA -$boxRA -$boxRA $boxRA $boxRA $boxRA
#  radius 2
#  nloop 20
#end structure

#structure $in
#  number 1
#  radius 5
#  fixed 0. 0. 0. 0. 0. 0. 0.
#  nloop 20
#end structure
#nloop 20
#" >  $out.inp
#fi

#if [ shape ==  "r" ];
#then
	echo -e "
	
## like-rectangular nanoparticle in an $atmos atmosphere

tolerance 2.0
filetype xyz
output $out.xyz
	
structure $DIR/simulationFiles/models/$atmos.xyz
  number $xMolecNP
  inside box -$(echo "scale=3;$x+$edge" | bc) -$(echo "scale=3;$y+$edge" | bc) -$(echo "scale=3;$z+$edge" | bc) $(echo "scale=3;$x+$edge" | bc) $(echo "scale=3;$y+$edge" | bc) $(echo "scale=3;$z+$edge" | bc)
  radius 2
  nloop 20
end structure

structure $in
  number 1
  radius 5
  fixed 0. 0. 0. 0. 0. 0. 0.
  nloop 20
end structure
nloop 20
" >  $out.inp

#fi


packmol < $out.inp
xyz2datCuG -i $out.xyz -m $material
