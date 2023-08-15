#!/bin/bash
DIR="$(cd "$(dirname "$0")" >/dev/null 2>&1 ; pwd -P)"

# Función para explicar el funcionamiento
printUsage() {
	echo -e "Usage:"
	printf "i: Mono-metallic NP models in xyz\n"
	printf "n: Number of flakes\n"
	printf "s: Side length for the square graphene flakes (in A)\n"
	printf "l: Layers per flake\n"
	printf "d: distance between the NP and the flakes (d>0) \n"
	printf "o: orientation (Tangent or Perpendicular)\n"
	printf "v: activa modo revisión de estructura\n"
	# Solo uno $flakeNum o $flakeSize puede ser > 0 
}

# Se reciben los argumentos
flakeNum=0; flakeSize=0; orientation="t"; visual=0; ff="G"; reaxCorr=0
while getopts "i:n:s:d:l:o:vf" flag; do
    case $flag in
	i) files="$OPTARG";;
	n) flakeNum="$OPTARG";;
	s) flakeSize="$OPTARG";;
	l) flakeLayer="$OPTARG";;
	d) flakeDist="$OPTARG";;
	o) orientation="$OPTARG";;
	v) visual=1;;
	f) ff="G_r";;
	*) printUsage
		exit 1;;
    esac
done

[ $ff = "G_r" ] && reaxCorr=2   
orientation=${orientation^^}

# Se evalúan los argumentos
isInt='^[0-9]+$'
isFloat='^[0-9]+([.][0-9]+)?$'
if [ $(echo $files | wc -w) -eq 0 ] || [ $flakeLayer -eq 0 ] || [ $flakeDist -eq 0  ] || ! [[ $flakeNum =~ $isInt ]]  || ! [[ $flakeSize =~ $isFloat ]] || ! [[ $flakeLayer =~ $isInt  ]] || ! [[ $flakeDist =~ $isFloat ]] || [[ $flakeNum -gt 0 && $flakeSize -gt 0 ]] || [[ $flakeNum -eq 0 && $flakeSize -eq 0 ]] || [[ "$orientation" != "T" && "$orientation" != "P" ]]  ; then
	echo -e "Malformed arguments."
	printUsage
	exit 1
fi

# Genera el .xyz y .dat
flakeNumBak=$flakeNum; flakeSizeBak=$flakeSize
for file in $files
do
if [ -f $file ]; then
		flakeNum=$flakeNumBak; flakeSize=$flakeSizeBak
		model=$(echo $file | cut -d "." -f 1)
		echo -e "\nCreating the $model@Graphene model. PLEASE WAIT...\n"
		echo -e "FOLDER SUMMARY\nOlguin-Orellana GJ (gabriel.olguin@utalca.cl)\n-----------------------------------------------\n\nMODELING COMMAND:\n   $0 $@" > folderInfo.log
		nohup vmd -dispdev text -e $DIR/flakesGen.tcl -args $model $flakeNum $flakeSize $flakeLayer $flakeDist $DIR $orientation $reaxCorr  &> modelFlakes.log
		echo -e "$(grep "^Para un área superficial" modelFlakes.log)"
		#continúa el procesamiento solo si los flakes fueron creados
		flakeNum=$(ls flake_* 2> /dev/null | wc -w) 
		# se realiza el procedimiento solo si los flakes fueron creados (no se pidieron más de 20)
		if [ $flakeNum -gt 0 ]; then 
			mAtom=$(sed '1q;d' ${model}_c.xyz)
			mElement=$(sed '3q;d' ${model}_c.xyz | awk '{print $1}')
			echo "$DIR/mergexyz -i '$(ls flake_*)' -o mergedFlake.xyz" | bash
			GAtom=$(awk 'NR>2{print $0}' mergedFlake.xyz | wc -l)
			flakeDistNew=$(grep "^     distance" folderInfo.log | awk '{print $NF}')
			[ $ff = "G" ] && newModel=$(echo "${mElement}${mAtom}_n${flakeNum}_s${flakeSize}_l${flakeLayer}_d${flakeDistNew}_o${orientation}")
			[ $ff = "G_r" ] && newModel=$(echo "${mElement}${mAtom}_n${flakeNum}_s${flakeSize}_l${flakeLayer}_d${flakeDistNew}_o${orientation}_Rff")
			vers=$(ls -d $newModel* 2> /dev/null | wc -w)
			newModel=$(echo "${newModel}_v$(echo "$vers+1" | bc)")
			echo "$DIR/mergexyz -i '${model}_c.xyz mergedFlake.xyz' -o $newModel.xyz" | bash
			[ $ff = "G"  ] && xyz2datCuG -i $newModel.xyz -m CuCH 
			[ $ff = "G_r" ] && xyz2datCuG -i $newModel.xyz -m CuCH  -f 
			[ $visual -eq 1 ] && nohup vmd $newModel.xyz &> thrash.log 
			mkdir $newModel
			mv  modelFlakes.log folderInfo.log $model.xyz $newModel.xyz $newModel.dat $newModel
			echo -e "Model $newModel was created."
		fi
		rm -f mainFlake_r00.xyz flake_* $model\_c.xyz 
		mv mergedFlake.xyz  $newModel/mergedFlake${mElement}${mAtom}.xyz; mv mainFlake.xyz $newModel/mainFlake${mElement}${mAtom}.xyz
		rm -f mergedFlake.xyz mainFlake.xyz
	else
		echo -e "\nFile $file doesn't exists." 
	fi
done
rm -f thrash.log 
echo -e "\n$(basename "$0") successfully finished."
