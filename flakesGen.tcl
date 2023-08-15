#!/usr/bin/env tclsh
package require nanotube

proc centerMol {molecule} {
	set centerMolecule [measure center $molecule]
	set xDifference [expr [lindex $centerMolecule 0]*-1]
	set yDifference [expr [lindex $centerMolecule 1]*-1]
	set zDifference [expr [lindex $centerMolecule 2]*-1]
	$molecule moveby "$xDifference $yDifference $zDifference"
}

set model [lindex $argv 0]; set flakeNum [lindex $argv 1]; set flakeSize [lindex $argv 2]; set flakeLayer [lindex $argv 3]; set flakeDist [lindex $argv 4]; set DIR [lindex $argv 5]; set orientation [lindex $argv 6]; set reaxCorr [lindex $argv 7] 

#flakeMargin1
set pi 3.14159265359; set circumMargin 3.5; set maxFlake 20; set flakeDistOld $flakeDist

puts $orientation
if {$orientation == "T"} {set flakeMargin 0
} else {set flakeMargin -1}

#### Preparación ####
# Carga el modelo base
mol load xyz $model.xyz
set structure [atomselect top all]
centerMol $structure
$structure writexyz $model\_c.xyz
set structElem [lsort -unique [$structure get name]]

# Genera la circunesfera considerando el eje más prolongado del modelo
set circumDim [measure minmax $structure]
set circumNP [lindex [lsort -real "[lindex $circumDim 1 0] [lindex $circumDim 1 1] [lindex $circumDim 1 2]"] 2]
set circumAreaNP [expr 4*$pi*pow($circumNP,2)]

# R es el tamaño de la circunesfera con que se calculan los flakes (se ajusta si es ReaxFF)
set circumR [expr $circumNP-$reaxCorr]
set shellR [expr $circumR+$circumMargin]
set circumArea [expr 4*$pi*pow($circumR,2)]
set shellArea [expr 4*$pi*pow($shellR,2)]
# Calcula flakeNum o flakeSize según corresponda 
if {$flakeNum>$flakeSize} {set flakeSize [expr sqrt($shellArea/$flakeNum)]; set flakeFlag 1; set sizeAdd 1
} else {set flakeNum [expr round($shellArea/pow($flakeSize,2))]; set flakeFlag 0}
set flakeArea [expr pow($flakeSize,2)]

#### Modelado de los flakes ####
# Se estima el tamaño que deben tener los flakes
while 1 {
	# Se divide por 10 porque Graphene trabaja con nm 
	graphene -lx [expr double($flakeSize)/10] -ly [expr double($flakeSize)/10] -type zigzag -nlayers 1 -b 1 -a 1 -d 1 -i 1 -cc 0.1418 -ma C-C
	set flake [atomselect top all]
	# Se calcula las dimension del flake 
	set xMin [lindex [measure minmax $flake] 0 0]; set xMax [lindex [measure minmax $flake] 1 0]
	set yMin [lindex [measure minmax $flake] 0 1]; set yMax [lindex [measure minmax $flake] 1 1]
	# Se calcula el largo del flake (promedio entre largo y ancho) 
	set flakeSize [expr ($xMax+$yMax)/2]
	set flakeArea [expr pow($flakeSize,2)]
	puts "\n\n\n el flakeSize es $flakeSize y el flakeArea $flakeArea y $flakeFlag" 
	# Si el tamaño o cantidad no es suficiente, se aumenta
	if {[expr $flakeArea*$flakeNum]>$circumArea} {break
	} elseif {$flakeFlag==1} { set flakeSize [expr $flakeSize+$sizeAdd]; set sizeAdd [expr $sizeAdd+1]
	} else {set flakeNum [expr $flakeNum+1]}
	mol delete top
	puts "\n\n\n el flakeSize es $flakeSize" 
}

# Convierte los átomos de borde -x y +x en átomos basura
set atomSel [atomselect top "x<[expr $xMin+0.5] or x>[expr $xMax-0.5]"]
$atomSel set name Tr; $atomSel set type Tr; $atomSel set element Tr
# Convierte átomos de esquina inferior basura
set atomSel [atomselect top "(x<[expr $xMin+1] and y<[expr $yMin+1]) or (x>[expr $xMax-1] and y<[expr $yMin+1])"]
$atomSel set name Tr; $atomSel set type Tr; $atomSel set element Tr
# Convierte lo átomos de esquina superior en basura
set atomSel [atomselect top "(x<[expr $xMin+3] and y>[expr $yMax-2]) or (x>[expr $xMax-3] and y>[expr $yMax-2])"]
$atomSel set name Tr; $atomSel set type Tr; $atomSel set element Tr

# Convierte los átomos de los nuevo borde en átomos de H
set atomSel [atomselect top "((x>[expr $xMin+0.5] and x<[expr $xMin+1]) or (x<[expr $xMax-0.5] and x>[expr $xMax-1]) or (y>[expr $yMin-1] and y<[expr $yMin+1]) or (y<[expr $yMax+1] and y>[expr $yMax-1])) and not ((x<[expr $xMin+1] and y<[expr $yMin+1]) or (x>[expr $xMax-1] and y<[expr $yMin+1])) and name C"]
$atomSel set name H; $atomSel set type H; $atomSel set element H
# Convierte los átomos no regulares de las esquinas superiores en H
set atomSel [atomselect top "(x>[expr $xMin+1] and x<[expr $xMax-1] and y>[expr $yMin+1] and y<[expr $yMax-1]) and name Tr"]
$atomSel set name H; $atomSel set type H; $atomSel set element H

set flake [atomselect top "name C H"]
centerMol $flake
$flake writexyz flake_0.xyz

#### Agrega las capas adicionaes
if {$flakeLayer>1} {
	mol new flake_0.xyz
	set tempFlake [atomselect top all]
	set xDesp 0.709; set yDesp 1.228; set zDesp 3.349
	for {set n 1} {$n<$flakeLayer} {incr n} {
		if {[expr $n+1]%2==0} {
			$tempFlake moveby "$xDesp $yDesp $zDesp"
			$tempFlake writexyz flake_$n.xyz
		} else { 
			$tempFlake moveby "[expr $xDesp *-1] [expr $yDesp*-1] $zDesp"
			$tempFlake writexyz flake_$n.xyz	
		}
	}	
}

set tempFlake [glob flake_*]
exec $DIR/mergexyz -i "$tempFlake" -o "mainFlake.xyz"
file delete {*}[glob flake_*]
mol delete all

while 1 {
	### Genera el flake centrado (r 0 0)
	mol load xyz mainFlake.xyz
	set flake [atomselect top all]
	set atomNum [$flake num]
	centerMol $flake
	set dimFlake [measure minmax $flake]
	if {$orientation=="T"} {$flake move [trans center "0 0 0" axis y 90]}
	set breakFlag 0
	# r = diametroCircuesfera + alturaFlake/2 + distanciaUsuario (o corregida) 
	set r [expr $circumR+(([lindex $dimFlake 1 2]-([lindex $dimFlake 0 2]))/2)+$flakeDist] 
	$flake moveby "$r 0 0"
	$flake writexyz mainFlake_r00.xyz
	mol delete top

	#### Genera todos los flakes de grafeno. Esta operación se realiza cuando $flakeNum<$maxFlake (lineas en el archivo de vectores)
	set vectorsLine [exec sed "${flakeNum}q;d" $DIR/simulationFiles/polyhedra/vectors]
	set vectors [split $vectorsLine ","]
	for {set n 0} {[lindex $vectors $n] ne ""} {incr n} {
		set v [lindex $vectors $n]
		set M [transvec $v]
		mol load xyz mainFlake_r00.xyz
		set flake [atomselect top all] 	
		$flake move $M
		$flake writexyz flake_$n.xyz
		# Aquí se puede devolver mainFlake a su posición anterior con M($flake)=M^-1($flake') en vez de borrar y cargar	
		mol delete top
	}
	
	# Se evalúa si hay superposiciones. Si las hay, se aumenta la distancia del usuario, se quiebra el ciclo y se prueba nuevamente.
	for {set n 0} {[lindex $vectors $n] ne ""} {incr n} {
		mol load xyz flake_$n.xyz
	 	set flake [atomselect top all]
		set flakeSize [measure minmax $flake]
		for {set m [expr $n+1]} {[lindex $vectors $m] ne ""} {incr m} {
			mol load xyz flake_$m.xyz
			set evalFlake [atomselect top "not x>[expr [lindex $flakeSize 0 0]-$flakeMargin] or not y>[expr [lindex $flakeSize 0 1]-$flakeMargin] or not z>[expr [lindex $flakeSize 0 2]-$flakeMargin] or not x<[expr [lindex $flakeSize 1 0]+$flakeMargin] or not y<[expr [lindex $flakeSize 1 1]+$flakeMargin] or not z<[expr [lindex $flakeSize 1 2]+$flakeMargin]"]
			set newAtomNum [$evalFlake num]
			# Al parecer, se puede simplificar esta operación cambiando la condición a "=="
			if {$newAtomNum<$atomNum} {
				set flakeDist [expr $flakeDist+1]
				file delete {*}[glob flake_*]
				mol delete top
				set breakFlag 1
				break
			}
			mol delete top
		}
		mol delete top	
		if {$newAtomNum < $atomNum} {break}
	}
	if {$breakFlag==0} {break}
	set breakFlag 0
}

# si los flakes solicitados no superan el máximo permitido
if {$flakeNum<=$maxFlake} {
	if {$flakeDist==$flakeDistOld} { puts "\n\nPara un área superficial de [expr round($circumArea)] A² en $model, se usaron $flakeNum flakes de grafeno de [expr round($flakeArea)] A² que cubrieron una superficie total de [expr round($flakeNum*$flakeArea)] A². Los flakes se posicionaron a $flakeDist A de $model.\n$flakeDist"
	} else { puts "\n\nPara un área superficial de [expr round($circumArea)] A² en $model, se usaron $flakeNum flakes de grafeno de [expr round($flakeArea)] A² que cubrieron una superficie total de [expr round($flakeNum*$flakeArea)] A². Los flakes se posicionaron a $flakeDist A de $model en lugar de $flakeDistOld A para evitar superposiciones.\n$flakeDist" }
} else { puts "\n\nPara un área superficial de [expr round($circumArea)] A² en $model, se requieren $flakeNum flakes de grafeno de [expr round($flakeArea)] A² que cubrirían una superficie total de [expr round($flakeNum*$flakeArea)] A²; sin embargo, la versión actual del script solo soporta un máximo de 20 flakes. El modelo no fue creado.\n" }


# Se escribe el resumen completo en el archivo
set fp [open "folderInfo.log" a]
puts $fp "\n   About $model"
puts $fp "     element: $structElem"
puts $fp "     atoms: [$structure num]"
puts $fp "     diameter of the circumsphere (A): [format "%.2f" $circumNP]"
puts $fp "     surface of the circumsphere (A²): [format "%.2f" $circumNP]"
puts $fp "\n   About the flakes" 
puts $fp "     number of flakes: $flakeNum"
puts $fp "     atoms per flake: $atomNum"
puts $fp "     surface of the flakes (A²): [format "%.2f" $flakeArea]"
puts $fp "     layers: $flakeLayer"
puts $fp "     distance $model-flake (A): [expr round ($flakeDist)]"
puts $fp "     orientation (Tangent/Perpendicular): $orientation"
puts $fp "     total atoms: [expr $atomNum*$flakeNum]"
puts $fp "     total covered surface (A²): [format "%.2f" [expr $flakeArea*$flakeNum]]\n"
close $fp

exit 1
