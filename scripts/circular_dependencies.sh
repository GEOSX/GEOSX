#!/bin/bash

# Writes to a dependency_matrix.txt file in the current directory.
# Each row has the component followed by a "yes" or "no" if the component
# has a header/dependency to the other component named at the top of the 
# column.
# Tab delimited, can be exported to Excel (Data --> From Test/CSV)

declare -ar coreComponents=(
 common
 codingUtilities
 dataRepository
 schema
 functions
 constitutive
 mesh
 linearAlgebra
 fieldSpecification
 finiteElement
 finiteVolume
 virtualElement
 discretizationMethods
 fileIO
 physicsSolvers
 events
 mainInterface
)

declare -a depMatrix

length=${#coreComponents[@]}

printf -v header '%s\t' "${coreComponents[@]}"
echo -e " x\t" "${header}" > dependency_matrix.txt


# For each component
for (( i=0; i<$length; i++))
do
  # Name of component, followed by dependency
  declare -a name_dependencies=( ${coreComponents[$i]} )
  echo "Checking dependencies derived from headers for component $name_dependencies..."

  # Check if header appears in component's files
  for (( j=0; j<$length; j++))
  do
    if grep -nr "#include" src/coreComponents/${coreComponents[$i]} | 
       grep -v "//#include" | 
       grep -v "<" | 
       grep -v "unitTests" |
       grep -q "${coreComponents[$j]}"; then

    	name_dependencies+=( yes )
    	depMatrix[($i * $length) + $j]=yes
    else 
    	name_dependencies+=( no )
    	depMatrix[($i * $length) + $j]=no
    fi
  done

  printf -v row '%s\t' "${name_dependencies[@]}"
  echo "${row}" >> dependency_matrix.txt
done

echo -e "\n\nFinish checking dependencies! Dependency matrix written to dependency_matrix.txt"
echo -e "\nNow checking for circular dependencies...\n\n"

# Check for circular dependencies
> circular_dependencies.txt
for (( x=0; x<$length; x++))
do
  # Check half of the matrix
  for (( y=($x + 1); y < $length; y++))
  do
    declare -i a=$x*$length+$y
    declare -i b=$y*$length+$x
  	if [ "${depMatrix[$a]}" == yes ] && [ "${depMatrix[$b]}" == yes ]; then
      echo -e "${coreComponents[$x]} & ${coreComponents[$y]} have a circular dependency" | tee -a circular_dependencies.txt
      
      echo -e "\n${coreComponents[$x]} ---> ${coreComponents[$y]}" >> circular_dependencies.txt
      grep -nr "#include" src/coreComponents/${coreComponents[$x]} | 
       grep -v "//#include" | 
       grep -v "<" | 
       grep -v "unitTests" |
       grep "${coreComponents[$y]}" >> circular_dependencies.txt

      echo -e "\n${coreComponents[$y]} ---> ${coreComponents[$x]}" >> circular_dependencies.txt
      grep -nr "#include" src/coreComponents/${coreComponents[$y]} | 
       grep -v "//#include" | 
       grep -v "<" | 
       grep -v "unitTests" |
       grep "${coreComponents[$x]}" >> circular_dependencies.txt

      echo -e "\n\n" >> circular_dependencies.txt
  	fi
  done
done

echo -e "\nFinished checking for circular dependecies! Any circular dependencies are written to circular_dependencies.txt"