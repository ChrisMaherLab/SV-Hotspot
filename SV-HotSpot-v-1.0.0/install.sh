#!/bin/bash

if [ "$#" -ne 2 ]; then
    echo
    echo "        ./intall.sh -o installation_directory"
    echo
    echo "        Options :"
    echo "        -o: string   [ the directory where SV-HotSpot is going to be installed ]"
    echo
    exit
fi

args=("$@")

dest_dir=${args[1]}

echo
echo "Creating installation folder..."

if [ ! -d $dest_dir ]; then
    mkdir -p $dest_dir
fi

### write the installation directory to a file 
echo
echo "Copying files to the installation folder..."
cp -r src/* $dest_dir
cd $dest_dir

#### make all files executable files 
echo "adding +x mode"
chmod +x *.r
chmod +x *.sh
chmod +x *.pl

#### add the installation path to the PATH variable 
export PATH="$PWD":$PATH
echo $PATH

#echo "Please add $(pwd)/$dest_dir to your PATH. See details at https://github.com/ChrisMaherLab/SV-Hotspot."  
#echo "Please make you you have installed the prerequisites. See details at https://github.com/ChrisMaherLab/SV-Hotspot."

#echo "---------------"
#echo "Prerequisites"
#echo "--------------"
#echo
#echo "Please make sure you have installed the following tools and they are available in the PATH:"
#echo
#echo "  1. BEDTools version 2.25.0"
#echo "  2. R version 3.1.2 or higher"
#echo "  3. Perl"
#echo 

#echo "-------------------"
#echo "Required R packages"
#echo "-------------------"

#echo
#echo "Please make sure you have installed the following R packages:"
#echo

#echo " 1. peakPick"
#echo " 2. ggplot2"
#echo " 3. reshape2"
#echo " 4. grid"
#echo " 5. gridBase"
#echo " 6. gridExtra"
#echo " 7. gtable"
#echo " 8. ggsignif"
#echo " 9. plyr"
#echo 

#echo "For more details, please refer to https://github.com/ChrisMaherLab/SV-Hotspot."
#echo 

#echo "-------------------"
#echo "Test the tool"
#echo "-------------------"
#echo 
#echo "You may run SV-HotSpot by: perl ./"$dest_dir"/src/sv-hotspot.pl --help" 
#echo


