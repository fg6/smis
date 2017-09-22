

errs=0
MYSMISDIR=`pwd`


sed -i '/export/d' src/makefile 
sed -i '/export MYSMISDIR/d' setup.sh 



setdone=`grep "export MYSMISDIR" setup.sh | wc -l`
if [[ $setdone < 2 ]]; then 
    setsmis=`echo export MYSMISDIR=$MYSMISDIR` 
    sed -i '2i\'"$setsmis" setup.sh
fi

export MYSMISDIR=$MYSMISDIR
setdone=0
setdone=`grep "export MYSMISDIR" src/makefile  | wc -l`
if [[ $setdone < 1 ]]; then
    setsmis=`echo export MYSMISDIR=$MYSMISDIR`
    sed -i '1i\'"$setsmis" src/makefile 
fi





cd $MYSMISDIR/src
if [[ ! -d bamtools ]] || [[ ! -f bamtools/bin/bamtools ]]; then
    rm -rf bamtools
    git clone https://github.com/pezmaster31/bamtools.git
    cd bamtools
    git checkout 71392a251aaea8b797661846b2878e16efd04285   #2.4.0
    mkdir build
    cd build
    cmake ..
    make
fi

cd $MYSMISDIR/src
if [[ ! -f bamtools/bin/bamtools ]]; then
    echo "  !! Error: bamtools not installed properly!"; 
    echo " is your CMake version >= 2.6.4 ? (Check with: $ cmake --version )" 
    errs=$(($errs+1))
    exit
fi




cd $MYSMISDIR/src
mkdir -p $MYSMISDIR/smissv-bin

updated=1
srcs=( smis_sort  smis_rename smis_shred  smis_select_chr spinner_pg )
for code in "${srcs[@]}"; do 
    if [[ ! -f $code ]] || [[ ! -f $MYSMISDIR/smissv-bin/$code ]] || [[ $code -ot $code.c ]] || [[ fasta.h -ot $code.c ]] || [[ smissv.h -ot $code.c ]]; then # || [[ $MYSMISDIR/smissv-bin/$code -ot $code ]] ; then
	updated=0
	rm -f $code
    fi
done

if [[ $updated == 0 ]]; then
    rm -f $MYSMISDIR/smissv-bin/*
    make
    cp smis_sort  smis_rename smis_shred  smis_select_chr  $MYSMISDIR/smissv-bin/.
    rm *.o
    cd ./spinner
    make
    cp spinner_pg ../
    cp spinner_pg $MYSMISDIR/smissv-bin/.
    rm *.o
fi
 
for code in "${srcs[@]}"; do 
    if [[ ! -f $MYSMISDIR/smissv-bin/$code ]]; then
	echo "  !! Error: " $code " not installed properly!"; 
	errs=$(($errs+1))
	exit
    fi
done


if [  $errs -gt 0 ]; then echo " ****  Errors occurred during installation! **** "; echo; exit; 
else echo " Congrats: installation successful!"; fi
