

mkdir -p smissv-bin
rm -f smissv-bin/*

cd src
make
mv smis_sort  smis_rename smis_shred  smis_select_chr ../smissv-bin/.
rm *.o

cd ./spinner
make
mv spinner_pg ../../smissv-bin/.
rm *.o
 
