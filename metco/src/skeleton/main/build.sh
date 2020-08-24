echo Construyendo METCO
echo 1- Generando
./gen.sh
echo 2- Configurando
./configure
echo 3- Creando directorios
make dirs
echo 4- Compilando
make