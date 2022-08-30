count=0

## INFLUX
for a in 0  1200  2000  8508  9000 12000 12800 19304 19600 20000 25200 31940 32120 33600 36400 42200
do

echo $a
echo $count

mkdir ./md/0/$count
mkdir ./md/0/$count/restrained

gmx trjconv -f OutOpen-InOpen.xtc  -s topology/start.gro  -dump $a -o ./md/0/$count/restrained/confout.gro << EOF
0
EOF


(( count++ ))

echo
echo
echo
done

