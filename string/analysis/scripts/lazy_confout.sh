max_iter=490

skip=50
for i in $(seq 0 $skip $max_iter); do
    prev=$(($i-$skip))
    
    if [ "$i" -eq "0" ]; then
	sed -i "s/iteration=1/iteration=$i/g" per_iteration_confout.sh
	bash per_iteration_confout.sh
    else
        sed -i "s/iteration=$prev/iteration=$i/g" per_iteration_confout.sh
        bash per_iteration_confout.sh
    fi
done

sed -i "s/iteration=$i/iteration=$max_iter/g" per_iteration_confout.sh
bash per_iteration_confout.sh
#reset back to 1
sed -i "s/iteration=$max_iter/iteration=1/g" per_iteration_confout.sh
bash per_iteration_confout.sh


