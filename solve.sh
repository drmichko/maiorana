job=$1
mode=20
num=0
file=pi-$job.txt
while read a q rem
do
	tmp=$(( num % mode  ))
	if [ $tmp == $job ] ; then
		echo affinity=$a multiplicity=$q
		./pi.exe $a $q 
	fi
	let num++
done < affinity.txt  > $file 

