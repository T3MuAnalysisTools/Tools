num=$(ls -1d Set*/ | wc -l)
echo $num
for i in $(seq 1 $num);
do
echo "------------ Set_$i --------------"
cd Set_$i;
nohup ./Set_$i\.sh &> nohup.txt &
#tail nohup.txt 
echo "----------------------------------"
cd ..
done;
