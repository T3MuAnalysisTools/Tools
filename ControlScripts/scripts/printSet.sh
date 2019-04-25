num=$(ls -1d Set*/ | wc -l)
echo $num
for i in $(seq 1 $num);
do
echo "------------ Set_$i --------------"
tail Set_$i/nohup.txt;
echo "----------------------------------"
done;
