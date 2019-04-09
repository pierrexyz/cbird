simname=ChallengeQuarter
marg=0 
kmaxbisp=0

cd /exports/pierre/cbird

nodenum=1

for stoch in 0 1 2 3 4
do
for number in A
do
for kmax in 0.2 0.25
do
qsub -l nodes=node0$nodenum:ppn=24 -N b$number.$kmax.$simname.$marg.$stoch -o log/output_b$number.$kmax.$simname.stoch$stoch.bisp$kmaxbisp.log -e log/error_b$number.$kmax.$simname.stoch$stoch.bisp$kmaxbisp.log -v aa=$number,b=$kmax,c=$simname,d=$marg,e=$stoch,f=$kmaxbisp mcbird-mercury.sh 
done
done
nodenum=$((nodenum+1))
done