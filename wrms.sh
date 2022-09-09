nline=`wc -l "rms_vel.out" | tr -d '[:alpha:][= =][=.=][=_=]'`
nline=`echo "$nline/2" | bc`
nu=`awk -v nl=$nline 'NR<=nl {x+=$2}END{print x/nl}' rms_vel.out`
echo "Half Urms: $nu"
nu=`awk '{x+=$2;next}END{print x/NR}' rms_vel.out`
echo "Full Urms: $nu"
nu=`awk -v nl=$nline 'NR<=nl {x+=$3}END{print x/nl}' rms_vel.out`
echo "Half Vrms: $nu"
nu=`awk '{x+=$3;next}END{print x/NR}' rms_vel.out`
echo "Full Vrms: $nu"
nu=`awk -v nl=$nline 'NR<=nl {x+=$4}END{print x/nl}' rms_vel.out`
echo "Half Wrms: $nu"
nu=`awk '{x+=$4;next}END{print x/NR}' rms_vel.out`
echo "Full Wrms: $nu"
nu=`awk -v nl=$nline 'NR<=nl {x+=$1}END{print x/nl}' rms_vel_axi.out`
echo "Half Urms k=0: $nu"
nu=`awk '{x+=$1;next}END{print x/NR}' rms_vel_axi.out`
echo "Full Urms k=0: $nu"
nu=`awk -v nl=$nline 'NR<=nl {x+=$2}END{print x/nl}' rms_vel_axi.out`
echo "Half Vrms k=0: $nu"
nu=`awk '{x+=$2;next}END{print x/NR}' rms_vel_axi.out`
echo "Full Vrms k=0: $nu"
nu=`awk -v nl=$nline 'NR<=nl {x+=$3}END{print x/nl}' rms_vel_axi.out`
echo "Half Wrms k=0: $nu"
nu=`awk '{x+=$3;next}END{print x/NR}' rms_vel_axi.out`
echo "Full Wrms k=0: $nu"
nu=`awk -v nl=$nline 'NR<=nl {x+=$1}END{print x/nl}' rms_vel_rnl.out`
echo "Half RNL Streak: $nu"
nu=`awk '{x+=$1;next}END{print x/NR}' rms_vel_rnl.out`
echo "Full RNL Streak: $nu"
nu=`awk -v nl=$nline 'NR<=nl {x+=$2}END{print x/nl}' rms_vel_rnl.out`
echo "Half RNL Roll : $nu"
nu=`awk '{x+=$2;next}END{print x/NR}' rms_vel_rnl.out`
echo "Full RNL Roll: $nu"
