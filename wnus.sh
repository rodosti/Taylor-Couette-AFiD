nline=`wc -l "torquenu.out" | tr -d '[:alpha:][= =][=.=]'`
nline=`echo "$nline/2" | bc`
nu=`awk -v nl=$nline 'NR<=nl {x+=$2}END{print x/nl}' torquenu.out`
echo "Half Inner torque: $nu"
nu=`awk '{x+=$2;next}END{print x/NR}' torquenu.out`
echo "Full Inner torque: $nu"
nu=`awk -v nl=$nline 'NR<=nl {x+=$3}END{print x/nl}' torquenu.out`
echo "Half Outer torque: $nu"
nu=`awk '{x+=$3;next}END{print x/NR}' torquenu.out`
echo "Full Outer torque: $nu"
nu=`awk '{x+=$2;next}END{print x/NR}' dissipnu.out`
echo "Full Dissp torque: $nu"
