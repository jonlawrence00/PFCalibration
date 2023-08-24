#!/bin/sh                                                                                                                        
##### It is a final step. This file is generally used for step3 or after ntuplising the sample.
rm *.txt
step=$1
filename=$2
path1=${step}/${filename}

if [ ! "$filename"  ] 
then
    echo "==============================="
    echo "No. of Directories in $step : "
    echo "==============================="
    xrdfs root://se01.indiacms.res.in ls -u /dpm/indiacms.res.in/home/cms/store/user/bkansal/${step}
    echo "==========================================================================================================="
    echo "!!!!!!!!!!! give 2nd argument which is the last required directory in the recent printed path !!!!!!!!!!!!!!!!!!"
    echo "==========================================================================================================="

else
    echo "==============================="
    echo "No. of Directories in $path1 : "
    echo "==============================="
    xrdfs root://se01.indiacms.res.in ls -u /dpm/indiacms.res.in/home/cms/store/user/bkansal/${path1}
    path=$3
    outfile=$4 
    if [ "$path" ]
    then
	dir=$outfile
	mkdir $dir
	for i in 0 1 2 #3 4 5
	do 
            xrdfs root://se01.indiacms.res.in ls -u /dpm/indiacms.res.in/home/cms/store/user/bkansal/${path1}/${path}/000${i}/
	    xrdfs root://se01.indiacms.res.in ls -u /dpm/indiacms.res.in/home/cms/store/user/bkansal/${path1}/${path}/000${i}/ > ${step}_input${i}.txt
#	    sed -i "s#root://se01.indiacms.res.in:1094//dpm/#/dpm/#" ${step}_input${i}.txt
#	    file=$5
	    file=${step}_input${i}
	    for line in $(cat ${file}.txt)
	    do
		echo "$line"
		xrdcp ${line} ${dir}/.
	    done
	done
	hadd -f $outfile.root ${dir}/*
#	mv $outfile.root root_files/.
	##Uncomment the below when you will be ready to delete
	#rm -rf ${dir}
    else
	echo "====================================================================================================================================="
	echo "!!!!!!!!!!! If the last directory in the recent printed path in the format of numeric XXXXX_XXXXX only, then "
	echo " 1. give it as 3rd argument "
        echo " 2. give the output root file name as 4th argument (expected in the format of PGun_RECO_<release name>_<energy range>)"
	echo "====================================================================================================================================="
    fi
fi

