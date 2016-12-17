#!/bin/bash

# shell script to run pyaftan for the cross-corelation results for JdF
# OCT 2016 Hongda

# Sta_dir_Lst=`ls -d /lustre/janus_scratch/howa1663/CC_JdF/Stack/*/`
Sta_dir_Lst=`find /lustre/janus_scratch/howa1663/CC_JdF/Stack -mindepth 1 -maxdepth 1 -type d`
Aim_dir=/lustre/janus_scratch/howa1663/CC_JdF/AFTAN_JdF
run_aftan=/projects/howa1663/Code/pyaftan/run_aftan_4_JdF.py # python script to run aftan
for sta_dir in $Sta_dir_Lst
do
	cd $sta_dir
	sta_1_Nam=`echo $PWD | gawk -F/ '{print $NF}'`
	CC_Lst=`ls COR'_'$sta_1_Nam'_'*'.'SAC`
	mkdir $Aim_dir/$sta_1_Nam
	cd $Aim_dir/$sta_1_Nam
	for CC_files in $CC_Lst
	do
		sta_2_name=`echo $CC_files | gawk -F[_.] '{print $3}'`
		python $run_aftan $sta_dir/$CC_files $sta_1_Nam $sta_2_name 2>>/lustre/janus_scratch/howa1663/CC_JdF/AFTAN_JdF/python.err 1>>/lustre/janus_scratch/howa1663/CC_JdF/AFTAN_JdF/python.out
		if [ $? -eq 0 ]
		then
			echo 'AFTAN finished for '$sta_dir/$CC_files >> /lustre/janus_scratch/howa1663/CC_JdF/AFTAN_JdF/run_aftan.log
		else
			echo 'AFTAN failed for '$sta_dir/$CC_files >> /lustre/janus_scratch/howa1663/CC_JdF/AFTAN_JdF/run_aftan.err
		fi
	done
done
