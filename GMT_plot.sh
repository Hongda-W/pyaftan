#!/bin/bash

Initialize() {
    los=-0.1
    gmtset FONT_TITLE 15
    gmtset MAP_TITLE_OFFSET -0.1
    gmtset FONT_LABEL 12
    gmtset MAP_LABEL_OFFSET 0.1

    FTANexe=/projects/howa1663/Code/FTAN_Ye/aftani_c_pgl_amp
    
    sacf=$1

    psout=Sac_FTAN_`echo ${sacf} | awk -F'.SAC' '{print $1}'`.ps
    rm -f $psout
}

FilterSac(){
    sac <<EOF
    read $1
    lp c 0.5 n 4 p 2
    write ${1}.filtered
    quit
EOF
}

PlotSac() {
    if [ -e ${1}.filtered ]; then
        sacfp=${1}.filtered
    else
        sacfp=$1
    fi
    Xoffset=$2
    Yoffset=$3
    if [ ! -e $sacfp ]; then echo $sacfp" not found!" ; exit; fi
    t_info=`saclst KEVNM KSTNM DIST USER0 f $sacfp | awk '{print $2" - "$3"   Distance="$4"km   TimeLength="$5"days"}'`
    #gmtset HEADER_OFFSET -1.2
    #gmtset LABEL_OFFSET -0.3
    xmrk=300
    xtic=100
    ymrk=10
    ytic=0.05
    tmax=400
    ### filter to 5 - 50 sec ###
    amax=`saclst DEPMIN DEPMAX f $sacfp | awk '{if(-$2>$3){print -$2*1.1}else{print $3*1.1}}'`
    REG=-R-${tmax}/${tmax}/-${amax}/${amax}
    SCA=-JX4.1i/1.i
    #gmtset LABEL_OFFSET -1.
    pssac $REG $SCA -Ba${xmrk}f${xtic}:"":/:$label:wESn -O -K -X${Xoffset} -Y${Yoffset} $sacfp >> $psout
    #gmtset LABEL_OFFSET $los
    ### label the sac info ###
    t_info=`saclst KEVNM KSTNM DIST USER0 f $sacfp | awk '{print $2" - "$3"   Distance="$4"   TimeLength="$5"days"}'`
    echo 0.5 11.7 $t_info | pstext -R0/10/0/10 -F+jLT -TO -Wthin,black -J -O -N -P -K >> $psout
    rm $sacfp
}

RunFTAN() {
    sacf=$1
    param=$2
    DispGrv="I_dont_believe_this_file_exists"
    DispPhv="./OBS_phvel.dat"
    echo $param" "$sacf" 1" > param_R.dat
    $FTANexe param_R.dat $DispPhv $DispGrv $3
    rm -f param_R.dat
}

PlotFTAN() {
    ### plot FTAN diagram
    sacf=$1
    Xoffset=$2
    Yoffset=$3
    SCAx=4.1i
    SCAy=3.1i
    ampf1=${sacf}_1_AMP
    disf2=${sacf}_2_DISP.1
    dirpair=.
    perl=3.
    perh=40.
    REG=`echo $perl $perh | awk '{print "-R"log($1)"/"log($2)"/0.2/4.5"}'`
    REG_L=`echo $perl $perh | awk '{print "-R"log($1*0.8)"/"log($2*1.2)"/0.4/5.4"}'`

    if [ ! -e ${dirpair}/${sacf}_amp_tmp ]; then
        awk '{print log($1), $2, $3}' $ampf1 > ${dirpair}/${sacf}_amp_tmp
    fi
    if [ ! -e ${dirpair}/${sacf}_amp.grd ]; then
        xyz2grd ${dirpair}/${sacf}_amp_tmp -G${dirpair}/amp.grd -I0.01 $REG
        grdfilter ${dirpair}/amp.grd -D0 -Fg0.5 -G${dirpair}/${sacf}_amp.grd $REG
        rm -f ${dirpair}/amp.grd
    fi
    grdimage $REG -JX${SCAx}/${SCAy} ${dirpair}/${sacf}_amp.grd \
        -C${dirpair}'/FTAN.cpt' -X${Xoffset} -Y${Yoffset} -K > $psout
    #grdimage $REG -JX${SCAx}/${SCAy} ${dirpair}/amp.grd \
    #    -C${dirpair}'/FTAN.cpt' -X${Xoffset} -Y${Yoffset} -K > $psout
    REG=-R${perl}/${perh}/0.2/4.5
    psbasemap $REG -JX${SCAx}l/${SCAy} -Ba3f3:"period (sec)":/a1.0f0.5:"group vel (km/s)"::."":WSen -V -O -K >> $psout
    awk '$3>3. && $3<34.5{print $3,$4}' $disf2 | psxy -R -J -Sc.15 -W0/0/0 -Gwhite -A -O -K >> $psout
    awk '$3>3. && $3<34.5{print $3,$5}' $disf2 | psxy -R -J -Sc.15 -W0/0/0 -Gpurple -A -O -K >> $psout
}

PlotSNR() {
    sacf=$1
    Xoffset=$2
    Yoffset=$3
    #snrf=${sacf}_amp_snr
    snrf=${sacf}_2_amp_snr
    SCA=-JX2.6il/4.1i
    perl=3.
    perh=40.
    SNR_Max=`awk '{print $3"\n"$5}' $snrf | sort -n | tail -1`
    YMax=`echo "($SNR_Max / 10 + 1) * 10" | bc`
    if (( $YMax > 100 )); then
         YMax=`echo "($YMax / 50 + 1) * 50" | bc`
    fi
    YTick=`echo $YMax / 10 | bc`
    REG=`echo $perl $perh $YMax | awk '{print "-R"$1"/"$2"/0./"$3}'`
    psbasemap $REG $SCA -Bxa3f3+l"period (sec)" -Bya${YTick}+lSNR \
                -BWSen -X${Xoffset} -Y${Yoffset} -V -O -K >> $psout
    pslegend -JX2.6i/4.1i -R -Dx.1i/3.5i+w1.i+jBL -F+pblack -t90 -O -K <<- EOF >> $psout
        S 0.05i - 12p - red 0.2i Positive
        S 0.05i - 12p - blue 0.2i Negative
EOF
    awk '{print $1,$3}' $snrf | psxy -R $SCA -Wthick,red -A -O -K >> $psout
    awk '{print $1,$5}' $snrf | psxy -R -J -Wthick,blue -A -O >> $psout

}

sacf=COR_J49A_WDC.SAC
Initialize ${sacf}
RunFTAN ${sacf} '-1 0.15 4.5 0.8 40. 30. 2.0 3. 0.99 12. 15. 1.5 3.0 0.6 3.' 0.08
FilterSac ${sacf}
PlotFTAN ${sacf} 2. 2.
PlotSac ${sacf} 0. 8.8
PlotSNR ${sacf} 12.3 -8.5

