set terminal qt size 1800,900
set lmargin 24
set rmargin 4
set bmargin 5
#set xlabel 'sweeps' font ',24'

#set title 'Histogram for λ=0, μ_0^2=1' font ',20'
set key font ',22'
set xtics font ',17'
set ytics font ',17'
set grid x,y
set xlabel 'T' font ',24' offset 0,-1


array column1[5]
array column2[5]
array column3[5]
array tit[5]

tit[1]='L=8'
tit[2]='L=16'
tit[3]='L=32'
tit[4]='L=64'


do for [i=1:4] {
column1[i]=1+3*(i-1)
column2[i]=2+3*(i-1)
column3[i]=3+3*(i-1)
set style line i pt 1 lw 2
}


set multiplot layout 2,2

        set yrange [-0.1:1.1]
        set title 'Magnetization' font ',24'
        set ylabel '|M|' font ',24' rotate by 0 offset -6,0
        plot for [i=1:4] 'coexistence/mag2.dat' u column1[i]:column2[i]:column3[i] w errorbars linestyle i title tit[i]
        
        unset yrange
        set key top left
        set title 'Energy' font ',24'
        set ylabel 'E' font ',24' rotate by 0 offset -6,0
        plot for [i=1:4] 'coexistence/ene2.dat' u column1[i]:column2[i]:column3[i] w errorbars linestyle i title tit[i]
        
        set title 'Specific Heat' font ',24'
        set ylabel 'C_V' font ',24' rotate by 0 offset -6,0
        plot for [i=1:4] 'coexistence/hea2.dat' u column1[i]:column2[i]:column3[i] w errorbars linestyle i title tit[i]
        
        set logscale y
        set format y "10^{%L}"
        set key bottom right
        set title 'Susceptibility' font ',24'
        set ylabel 'χ_M' font ',24' rotate by 0 offset -6,0
        plot for [i=1:4] 'coexistence/sus2.dat' u column1[i]:column2[i]:column3[i] w errorbars linestyle i title tit[i]

        #unset yrange


pause -1

unset multiplot
