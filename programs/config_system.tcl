for {set i 0} {$i < 4} {incr i} {
    mol new amyloid$i.pdb
    set sel [atomselect top "all"]
    set vector [measure center $sel weight mass]
    set gx [lindex $vector 0]
    set gy [lindex $vector 1]
    set gz [lindex $vector 2]
    set inverse [list [expr -$gx] [expr -$gy] [expr -$gz]]
    $sel moveby $inverse

    set min 0
    set max 360
    set x_dig [expr ($max - $min) * rand() + $min]
    set y_dig [expr ($max - $min) * rand() + $min]
    set z_dig [expr ($max - $min) * rand() + $min]
    $sel move [transaxis x $x_dig deg]
    $sel move [transaxis y $y_dig deg]
    $sel move [transaxis z $z_dig deg]

    set a 40
    if {$i == 0} {
        set x_pal 0
        set y_pal 0
        set z_pal 0
    } elseif {$i == 1} {
        set x_pal $a
        set y_pal 0
        set z_pal 0
    } elseif {$i == 2} {
        set x_pal [expr 0.5 * $a]
        set y_pal [expr 0.866025 * $a]
        set z_pal 0
    } else {
        set x_pal [expr 0.5 * $a]
        set y_pal [expr 0.288675 * $a]
        set z_pal [expr 0.816496 * $a]
    }
    set pallarel [list $x_pal $y_pal $z_pal]
    $sel moveby $pallarel
    $sel writepdb amyloid_transed$i.pdb
}
exit
