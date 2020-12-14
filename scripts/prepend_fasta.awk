BEGIN{ 
    if(x == "M"){x = 25}
    if(x == "Y"){x = 24}
    if(x == "X"){x = 23} 
}{ 
    print ">"$name"0000000000"x, $comment
    print $seq
}
