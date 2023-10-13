BEGIN{
    OFS="," 
    if(x == "M"){x = 25}
    if(x == "Y"){x = 24}
    if(x == "X"){x = 23} 
}NR ==1{
    print $0
}NR > 1{
   print $1"0000000000"x, $2,$3
}
