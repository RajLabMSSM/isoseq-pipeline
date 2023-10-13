# remove monoexons (unspliced reads: CIGAR string $6 lacks "N"
# remove reads with fA above threshold
in=$1
out=$2
threshold=$3

awk -v threshold=$threshold '
    $1 ~ /@/{
        print $0
    }
    $1 !~ /@/ && $6 ~ /N/ {
        # freqA is always final tag
        x = split($NF, a, ":");
        if(a[3] <= threshold){
            print($0)
        } 
    }
' $in > $out
