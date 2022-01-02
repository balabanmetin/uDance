mysponge () (
    append=false

    while getopts 'a' opt; do
        case $opt in
            a) append=true ;;
            *) echo error; exit 1
        esac
    done
    shift "$(( OPTIND - 1 ))"

    # if no outfile is provided, use stdout as output file
    if [ $# -ge 1 ] && [ -n "$1" ] ; then
      outfile=$1
    else
      outfile=""
    fi

    tmpfile=$(mktemp -t tmp-sponge.XXXXXXXX) &&
    cat >"$tmpfile" &&
    if "$append"; then
        cat "$tmpfile" >>"$outfile"
    else
#        if [ -f "$outfile" ]; then
#            chmod --reference="$outfile" "$tmpfile"
#        fi
        if [ -f "$outfile" ]; then
            mv "$tmpfile" "$outfile"
        elif [ -n "$outfile" ] && [ ! -e "$outfile" ]; then
            cat "$tmpfile" >"$outfile"
        else
            cat "$tmpfile"
        fi
    fi &&
    rm -f "$tmpfile"
)
export -f mysponge