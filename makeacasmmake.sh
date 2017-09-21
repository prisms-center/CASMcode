#Strip all the '/' characters from a string
_smushdir()
{
    local name=$1
    echo "${name///}"
}

#Print each argument passed followed by '\', except the last member
_echo_with_backslash()
{
    local last="${@: -1}"
    for e in $@; do
        echo -ne $e
        #No \ if it's the last file (end of the set)
        if [[ $e != $last ]]; then
            echo "\\"
        fi
    done
}

#Create Makemodule.am file for headers, pass clex, monte_carlo, etc as parameter
_make_for_include()
{
    local target=$1
    #find all subdirectories inside a particular include set, such as monte_carlo
    searchdir=./include/casm/$target
    for d in $(find $searchdir -type d); do
        #Strip the path from all the slashes to get a unique name
        long=$(_smushdir $d)
        #Shorten that name so that it doesn't have "includecasm"
        basename=${long#*includecasm}include

        #If there's no hh files carry on
        headers=$(find $d -type f -name "*.hh")
        if [[ $(echo $headers | wc -w) == 0 ]]; then
            continue
        fi

        #Declare the directory a set of header files has to go in
        echo "${basename}dir=\$(includedir)"/${d#*include\/}
        #Begin listing all the header files that go in the directory
        echo "${basename}_HEADERS=\\"

        _echo_with_backslash $headers
        echo -e
        echo -e

    done
}

_make_for_src()
{
    local target=$1
    #find all files inside a particular src folder
    searchdir=./src/casm/$target

    srcs=$(find $searchdir -type f -name "*.cc")
    echo "libcasm_la_SOURCES +=\\"
    _echo_with_backslash $srcs
    echo -e
}

###############################################################################


#You must run this script from the git root directory
amtype=$1   #src vs inc
target=$2   #clex, monte_carlo, etc

#Deal with the src files here
if [[ "$amtype" == "src" ]]; then
    _make_for_src $target

#Deal with the header files here
elif [[ "$amtype" == "inc" ]]; then
    _make_for_include $target

#Invalid request
else
    echo "Usage:"
    echo "bash makecasmmake.sh amtype target"
    echo "amtype: src or inc"
    echo "target: one of clex, monte_carlo, etc"
fi
