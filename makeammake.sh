searchdir="include/casm/external"

for d in $(find "$searchdir" -type d); do
    dirname="${d///}"
    stripdir=${d#include\/}

    echo "${dirname}dir=\$(includedir)/${stripdir}"
    echo -ne "${dirname}_HEADERS="
    for f in $(find $d -maxdepth 1 -type f -not -name "*.txt"); do
        echo "${f}\\"
    done
    echo ""
done
