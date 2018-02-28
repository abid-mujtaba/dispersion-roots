#! /bin/bash
# This script simply adds wrapping lines to the plot tex files to enable standalone compilation and then compiles it

# Copy TeX file to /tmp/ for encapsulation
cp $1 /tmp/

# Prepend lines using sed. Source: https://www.cyberciti.biz/faq/bash-prepend-text-lines-to-file/
sed -i '1s;^;\\\documentclass{standalone}\n\\usepackage{tikz}\n\\\begin{document}\n;' /tmp/$1

# Append line at end to close document
sed -i '$ a \\\end{document}' /tmp/$1

# Compile the file
pdflatex /tmp/$1

# Construct pdf file name from tex file name and show it
pdf=$(echo $1 | sed 's/\..*$/\.pdf/')
# vupdf $pdf
